#!/usr/bin/env python

import os
import sys
import json
import numpy as np
from xml.etree import cElementTree as ET
from urllib.error import HTTPError

import pandas as pd
from Bio import Entrez

from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich import box

console = Console()

# Allowed read types
ALLOWED_READ_TYPES = {"paired", "single", "pacbio", "nanopore"}


# ----------------------------------------------------------------------
# 1. Local file validation
# ----------------------------------------------------------------------
def check_local_files(samples, verbose=False):
    found_local = 0
    notfound_acc = []
    err_acc = []

    for sample, items in samples.items():
        R1 = items["R1"]
        R2 = items.get("R2", "")
        if os.path.exists(R1) and (not R2 or os.path.exists(R2)):
            found_local += 1
            if verbose:
                console.log(
                    f"[green]✓[/] {sample}: found locally (R1={R1}, R2={R2 or 'N/A'})"
                )
            continue
        if items.get("accession", ""):
            notfound_acc.append(items["accession"])
            if verbose:
                console.log(
                    f"[yellow]?[/] {sample}: not local, will try SRA accession {items['accession']}"
                )
        else:
            err_acc.append(sample)
            if verbose:
                console.log(
                    f"[red]✗[/] {sample}: no local file and no accession provided"
                )

    return found_local, notfound_acc, err_acc


# ----------------------------------------------------------------------
# 2. SRA / remote validation
# ----------------------------------------------------------------------
def fetch_sra_info(accessions, verbose=False):
    found_accessions = []
    sizes_gb = []

    if verbose:
        console.log(f"[cyan]Fetching SRA info for {len(accessions)} accessions...[/]")

    for i in range(0, len(accessions), 500):
        batch = accessions[i : i + 500]
        if verbose:
            console.log(f"  Querying batch {i+1}–{min(i+500, len(accessions))}")

        try:
            handle = Entrez.efetch(
                db="sra", id=batch, retmax=1000, rettype="full", retmode="xml"
            )
        except HTTPError as err:
            console.print(
                Panel.fit(
                    f"HTTP Error {err.code} during Entrez.efetch() on batch {batch[:3]}...\n"
                    "Check your NCBI API key or try again later.",
                    title="SRA Fetch Error",
                    border_style="red",
                )
            )
            sys.exit(1)

        record = handle.read()
        if not isinstance(record, str):
            record = record.decode("utf-8")
        root = ET.fromstring(record)

        for run in root.findall(".//RUN"):
            acc = run.attrib.get("accession")
            sra_file = run.find(".//SRAFile")
            if sra_file is not None:
                size = sra_file.attrib.get("size")
                if size:
                    size_gb = round(int(size) / 1024**3, 2)
                    found_accessions.append(acc)
                    sizes_gb.append(size_gb)
                    if verbose:
                        console.log(f"    Found {acc} ({size_gb} GB)")

    return found_accessions, sizes_gb


def validate_remote_samples(notfound_acc, verbose=False):
    if not notfound_acc:
        return 0.0

    if verbose:
        console.log(f"Performing remote SRA search on {len(notfound_acc)} samples...")
    else:
        console.print(f"Performing remote SRA search on {len(notfound_acc)} samples...")

    found_acc, sizes = fetch_sra_info(notfound_acc, verbose)
    missing = set(notfound_acc) - set(found_acc)

    if missing:
        console.print(
            Panel.fit(
                f"Accessions not found:\n{missing}\n\n"
                "Check that these are valid SRA run accessions.",
                title="SRA Accession Error",
                border_style="red",
            )
        )
        sys.exit(1)

    total_gb = round(sum(sizes))
    if verbose:
        console.log(
            f"Found all {len(found_acc)} accessions. Total download: ~{total_gb} GB"
        )
    else:
        console.print(
            f"Found all {len(found_acc)} accessions. Total download: ~{total_gb} GB"
        )
    return total_gb


# ----------------------------------------------------------------------
# 3. Main validation orchestrator
# ----------------------------------------------------------------------
def validate_samples(samples, quiet=False, verbose=False):
    if not quiet:
        console.print(
            Panel.fit(
                "[dim]Validating sample availability (local files or SRA).\n"
                "For local files, provide full paths in R1/R2 columns or place files in datadir with <sample>_{1,2}.fastq.gz[/]",
                title="Sample Validation",
                subtitle="In Progress...",
                border_style="cyan",
            )
        )

    # Display read-type summary
    read_types = {}
    for s, d in samples.items():
        rt = d.get("read_type", "paired")
        read_types[rt] = read_types.get(rt, 0) + 1

    table = Table(title="Sample Read Types", box=box.ROUNDED, title_style="bold cyan")
    table.add_column("Sample ID", style="cyan")
    table.add_column("Read Type", style="green")
    table.add_column("R1", style="white")
    table.add_column("R2", style="white")
    for s, d in samples.items():
        table.add_row(s, d.get("read_type", "paired"), d.get("R1", ""), d.get("R2", ""))
    console.print(table)

    if verbose:
        console.print("\n[cyan]Read type distribution:[/]")
        for rt, count in read_types.items():
            console.print(f"  {rt}: {count} samples")

    # Local check
    found_local, notfound_acc, err_acc = check_local_files(samples, verbose)
    if verbose:
        console.log(f"{found_local} samples pre-downloaded locally.")
    else:
        console.print(f"{found_local} samples pre-downloaded locally.")

    if err_acc:
        console.print(
            Panel.fit(
                f"Samples without local files and no accession:\n{err_acc}\n\n"
                "Please provide either the files or valid SRA accessions.",
                title="Local File Error",
                border_style="red",
            )
        )
        sys.exit(1)

    # Remote check for missing ones
    if notfound_acc:
        total_gb = validate_remote_samples(notfound_acc, verbose)
        if not quiet:
            if verbose:
                console.log(f"Will download ~{total_gb} GB of data.")
            else:
                console.print(f"Will download ~{total_gb} GB of data.")

    if verbose:
        console.log("[green]Validation complete.[/]")


# ----------------------------------------------------------------------
# 4. CSV parsing helpers
# ----------------------------------------------------------------------
def read_sample_csv(filepath, datadir, verbose=False):
    """
    Read sample CSV using pandas.
    - Requires a header row with the exact required columns.
    - No extra columns are allowed.
    - Raises a clear error if extra columns or invalid read_type are present.
    """
    if not os.path.isfile(filepath) or not filepath.endswith(".csv"):
        console.print(
            Panel.fit(
                f"File '{filepath}' does not exist or is not a .csv file.",
                title="CSV File Error",
                border_style="red",
            )
        )
        sys.exit(1)

    try:
        df = pd.read_csv(filepath, comment="#", dtype=str)
    except Exception as e:
        console.print(
            Panel.fit(
                f"Failed to read CSV file '{filepath}':\n{e}",
                title="CSV Read Error",
                border_style="red",
            )
        )
        if verbose:
            console.print_exception()
        sys.exit(1)

    if df.empty:
        console.print(
            Panel.fit(
                f"CSV file '{filepath}' is empty or contains only comments.",
                title="Empty CSV",
                border_style="red",
            )
        )
        sys.exit(1)

    # Strip whitespace from column names
    df.columns = df.columns.str.strip()
    required_cols = ["sample_id", "accession", "assembly", "R1", "R2", "read_type"]
    required_lower = [c.lower() for c in required_cols]

    # Check for missing required columns
    missing = [c for c in required_lower if c not in df.columns.str.lower()]
    if missing:
        console.print(
            Panel.fit(
                f"Missing required columns: {missing}\n"
                f"Required columns are: {required_cols}",
                title="Missing Columns",
                border_style="red",
            )
        )
        sys.exit(1)

    # Check for extra columns
    extra = [c for c in df.columns if c.lower() not in required_lower]
    if extra:
        console.print(
            Panel.fit(
                f"Extra columns found: {extra}\n"
                f"Only these columns are accepted: {required_cols}",
                title="Extra Columns",
                border_style="red",
            )
        )
        sys.exit(1)

    # Normalise column names to exact case
    df.columns = required_cols

    # Normalise data
    df["read_type"] = df["read_type"].fillna("paired")
    df["read_type"] = df["read_type"].astype(str).str.lower().str.strip()
    df["read_type"] = df["read_type"].replace("", "paired")

    # --- Validate read_type ---
    invalid = df[~df["read_type"].isin(ALLOWED_READ_TYPES)]
    if not invalid.empty:
        console.print(
            Panel.fit(
                f"Invalid read_type(s) found:\n{invalid[['read_type']].drop_duplicates().to_string(index=False)}\n"
                f"Allowed values: {sorted(ALLOWED_READ_TYPES)}",
                title="Invalid Read Type",
                border_style="red",
            )
        )
        sys.exit(1)

    df["sample_id"] = df["sample_id"].fillna(df["accession"])
    df["sample_id"] = df["sample_id"].astype(str).str.strip()
    df["sample_id"] = df["sample_id"].replace("", np.nan)
    df = df.dropna(subset=["sample_id"])

    df["assembly"] = df["assembly"].fillna(df["sample_id"])

    if not datadir.endswith(os.sep):
        datadir = datadir + os.sep

    df["R1"] = df["R1"].fillna(datadir + df["sample_id"] + "_1.fastq.gz")
    df["R2"] = df["R2"].fillna(
        df.apply(
            lambda row: (
                datadir + row["sample_id"] + "_2.fastq.gz"
                if row["read_type"] == "paired"
                else ""
            ),
            axis=1,
        )
    )

    df["accession"] = df["accession"].fillna("")

    df.set_index("sample_id", inplace=True, verify_integrity=True)

    return df


def check_duplicates(df):
    if df.index.duplicated().any():
        console.print(
            Panel.fit(
                f"Duplicate sample_ids found: {df.index[df.index.duplicated()].tolist()}",
                title="Duplicate Sample IDs",
                border_style="red",
            )
        )
        sys.exit(1)

    if df.duplicated().any():
        console.print(
            Panel.fit(
                "Duplicate rows found. Please check your sample list.",
                title="Duplicate Rows",
                border_style="red",
            )
        )
        sys.exit(1)


def build_assemblies_dict(df, verbose=False):
    assemblies = {}
    for assembly in df["assembly"].unique():
        sub = df[df["assembly"] == assembly]
        types = sub["read_type"].unique()
        if len(types) > 1:
            console.print(
                Panel.fit(
                    f"Assembly '{assembly}' has mixed read_types: {types}. "
                    "All samples in the same assembly must have the same read_type.",
                    title="Assembly Conflict",
                    border_style="red",
                )
            )
            sys.exit(1)

        assemblies[assembly] = {
            "R1": sub["R1"].tolist(),
            "R2": sub["R2"].tolist(),
            "sample_id": sub.index.tolist(),
            "accession": sub["accession"].tolist(),
            "read_type": types[0] if len(types) == 1 else "paired",
        }
        if verbose:
            console.log(
                f"Assembly '{assembly}' -> {len(sub)} samples, type={types[0] if types else 'paired'}"
            )
    return assemblies


def build_samples_dict(df, verbose=False):
    samples = {}
    for sample_id, row in df.iterrows():
        samples[sample_id] = {
            "R1": row["R1"],
            "R2": row["R2"],
            "accession": row["accession"],
            "assembly": row["assembly"],
            "read_type": row["read_type"],
        }
        if verbose:
            console.log(
                f"Sample '{sample_id}': read_type={row['read_type']}, R1={row['R1']}, R2={row['R2']}"
            )
    return samples


def write_json_files(samples, assemblies, outdir, nowstr):
    logdir = os.path.join(outdir, ".vomix", "log", f"vomix{nowstr}")
    try:
        os.makedirs(logdir, exist_ok=True)
    except Exception as e:
        console.print(
            Panel.fit(
                f"Could not create log directory '{logdir}':\n{e}",
                title="Directory Creation Error",
                border_style="red",
            )
        )
        sys.exit(1)

    sample_path = os.path.join(logdir, "samples.json")
    assembly_path = os.path.join(logdir, "assemblies.json")

    try:
        with open(sample_path, "w") as f:
            json.dump(samples, f, indent=2)
        with open(assembly_path, "w") as f:
            json.dump(assemblies, f, indent=2)
    except Exception as e:
        console.print(
            Panel.fit(
                f"Could not write JSON files to '{logdir}':\n{e}",
                title="JSON Write Error",
                border_style="red",
            )
        )
        sys.exit(1)

    return sample_path, assembly_path


def check_existing_jsons(samples, outdir):
    samplejson = os.path.join(outdir, ".vomix", "samples.json")
    if not os.path.exists(samplejson):
        return False

    try:
        with open(samplejson, "r") as f:
            old = json.load(f)
    except Exception:
        return False

    return old == samples


# ----------------------------------------------------------------------
# 5. Main parser orchestrator
# ----------------------------------------------------------------------
def parse_sample_list(
    f, datadir, outdir, email, api_key, nowstr, quiet=False, verbose=False
):
    try:
        df = read_sample_csv(f, datadir, verbose)
        check_duplicates(df)

        if verbose:
            console.log("[cyan]DataFrame after cleaning:[/]")
            console.log(df.to_string())

        samples = build_samples_dict(df, verbose)
        assemblies = build_assemblies_dict(df, verbose)

        write_json_files(samples, assemblies, outdir, nowstr)

        if check_existing_jsons(samples, outdir):
            if not quiet:
                if verbose:
                    console.log(
                        f"[yellow]Skipping validation: samples.json unchanged from {f}.[/]\n"
                        f"To redo validation, delete [cyan]{os.path.join(outdir, '.vomix', 'samples.json')}[/] and re-run."
                    )
                else:
                    console.print(
                        Panel.fit(
                            f"Skipping validation: samples.json unchanged from {f}.\n"
                            f"To redo validation, delete {os.path.join(outdir, '.vomix', 'samples.json')} and re-run.",
                            title="Skipping Validation",
                            border_style="yellow",
                        )
                    )
            return samples, assemblies

        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key

        validate_samples(samples, quiet=quiet, verbose=verbose)

        final_sample_json = os.path.join(outdir, ".vomix", "samples.json")
        final_assembly_json = os.path.join(outdir, ".vomix", "assemblies.json")

        with open(final_sample_json, "w") as f:
            json.dump(samples, f, indent=2)
        with open(final_assembly_json, "w") as f:
            json.dump(assemblies, f, indent=2)

        if verbose:
            console.log(
                f"[green]✓[/] Final samples.json written to {final_sample_json}"
            )
            console.log(
                f"[green]✓[/] Final assemblies.json written to {final_assembly_json}"
            )
        else:
            console.print(
                Panel.fit(
                    f"[dim]Sample list parsed successfully.\n"
                    f"samples.json written to {final_sample_json}\n"
                    f"assemblies.json written to {final_assembly_json}[/]",
                    title="Success",
                    border_style="green",
                )
            )

        return samples, assemblies

    except Exception as e:
        console.print(
            Panel.fit(
                f"Unexpected error during sample list parsing:\n{type(e).__name__}: {e}",
                title="Fatal Error",
                border_style="red",
            )
        )
        if verbose:
            console.print_exception()
        sys.exit(1)
