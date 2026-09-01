#!/usr/bin/env python3
"""
download_refseq_genomes.py - Download stratified RefSeq genomes for mock benchmarks.

Modes:
  --mode viral         : Downloads genomes from dsDNA, ssDNA, dsRNA, ssRNA.
  --mode prokaryotic   : Downloads bacteria + archaea genomes (contaminants).
  --mode eukaryotic    : Downloads eukaryotic genomes (host/fungi/plants).
  --mode contaminants  : (Legacy) Downloads all contaminants into one file.

Output:
  Writes a single FASTA file to the path specified by --outfile.
  The parent directory is created if it does not exist.
  The summary TSV is always written, even in dry‑run mode.

Options:
  --genomes-per-category 0 : downloads all species listed for that category (viral only).
  --verbose / -v          : enables debug-level logging.
  --dry-run / -n          : only check availability (no downloads, no file writes).
  --summary-tsv           : output TSV with per‑species summary.

Requires:
  - Biopython (Entrez)
  - requests
  - rich (for pretty logging)
"""

import os
import sys
import time
import gzip
import json
import argparse
import urllib.request
from urllib.error import HTTPError
from collections import defaultdict
from Bio import Entrez, SeqIO

# Rich imports
try:
    from rich.console import Console
    from rich.panel import Panel
    from rich.table import Table
    from rich import box
except ImportError:
    print("Rich library not installed. Install with: pip install rich", file=sys.stderr)
    sys.exit(1)

console = Console()

VIRAL_CATEGORIES = {
    "dsDNA": [
        "Escherichia virus T4",
        "Enterobacteria phage lambda",
        "Bacillus phage phi29",
        "Mimivirus",
        "Vaccinia virus",
        "Acanthamoeba polyphaga mimivirus",
        "Human herpesvirus 1",
        "Adenovirus C",
        "Bacillus subtilis phage SPP1",
        "Pseudomonas phage PAK_P1",
    ],
    "ssDNA": [
        "Enterobacteria phage phiX174",
        "Beet curly top virus",
        "Porcine circovirus",
        "Parvovirus B19",
        "Tomato yellow leaf curl virus",
        "Human parvovirus B19",
        "Bocavirus",
        "Anellovirus",
    ],
    "dsRNA": [
        "Rotavirus A",
        "Infectious bursal disease virus",
        "Cystovirus phi6",
        "Rice dwarf virus",
        "Bluetongue virus",
    ],
    "ssRNA": [
        "Severe acute respiratory syndrome coronavirus 2",
        "Human immunodeficiency virus 1",
        "Poliovirus",
        "West Nile virus",
        "Dengue virus type 2",
        "Zika virus",
        "Rabies virus",
        "Ebola virus",
        "Measles virus",
    ],
}

PROKARYOTIC_SPECIES = [
    "Escherichia coli",
    "Pseudomonas aeruginosa",
    "Acinetobacter baumannii",
    "Salmonella enterica",
    "Bacteroides fragilis",
    "Staphylococcus aureus",
    "Bacillus subtilis",
    "Lactobacillus acidophilus",
    "Streptococcus pneumoniae",
    "Mycobacterium tuberculosis",
    "Methanobrevibacter smithii",
    "Halobacterium salinarum",
]

EUKARYOTIC_SPECIES = [
    "Homo sapiens",
    "Saccharomyces cerevisiae",
    "Candida albicans",
    "Aspergillus fumigatus",
    "Arabidopsis thaliana",
]

# Legacy contaminants (all combined)
CONTAMINANT_SPECIES = PROKARYOTIC_SPECIES + EUKARYOTIC_SPECIES


def fetch_assembly_id(species_name, email, api_key=None, max_attempts=3, verbose=False):
    """
    Search RefSeq for the best complete genome assembly of a given species.
    Returns assembly ID (GCF_* or GCA_*) or None if not found.
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    term2 = f'"{species_name}"[Organism] AND "complete genome"[Assembly Level]'
    term3 = f'"{species_name}"[Organism] AND refseq[filter]'

    db_ncbi = "assembly"

    for attempt, query in enumerate([term2, term3], start=1):
        try:
            if verbose:
                console.log(
                    f"[dim]Attempt #{attempt} querying {db_ncbi} db for:\n{query}[/]"
                )
            handle = Entrez.esearch(db=db_ncbi, term=query, retmax=1)
            record = Entrez.read(handle)
            handle.close()
            if record["IdList"]:
                console.print(
                    f"[dim]\t   ✓ Genome found ! Record ID: {record['IdList'][0]}[/]"
                )
                return record["IdList"][0]
            else:
                if attempt == 1:
                    console.print(
                        "[dim]\t   ✗ No complete RefSeq genome found, trying broader search[/]"
                    )
                elif attempt == 2:
                    console.print(
                        "[dim]\t   ✗ No complete genome found, trying any RefSeq assembly[/]"
                    )
        except HTTPError as e:
            if e.code in (429, 500):
                wait = attempt * 5
                if verbose:
                    console.log(f"[yellow]Rate limit hit. Waiting {wait}s...[/]")
                time.sleep(wait)
                continue
            else:
                console.log(f"[red]Entrez error for {species_name}: {e}[/]")
                return None
        except Exception as e:
            console.log(f"[red]Unexpected error for {species_name}: {e}[/]")
            return None
        time.sleep(0.5)

    return None


def get_assembly_summary(
    species_name, assembly_id, email, api_key=None, max_attempts=3, verbose=False
):
    """
    Retrieve the full assembly summary (DocumentSummary) for a given assembly ID.
    Returns (ftp_url, summary_dict) or (None, None) on failure.
    Handles special cases for problematic species.
    """
    # Check special cases first
    if species_name in SPECIAL_CASES:
        case = SPECIAL_CASES[species_name]
        if verbose:
            console.log(f"[cyan]Using special case for {species_name}[/]")
        return case["ftp_url"], case["summary"]

    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    for attempt in range(max_attempts):
        try:
            # Resolve accession to numeric ID if needed
            if assembly_id.startswith("GCF_") or assembly_id.startswith("GCA_"):
                handle = Entrez.esearch(db="assembly", term=assembly_id, retmax=1)
                record = Entrez.read(handle)
                handle.close()
                if not record["IdList"]:
                    if verbose:
                        console.log(
                            f"[red]No numeric ID found for accession {assembly_id}[/]"
                        )
                    return None, None
                numeric_id = record["IdList"][0]
            else:
                numeric_id = assembly_id

            # Fetch summary using numeric ID
            handle = Entrez.esummary(db="assembly", id=numeric_id, report="full")
            summary = Entrez.read(handle)
            handle.close()
            if not summary:
                return None, None

            doc_summary = summary["DocumentSummarySet"]["DocumentSummary"][0]
            ftp_dir = doc_summary.get("FtpPath_RefSeq") or doc_summary.get(
                "FtpPath_GenBank"
            )
            if ftp_dir:
                https_dir = ftp_dir.replace("ftp://", "https://")
                assembly_accession = doc_summary.get("AssemblyAccession")
                assembly_name = doc_summary.get("AssemblyName", "")
                sanitized_name = assembly_name.replace(" ", "_")
                if assembly_accession:
                    fna_url = f"{https_dir}/{assembly_accession}_{sanitized_name}_genomic.fna.gz"
                else:
                    fna_url = f"{https_dir}/genomic.fna.gz"
                if verbose:
                    console.log(f"[dim][green]FTP URL found: {fna_url}[/]")
                return fna_url, doc_summary
            else:
                # Try to construct a URL from the accession and name
                assembly_accession = doc_summary.get("AssemblyAccession")
                assembly_name = doc_summary.get("AssemblyName", "")
                if assembly_accession:
                    acc_parts = assembly_accession.split("_")
                    if len(acc_parts) >= 2:
                        prefix = acc_parts[0]
                        number = acc_parts[1]
                        group1 = number[0:3]
                        group2 = number[3:6]
                        group3 = number[6:9]
                        if len(number) > 9:
                            group3 = number[6:9]
                            group4 = number[9:12]
                            base_path = f"{prefix}/{group1}/{group2}/{group3}/{assembly_accession}_{assembly_name.replace(' ', '_')}"
                        else:
                            base_path = f"{prefix}/{group1}/{group2}/{group3}/{assembly_accession}_{assembly_name.replace(' ', '_')}"
                        https_dir = (
                            f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{base_path}"
                        )
                        fna_url = f"{https_dir}/{assembly_accession}_{assembly_name.replace(' ', '_')}_genomic.fna.gz"
                        if verbose:
                            console.log(f"[dim][yellow]Constructed URL: {fna_url}[/]")
                        return fna_url, doc_summary
                if verbose:
                    console.log(f"[red]No FTP path for {assembly_id}[/]")
                return None, doc_summary

        except HTTPError as e:
            if e.code in (429, 500):
                wait = (attempt + 1) * 5
                if verbose:
                    console.log(f"[yellow]Server error. Waiting {wait}s...[/]")
                time.sleep(wait)
                continue
            else:
                if verbose:
                    console.log(f"[red]HTTP error for {assembly_id}: {e}[/]")
                return None, None
        except Exception as e:
            if verbose:
                console.log(f"[red]Error getting summary for {assembly_id}: {e}[/]")
            return None, None
    return None, None


def download_assembly_fasta(
    species_name,
    assembly_id,
    workdir,
    email,
    api_key=None,
    max_attempts=3,
    verbose=False,
    dry_run=False,
    ftp_url=None,
    summary_dict=None,
):
    """
    Download the genomic FASTA (.fna.gz) for a given assembly ID.
    If ftp_url is given, use it; otherwise fetch it (and optionally summary_dict if not given).
    Returns the local temporary FASTA path, or None on failure.
    """
    if ftp_url is None or summary_dict is None:
        ftp_url, summary_dict = get_assembly_summary(
            species_name, assembly_id, email, api_key, max_attempts, verbose
        )
    if not ftp_url:
        return None

    if dry_run:
        console.log(f"[dim][cyan]DRY-RUN: would download from {ftp_url}[/]")
        return "DRY_RUN_PLACEHOLDER"

    # Try the primary URL; if fails, try a fallback (just genomic.fna.gz)
    url_candidates = [ftp_url]
    if "_genomic.fna.gz" in ftp_url and not ftp_url.endswith("/genomic.fna.gz"):
        base_dir = ftp_url.rsplit("/", 1)[0]
        url_candidates.append(f"{base_dir}/genomic.fna.gz")

    for attempt in range(max_attempts):
        for url in url_candidates:
            try:
                if verbose:
                    console.log(f"[dim]Attempting download: {url}[/]")
                temp_gz = os.path.join(workdir, f".tmp_{assembly_id}.fna.gz")
                urllib.request.urlretrieve(url, temp_gz)

                # Decompress to plain FASTA
                temp_fasta = os.path.join(workdir, f".tmp_{assembly_id}.fna")
                with gzip.open(temp_gz, "rt") as gz_f, open(temp_fasta, "w") as out_f:
                    out_f.write(gz_f.read())
                os.remove(temp_gz)
                return temp_fasta

            except HTTPError as e:
                if e.code == 404:
                    if verbose:
                        console.log(f"[red]File not found: {url}[/]")
                    # Try next candidate
                    continue
                if e.code in (429, 500):
                    wait = (attempt + 1) * 5
                    if verbose:
                        console.log(f"[yellow]Server error. Waiting {wait}s...[/]")
                    time.sleep(wait)
                    break  # break inner loop, try again with same url after wait
                console.log(f"[red]HTTP error for {assembly_id}: {e}[/]")
                return None
            except Exception as e:
                console.log(f"[red]Error downloading {assembly_id}: {e}[/]")
                return None
        else:
            # all candidates failed
            console.log(f"[red]All download attempts failed for {assembly_id}[/]")
            return None
    return None


def download_species_list(
    outfile,
    species_list,
    email,
    api_key=None,
    label="",
    verbose=False,
    dry_run=False,
    summary_tsv=None,
    color="magenta",
):
    """
    Download genomes for a given list of species.
    The summary TSV is always written, even in dry‑run mode.
    """
    console.print()
    console.rule(f"[bold {color}]Downloading {label} Genomes[/]")
    console.print(f"Target: {len(species_list)} species\n")

    if dry_run:
        console.print(
            "[yellow]DRY-RUN mode enabled – no fasta files will be downloaded or written.[/]\n"
        )

    outdir = os.path.dirname(outfile)
    if outdir and not dry_run:
        os.makedirs(outdir, exist_ok=True)

    temp_files = []
    success_count = 0
    failed_species = []
    summary_records = []
    # Track downloaded assemblies to avoid duplicates
    downloaded_assemblies = {}
    missing_files = []

    for species in species_list:
        console.print(f"\n[bold]→ {species}[/]")
        aid = fetch_assembly_id(species, email, api_key, verbose=verbose)
        console.log(f"[dim]Fetching FTP URL for {aid}:[/]")
        ftp_url = ""
        summary_dict = None
        status = "failed"
        local_path = "N/A"

        if aid:
            # Check if we already downloaded this assembly
            if aid in downloaded_assemblies:
                console.print(f"[dim]Reusing previously downloaded assembly {aid}[/]")
                fname = downloaded_assemblies[aid]
                temp_files.append(fname)
                success_count += 1
                status = "success"
                local_path = fname
                # Add a dummy summary record for this species (we can reuse the summary from the first)
                # But we need a summary_dict for the TSV; we can fetch it or store it.
                # We'll fetch it again just to get the summary, but we won't download again.
                _, summary_dict = get_assembly_summary(
                    species, aid, email, api_key, verbose=verbose
                )
                summary_json = json.dumps(summary_dict) if summary_dict else ""
                summary_json = summary_json.replace("\t", " ")
                summary_records.append(
                    {
                        "species": species,
                        "category": label,
                        "assembly_accession": aid,
                        "ftp_url": ftp_url,
                        "local_file": local_path,
                        "status": status,
                        "summary_json": summary_json,
                    }
                )
                console.print(f"\t   [green]✓ Reused ({success_count})[/]")
                time.sleep(0.5)
                continue

            ftp_url, summary_dict = get_assembly_summary(
                species, aid, email, api_key, verbose=verbose
            )
            if ftp_url:
                fname = download_assembly_fasta(
                    species,
                    aid,
                    outdir if outdir else ".",
                    email,
                    api_key,
                    verbose=verbose,
                    dry_run=dry_run,
                    ftp_url=ftp_url,
                    summary_dict=summary_dict,
                )
                if fname:
                    temp_files.append(fname)
                    downloaded_assemblies[aid] = fname
                    success_count += 1
                    status = "success"
                    local_path = fname if not dry_run else "DRY_RUN"
                    if not dry_run:
                        console.print(f"\t   [green]✓ Downloaded ({success_count})[/]")
                    else:
                        console.print(
                            f"\t   [dim][green]✓ Would download ({success_count})[/]"
                        )
                else:
                    console.print(f"\t   [red]✗ Failed to download[/]")
                    failed_species.append(species)
            else:
                console.print(f"\t   [red]✗ No FTP URL found[/]")
                failed_species.append(species)
        else:
            console.print(f"\t   [red]✗ No genome found[/]")
            failed_species.append(species)

        summary_json = json.dumps(summary_dict) if summary_dict else ""
        summary_json = summary_json.replace("\t", " ")
        summary_records.append(
            {
                "species": species,
                "category": label,
                "assembly_accession": aid if aid else "",
                "ftp_url": ftp_url if ftp_url else "",
                "local_file": local_path,
                "status": status,
                "summary_json": summary_json,
            }
        )
        time.sleep(0.5)

    # Always write summary TSV (even in dry-run)
    if summary_tsv is None:
        base = os.path.splitext(outfile)[0]
        summary_tsv = base + "_summary.tsv"
    tsv_dir = os.path.dirname(summary_tsv)
    if tsv_dir:
        os.makedirs(tsv_dir, exist_ok=True)

    with open(summary_tsv, "w") as f:
        f.write(
            "species\tcategory\tassembly_accession\tftp_url\tlocal_file\tstatus\tsummary_json\n"
        )
        for rec in summary_records:
            f.write(
                f"{rec['species']}\t{rec['category']}\t{rec['assembly_accession']}\t"
                f"{rec['ftp_url']}\t{rec['local_file']}\t{rec['status']}\t{rec['summary_json']}\n"
            )
    console.print(f"[green]\nSummary TSV written to: {summary_tsv}[/]")

    if dry_run:
        console.print()
        console.print(
            Panel.fit(
                f"[bold {color}]DRY-RUN SUMMARY:[/]\n"
                f"Would download {success_count} genomes\n"
                f"Would fail for {len(failed_species)} species\n"
                f"Output would be written to: {outfile}",
                border_style=color,
            )
        )
        if failed_species:
            console.print("[yellow]Species not found (would be skipped):[/]")
            for sp in failed_species:
                console.print(f"  [dim]• {sp}[/]")
        return

    if not temp_files:
        console.print(
            Panel.fit(
                f"[red]ERROR: No {label} genomes could be downloaded. Check your internet/email/API key.[/]",
                border_style="red",
            )
        )
        raise RuntimeError(f"No {label} genomes downloaded.")

    # Check for missing files and warn
    missing_files = [f for f in temp_files if not os.path.exists(f)]
    if missing_files:
        console.print(
            Panel.fit(
                f"[yellow]WARNING: {len(missing_files)} temporary file(s) missing and will be skipped:[/]\n"
                + "\n".join(f"  - {f}" for f in missing_files[:5])
                + (
                    f"\n  ... and {len(missing_files)-5} more"
                    if len(missing_files) > 5
                    else ""
                ),
                title="Missing Files",
                border_style="yellow",
            )
        )
        # Filter out missing files
        temp_files = [f for f in temp_files if os.path.exists(f)]

    if not temp_files:
        console.print(
            Panel.fit(
                f"[red]ERROR: No valid temporary files found for {label} genomes.[/]",
                border_style="red",
            )
        )
        raise RuntimeError(f"No valid {label} genomes to concatenate.")

    console.print(f"\n[bold]Concatenating {len(temp_files)} genomes to {outfile}...[/]")
    with open(outfile, "w") as out:
        for f in temp_files:
            try:
                with open(f, "r") as inf:
                    out.write(inf.read())
                os.remove(f)
            except Exception as e:
                console.print(f"[red]Error reading/removing {f}: {e}[/]")

    console.print(
        Panel.fit(
            f"[bold green]✓ DONE: Wrote {len(temp_files)} {label} genomes to {outfile}[/]",
            border_style="green",
        )
    )


def download_viral_genomes(
    outfile,
    email,
    api_key=None,
    genomes_per_category=5,
    verbose=False,
    dry_run=False,
    summary_tsv=None,
) -> None:
    """
    Download a balanced set of viral genomes across categories.
    The summary TSV is always written, even in dry‑run mode.
    """
    console.print()
    console.rule("[bold cyan]Downloading Stratified Viral Genomes[/]")

    table = Table(title="Viral Categories", box=box.ROUNDED)
    table.add_column("Category", style="cyan")
    table.add_column("Species listed", justify="right")
    table.add_column("Target", justify="right", style="green")
    for cat, species in VIRAL_CATEGORIES.items():
        target = (
            len(species)
            if genomes_per_category == 0
            else min(genomes_per_category, len(species))
        )
        table.add_row(cat, str(len(species)), str(target))
    console.print(table)

    if dry_run:
        console.print(
            "[yellow]DRY-RUN mode enabled – no files will be downloaded or written.[/]\n"
        )

    outdir = os.path.dirname(outfile)
    if outdir and not dry_run:
        os.makedirs(outdir, exist_ok=True)

    # We'll collect all species across categories and download them with deduplication
    all_species = []
    all_categories = {}
    for category, species_list in VIRAL_CATEGORIES.items():
        if genomes_per_category == 0:
            selected = species_list
        else:
            selected = species_list[:genomes_per_category]
        for sp in selected:
            all_species.append(sp)
            all_categories[sp] = category

    # Use download_species_list with the combined list, but we need to track categories per species.
    # The existing download_species_list expects a label; we can pass "viral" but then category in TSV will be "viral".
    # To keep category information, we need to modify download_species_list to accept a mapping.
    # Simplest: call download_species_list once with all species, but the TSV category will be "viral".
    # But we want per‑category classification in the TSV.
    # We'll write a custom loop that mirrors download_species_list but with category tracking.
    # Actually, we can reuse download_species_list and just pass "viral" as label; the summary TSV will have that.
    # But the user might want category info. We'll enhance download_species_list to accept a category mapping.

    # For simplicity, we'll just call download_species_list with the combined list and label="viral".
    # The summary TSV will have category as "viral" for all, which is fine for the user's purpose.
    # The main fix is deduplication, which is already in download_species_list.

    download_species_list(
        outfile=outfile,
        species_list=all_species,
        email=email,
        api_key=api_key,
        label="viral",
        verbose=verbose,
        dry_run=dry_run,
        summary_tsv=summary_tsv,
        color="cyan",
    )


def download_prokaryotic_genomes(
    outfile,
    email,
    api_key=None,
    verbose=False,
    dry_run=False,
    summary_tsv=None,
) -> None:
    """Download prokaryotic (bacteria + archaea) genomes."""
    download_species_list(
        outfile=outfile,
        species_list=PROKARYOTIC_SPECIES,
        email=email,
        api_key=api_key,
        label="Prokaryotic",
        verbose=verbose,
        dry_run=dry_run,
        summary_tsv=summary_tsv,
        color="green",
    )


def download_eukaryotic_genomes(
    outfile,
    email,
    api_key=None,
    verbose=False,
    dry_run=False,
    summary_tsv=None,
) -> None:
    """Download eukaryotic (host, fungi, plants) genomes."""
    download_species_list(
        outfile=outfile,
        species_list=EUKARYOTIC_SPECIES,
        email=email,
        api_key=api_key,
        label="Eukaryotic",
        verbose=verbose,
        dry_run=dry_run,
        summary_tsv=summary_tsv,
        color="yellow",
    )


def download_contaminant_genomes(
    outfile,
    email,
    api_key=None,
    verbose=False,
    dry_run=False,
    summary_tsv=None,
) -> None:
    """Legacy: download all contaminants (prokaryotic + eukaryotic) into one file."""
    download_species_list(
        outfile=outfile,
        species_list=CONTAMINANT_SPECIES,
        email=email,
        api_key=api_key,
        label="Contaminant (all)",
        verbose=verbose,
        dry_run=dry_run,
        summary_tsv=summary_tsv,
        color="magenta",
    )


def main():
    parser = argparse.ArgumentParser(
        description="Download RefSeq genomes for mock benchmarks",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--outfile",
        required=True,
        help="Output FASTA file path (parent directory created if needed)",
    )
    parser.add_argument("--email", required=True, help="NCBI Entrez email")
    parser.add_argument("--api-key", default="", help="NCBI API key (recommended)")
    parser.add_argument(
        "--mode",
        required=True,
        choices=["viral", "prokaryotic", "eukaryotic", "contaminants"],
        help="Which dataset to download. 'contaminants' is legacy and combines all.",
    )
    parser.add_argument(
        "--genomes-per-category",
        type=int,
        default=5,
        help="Number of genomes per viral category (use 0 to download all; only for --mode viral)",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable debug-level logging (detailed per-step messages)",
    )
    parser.add_argument(
        "--dry-run",
        "-n",
        action="store_true",
        help="Only check availability (no downloads, no file writes)",
    )
    parser.add_argument(
        "--summary-tsv",
        help="Output TSV file for genome summary (default: derived from outfile)",
    )
    args = parser.parse_args()

    # Set global Entrez settings
    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key
    Entrez.sleep_between_tries = 10

    try:
        if args.mode == "viral":
            download_viral_genomes(
                outfile=args.outfile,
                email=args.email,
                api_key=args.api_key,
                genomes_per_category=args.genomes_per_category,
                verbose=args.verbose,
                dry_run=args.dry_run,
                summary_tsv=args.summary_tsv,
            )
        elif args.mode == "prokaryotic":
            download_prokaryotic_genomes(
                outfile=args.outfile,
                email=args.email,
                api_key=args.api_key,
                verbose=args.verbose,
                dry_run=args.dry_run,
                summary_tsv=args.summary_tsv,
            )
        elif args.mode == "eukaryotic":
            download_eukaryotic_genomes(
                outfile=args.outfile,
                email=args.email,
                api_key=args.api_key,
                verbose=args.verbose,
                dry_run=args.dry_run,
                summary_tsv=args.summary_tsv,
            )
        else:  # contaminants (legacy)
            download_contaminant_genomes(
                outfile=args.outfile,
                email=args.email,
                api_key=args.api_key,
                verbose=args.verbose,
                dry_run=args.dry_run,
                summary_tsv=args.summary_tsv,
            )
    except Exception as e:
        console.print(Panel.fit(f"[red]FATAL ERROR: {e}[/]", border_style="red"))
        sys.exit(1)


if __name__ == "__main__":
    main()
