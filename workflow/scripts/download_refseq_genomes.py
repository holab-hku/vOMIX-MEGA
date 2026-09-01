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

# ----------------------------------------------------------------------
# Species lists
# ----------------------------------------------------------------------
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

CONTAMINANT_SPECIES = PROKARYOTIC_SPECIES + EUKARYOTIC_SPECIES


# ----------------------------------------------------------------------
# Entrez helpers
# ----------------------------------------------------------------------
def fetch_assembly_id(species_name, email, api_key=None, verbose=False):
    """Search RefSeq for the best complete genome assembly."""
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    queries = [
        f'"{species_name}"[Organism] AND "complete genome"[Assembly Level]',
        f'"{species_name}"[Organism] AND refseq[filter]',
    ]

    for attempt, query in enumerate(queries, 1):
        try:
            if verbose:
                console.log(f"[dim]Attempt #{attempt}: {query}[/]")
            handle = Entrez.esearch(db="assembly", term=query, retmax=1)
            record = Entrez.read(handle)
            handle.close()
            if record["IdList"]:
                if verbose:
                    console.log(f"[dim]   ✓ Found ID: {record['IdList'][0]}[/]")
                return record["IdList"][0]
            else:
                if verbose and attempt == 1:
                    console.log(
                        "[dim]   ✗ No complete genome, trying broader search[/]"
                    )
        except HTTPError as e:
            if e.code in (429, 500):
                time.sleep(attempt * 5)
                continue
            console.log(f"[red]Entrez error for {species_name}: {e}[/]")
            return None
        except Exception as e:
            console.log(f"[red]Error: {e}[/]")
            return None
        time.sleep(0.5)
    return None


def build_ftp_url(ftp_dir, assembly_accession, assembly_name):
    """Build the genomic.fna.gz URL from assembly metadata."""
    if not ftp_dir:
        return None
    https_dir = ftp_dir.replace("ftp://", "https://")
    if assembly_accession:
        # Build a clean name: remove spaces, underscores okay
        name_clean = assembly_name.replace(" ", "_")
        return f"{https_dir}/{assembly_accession}_{name_clean}_genomic.fna.gz"
    else:
        return f"{https_dir}/genomic.fna.gz"


def get_assembly_summary(assembly_id, email, api_key=None, verbose=False):
    """
    Retrieve assembly summary and FTP URL.
    Returns (ftp_url, summary_dict) or (None, None).
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    try:
        # Resolve accession to numeric ID if needed
        if assembly_id.startswith("GCF_") or assembly_id.startswith("GCA_"):
            handle = Entrez.esearch(db="assembly", term=assembly_id, retmax=1)
            record = Entrez.read(handle)
            handle.close()
            if not record["IdList"]:
                return None, None
            numeric_id = record["IdList"][0]
        else:
            numeric_id = assembly_id

        # Fetch summary
        handle = Entrez.esummary(db="assembly", id=numeric_id, report="full")
        summary = Entrez.read(handle)
        handle.close()
        if not summary:
            return None, None

        doc = summary["DocumentSummarySet"]["DocumentSummary"][0]
        ftp_dir = doc.get("FtpPath_RefSeq") or doc.get("FtpPath_GenBank")
        if ftp_dir:
            ftp_url = build_ftp_url(
                ftp_dir, doc.get("AssemblyAccession"), doc.get("AssemblyName", "")
            )
        else:
            # Fallback: try to construct URL from accession
            acc = doc.get("AssemblyAccession")
            name = doc.get("AssemblyName", "").replace(" ", "_")
            if acc:
                # Build base path (example: GCF_000865725.1 -> /GCF/000/865/725)
                parts = acc.split("_")
                if len(parts) >= 2:
                    number = parts[1]
                    groups = [number[i : i + 3] for i in range(0, len(number), 3)]
                    base = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{parts[0]}/{groups[0]}/{groups[1]}/{groups[2]}/{acc}_{name}"
                    ftp_url = f"{base}/{acc}_{name}_genomic.fna.gz"
                else:
                    ftp_url = None
            else:
                ftp_url = None

        if verbose and ftp_url:
            console.log(f"[dim][green]FTP URL: {ftp_url}[/]")
        return ftp_url, doc
    except Exception as e:
        if verbose:
            console.log(f"[red]Error fetching summary: {e}[/]")
        return None, None


# ----------------------------------------------------------------------
# Download helper
# ----------------------------------------------------------------------
def download_file(url, local_path, max_attempts=3, verbose=False):
    """Download a file with retries. Returns True on success."""
    for attempt in range(max_attempts):
        try:
            if verbose:
                console.log(f"[dim]Downloading: {url}[/]")
            urllib.request.urlretrieve(url, local_path)
            return True
        except HTTPError as e:
            if e.code == 404:
                if verbose:
                    console.log(f"[red]File not found: {url}[/]")
                return False
            if e.code in (429, 500):
                wait = (attempt + 1) * 5
                if verbose:
                    console.log(f"[yellow]Server error. Waiting {wait}s...[/]")
                time.sleep(wait)
                continue
            console.log(f"[red]HTTP error: {e}[/]")
            return False
        except Exception as e:
            console.log(f"[red]Download error: {e}[/]")
            return False
    return False


def decompress_gz(gz_path, fasta_path):
    """Decompress .gz to plain FASTA."""
    try:
        with gzip.open(gz_path, "rt") as gz_f, open(fasta_path, "w") as out_f:
            out_f.write(gz_f.read())
        os.remove(gz_path)
        return True
    except Exception as e:
        console.log(f"[red]Decompression error: {e}[/]")
        return False


def download_assembly_fasta(
    assembly_id, workdir, ftp_url, verbose=False, dry_run=False
):
    """
    Download and decompress a genome FASTA.
    Returns local FASTA path or None.
    """
    if dry_run:
        console.log(f"[dim][cyan]DRY-RUN: would download from {ftp_url}[/]")
        return "DRY_RUN_PLACEHOLDER"

    # Attempt primary URL, then fallback to genomic.fna.gz in same directory
    urls = [ftp_url]
    base_dir = ftp_url.rsplit("/", 1)[0]
    if not ftp_url.endswith("/genomic.fna.gz"):
        urls.append(f"{base_dir}/genomic.fna.gz")

    temp_gz = os.path.join(workdir, f".tmp_{assembly_id}.fna.gz")
    for url in urls:
        if download_file(url, temp_gz, verbose=verbose):
            temp_fasta = os.path.join(workdir, f".tmp_{assembly_id}.fna")
            if decompress_gz(temp_gz, temp_fasta):
                return temp_fasta
    return None


# ----------------------------------------------------------------------
# Species processing
# ----------------------------------------------------------------------
def process_species(
    species,
    email,
    api_key,
    outdir,
    verbose,
    dry_run,
    downloaded_assemblies,
    summary_records,
):
    """
    Process a single species: fetch ID, download if needed.
    Returns (success, temp_file_path) where temp_file_path is None on failure.
    """
    aid = fetch_assembly_id(species, email, api_key, verbose=verbose)
    if not aid:
        summary_records.append(
            {
                "species": species,
                "category": "",
                "assembly_accession": "",
                "ftp_url": "",
                "local_file": "N/A",
                "status": "failed",
                "summary_json": "",
            }
        )
        return False, None

    # Check if already downloaded
    if aid in downloaded_assemblies:
        fpath = downloaded_assemblies[aid]
        summary_records.append(
            {
                "species": species,
                "category": "",
                "assembly_accession": aid,
                "ftp_url": "",
                "local_file": fpath,
                "status": "success",
                "summary_json": "",
            }
        )
        return True, fpath

    ftp_url, summary_dict = get_assembly_summary(aid, email, api_key, verbose=verbose)
    if not ftp_url:
        summary_records.append(
            {
                "species": species,
                "category": "",
                "assembly_accession": aid,
                "ftp_url": "",
                "local_file": "N/A",
                "status": "failed",
                "summary_json": json.dumps(summary_dict) if summary_dict else "",
            }
        )
        return False, None

    fpath = download_assembly_fasta(
        aid, outdir, ftp_url, verbose=verbose, dry_run=dry_run
    )
    if fpath:
        downloaded_assemblies[aid] = fpath
        summary_records.append(
            {
                "species": species,
                "category": "",
                "assembly_accession": aid,
                "ftp_url": ftp_url,
                "local_file": fpath,
                "status": "success",
                "summary_json": json.dumps(summary_dict) if summary_dict else "",
            }
        )
        return True, fpath
    else:
        summary_records.append(
            {
                "species": species,
                "category": "",
                "assembly_accession": aid,
                "ftp_url": ftp_url,
                "local_file": "N/A",
                "status": "failed",
                "summary_json": json.dumps(summary_dict) if summary_dict else "",
            }
        )
        return False, None


# ----------------------------------------------------------------------
# Main download orchestrator
# ----------------------------------------------------------------------
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
    """Download genomes for a list of species."""
    console.print()
    console.rule(f"[bold {color}]Downloading {label} Genomes[/]")
    console.print(f"Target: {len(species_list)} species\n")

    if dry_run:
        console.print("[yellow]DRY-RUN mode – no files will be downloaded.[/]\n")

    outdir = os.path.dirname(outfile)
    if outdir and not dry_run:
        os.makedirs(outdir, exist_ok=True)

    temp_files = []
    success_count = 0
    downloaded_assemblies = {}
    summary_records = []
    failed_species = []

    for species in species_list:
        console.print(f"\n[bold]→ {species}[/]")
        console.log("[dim]Fetching...[/]")
        success, fpath = process_species(
            species,
            email,
            api_key,
            outdir or ".",
            verbose,
            dry_run,
            downloaded_assemblies,
            summary_records,
        )
        if success and fpath:
            # Avoid duplicates
            if fpath not in temp_files:
                temp_files.append(fpath)
            success_count += 1
            if not dry_run:
                console.print(f"\t   [green]✓ Downloaded ({success_count})[/]")
            else:
                console.print(f"\t   [green]✓ Would download ({success_count})[/]")
        else:
            failed_species.append(species)
            console.print(f"\t   [red]✗ Failed[/]")
        time.sleep(0.5)

    # Write summary TSV
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
        console.print(
            Panel.fit(
                f"[bold {color}]DRY-RUN: Would download {success_count} genomes, {len(failed_species)} failures[/]",
                border_style=color,
            )
        )
        return

    # Check for missing files
    valid_files = [f for f in temp_files if os.path.exists(f)]
    if len(valid_files) < len(temp_files):
        missing = set(temp_files) - set(valid_files)
        console.print(
            f"[yellow]WARNING: {len(missing)} temporary file(s) missing and skipped.[/]"
        )

    if not valid_files:
        console.print(
            Panel.fit(
                f"[red]ERROR: No {label} genomes downloaded.[/]", border_style="red"
            )
        )
        raise RuntimeError(f"No {label} genomes downloaded.")

    # Concatenate
    console.print(
        f"\n[bold]Concatenating {len(valid_files)} genomes to {outfile}...[/]"
    )
    with open(outfile, "w") as out:
        for f in valid_files:
            try:
                with open(f, "r") as inf:
                    out.write(inf.read())
                os.remove(f)
            except Exception as e:
                console.print(f"[red]Error reading/removing {f}: {e}[/]")

    console.print(
        Panel.fit(
            f"[bold green]✓ DONE: Wrote {len(valid_files)} {label} genomes to {outfile}[/]",
            border_style="green",
        )
    )


# ----------------------------------------------------------------------
# Wrapper functions for each mode
# ----------------------------------------------------------------------
def download_viral_genomes(
    outfile,
    email,
    api_key=None,
    genomes_per_category=5,
    verbose=False,
    dry_run=False,
    summary_tsv=None,
):
    """Download stratified viral genomes."""
    console.print()
    console.rule("[bold cyan]Downloading Stratified Viral Genomes[/]")

    table = Table(title="Viral Categories", box=box.ROUNDED)
    table.add_column("Category", style="cyan")
    table.add_column("Species listed", justify="right")
    table.add_column("Target", justify="right", style="green")
    all_species = []
    for cat, species in VIRAL_CATEGORIES.items():
        target = (
            len(species)
            if genomes_per_category == 0
            else min(genomes_per_category, len(species))
        )
        selected = species[:target]
        all_species.extend(selected)
        table.add_row(cat, str(len(species)), str(target))
    console.print(table)

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
    outfile, email, api_key=None, verbose=False, dry_run=False, summary_tsv=None
):
    download_species_list(
        outfile,
        PROKARYOTIC_SPECIES,
        email,
        api_key,
        "Prokaryotic",
        verbose,
        dry_run,
        summary_tsv,
        "green",
    )


def download_eukaryotic_genomes(
    outfile, email, api_key=None, verbose=False, dry_run=False, summary_tsv=None
):
    download_species_list(
        outfile,
        EUKARYOTIC_SPECIES,
        email,
        api_key,
        "Eukaryotic",
        verbose,
        dry_run,
        summary_tsv,
        "yellow",
    )


def download_contaminant_genomes(
    outfile, email, api_key=None, verbose=False, dry_run=False, summary_tsv=None
):
    download_species_list(
        outfile,
        CONTAMINANT_SPECIES,
        email,
        api_key,
        "Contaminants (all)",
        verbose,
        dry_run,
        summary_tsv,
        "magenta",
    )


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Download RefSeq genomes for mock benchmarks"
    )
    parser.add_argument("--outfile", required=True, help="Output FASTA file path")
    parser.add_argument("--email", required=True, help="NCBI Entrez email")
    parser.add_argument("--api-key", default="", help="NCBI API key (recommended)")
    parser.add_argument(
        "--mode",
        required=True,
        choices=["viral", "prokaryotic", "eukaryotic", "contaminants"],
        help="Dataset to download",
    )
    parser.add_argument(
        "--genomes-per-category",
        type=int,
        default=5,
        help="Number per viral category (0=all; only for --mode viral)",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable debug logging"
    )
    parser.add_argument("-n", "--dry-run", action="store_true", help="Dry run only")
    parser.add_argument("--summary-tsv", help="Output summary TSV path")
    args = parser.parse_args()

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
                args.outfile,
                args.email,
                args.api_key,
                args.verbose,
                args.dry_run,
                args.summary_tsv,
            )
        elif args.mode == "eukaryotic":
            download_eukaryotic_genomes(
                args.outfile,
                args.email,
                args.api_key,
                args.verbose,
                args.dry_run,
                args.summary_tsv,
            )
        elif args.mode == "contaminants":
            download_contaminant_genomes(
                args.outfile,
                args.email,
                args.api_key,
                args.verbose,
                args.dry_run,
                args.summary_tsv,
            )
    except Exception as e:
        console.print(Panel.fit(f"[red]FATAL ERROR: {e}[/]", border_style="red"))
        sys.exit(1)


if __name__ == "__main__":
    main()
