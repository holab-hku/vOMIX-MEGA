#!/usr/bin/env python3
"""
download_refseq_genomes.py - Download stratified RefSeq genomes for mock benchmarks.

Modes:
  --mode viral         : Downloads genomes from dsDNA, ssDNA, dsRNA, ssRNA.
  --mode contaminants  : Downloads a curated list of bacterial, archaeal, and eukaryotic genomes.

Output:
  Writes a single FASTA file to the path specified by --outfile.
  The parent directory is created if it does not exist.

Options:
  --genomes-per-category 0 : downloads all species listed for that category.
  --verbose / -v          : enables debug-level logging.

Requires:
  - Biopython (Entrez)
  - requests
  - rich (for pretty logging)
"""

import os
import sys
import time
import gzip
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

# ============================================================
# 1. CURATED VIRAL SPECIES LISTS (unique, representative)
# ============================================================

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
        "Influenza A virus",
        "Human immunodeficiency virus 1",
        "Severe acute respiratory syndrome coronavirus 2",
        "Poliovirus",
        "West Nile virus",
        "Dengue virus type 2",
        "Zika virus",
        "Rabies virus",
        "Ebola virus",
        "Measles virus",
    ],
}

# Common metagenomic contaminants (bacteria, archaea, eukaryotes)
CONTAMINANT_SPECIES = [
    # Bacteria (gram-negative)
    "Escherichia coli",
    "Pseudomonas aeruginosa",
    "Acinetobacter baumannii",
    "Salmonella enterica",
    "Bacteroides fragilis",
    # Bacteria (gram-positive)
    "Staphylococcus aureus",
    "Bacillus subtilis",
    "Lactobacillus acidophilus",
    "Streptococcus pneumoniae",
    "Mycobacterium tuberculosis",
    # Archaea
    "Methanobrevibacter smithii",
    "Halobacterium salinarum",
    # Eukaryotes (host & fungi)
    "Homo sapiens",
    "Saccharomyces cerevisiae",
    "Candida albicans",
    "Aspergillus fumigatus",
    "Arabidopsis thaliana",
]

# ============================================================
# 2. NCBI Entrez Helper Functions
# ============================================================


def fetch_assembly_id(species_name, email, api_key=None, max_attempts=3, verbose=False):
    """
    Search RefSeq for the best complete genome assembly of a given species.
    Returns assembly ID (string) or None if not found.
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    query = f'"{species_name}"[Organism] AND "complete genome"[Assembly Level] AND refseq[filter]'

    for attempt in range(max_attempts):
        try:
            if verbose:
                console.log(f"[dim]Querying Entrez for: {species_name}[/]")
            handle = Entrez.esearch(db="assembly", term=query, retmax=1)
            record = Entrez.read(handle)
            handle.close()
            if record["IdList"]:
                return record["IdList"][0]
            else:
                # Retry without complete genome filter
                if attempt == 0:
                    query = f'"{species_name}"[Organism] AND refseq[filter]'
                    continue
                return None
        except HTTPError as e:
            if e.code in (429, 500):
                wait = (attempt + 1) * 5
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
    return None


def download_assembly_fasta(
    assembly_id, workdir, email, api_key=None, max_attempts=3, verbose=False
):
    """
    Download the genomic FASTA (.fna.gz) for a given assembly ID.
    Writes a temporary file in workdir with a .tmp_ prefix.
    Returns path to that temporary file or None on failure.
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    for attempt in range(max_attempts):
        try:
            handle = Entrez.esummary(db="assembly", id=assembly_id)
            summary = Entrez.read(handle)
            handle.close()
            if not summary:
                if verbose:
                    console.log(f"[red]No summary for {assembly_id}[/]")
                return None

            ftp_path = summary[0].get("FtpPath_RefSeq", "")
            if not ftp_path:
                ftp_path = summary[0].get("FtpPath_GenBank", "")
            if not ftp_path:
                if verbose:
                    console.log(f"[red]No FTP path for {assembly_id}[/]")
                return None

            base = ftp_path.split("/")[-1]
            fna_url = f"{ftp_path}/{base}_genomic.fna.gz"
            temp_gz = os.path.join(workdir, f".tmp_{assembly_id}.fna.gz")
            if verbose:
                console.log(f"[dim]Downloading: {fna_url}[/]")
            urllib.request.urlretrieve(fna_url, temp_gz)

            # Decompress to plain FASTA
            temp_fasta = os.path.join(workdir, f".tmp_{assembly_id}.fna")
            with gzip.open(temp_gz, "rt") as gz_f, open(temp_fasta, "w") as out_f:
                out_f.write(gz_f.read())
            os.remove(temp_gz)
            return temp_fasta

        except HTTPError as e:
            if e.code == 404:
                if verbose:
                    console.log(f"[red]File not found for {assembly_id}[/]")
                return None
            if e.code in (429, 500):
                wait = (attempt + 1) * 5
                if verbose:
                    console.log(f"[yellow]Server error. Waiting {wait}s...[/]")
                time.sleep(wait)
                continue
            console.log(f"[red]HTTP error for {assembly_id}: {e}[/]")
            return None
        except Exception as e:
            console.log(f"[red]Error downloading {assembly_id}: {e}[/]")
            return None
    return None


# ============================================================
# 3. Main download functions
# ============================================================


def download_viral_genomes(
    outfile, email, api_key=None, genomes_per_category=5, verbose=False
):
    """
    Download a balanced set of viral genomes across categories.
    If genomes_per_category == 0, download all species in each category.
    Writes final file to outfile, creating parent directory if needed.
    """
    console.print()
    console.rule("[bold cyan]Downloading Stratified Viral Genomes[/]")

    # Show summary table
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

    # Ensure output directory exists
    outdir = os.path.dirname(outfile)
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    temp_files = []
    success_count = 0

    for category, species_list in VIRAL_CATEGORIES.items():
        if genomes_per_category == 0:
            selected = species_list
        else:
            selected = species_list[:genomes_per_category]

        console.print(f"\n[bold]{category}[/] – targeting {len(selected)} genomes")

        for species in selected:
            console.print(f"  [dim]→ {species}[/]")
            aid = fetch_assembly_id(species, email, api_key, verbose=verbose)
            if aid:
                fname = download_assembly_fasta(
                    aid, outdir, email, api_key, verbose=verbose
                )
                if fname:
                    temp_files.append(fname)
                    success_count += 1
                    console.print(f"    [green]✓ Downloaded ({success_count})[/]")
                else:
                    console.print(f"    [red]✗ Failed to download[/]")
            else:
                console.print(f"    [red]✗ No complete assembly found[/]")
            time.sleep(0.5)  # NCBI rate limit

    if not temp_files:
        console.print(
            Panel.fit(
                "[red]ERROR: No viral genomes could be downloaded. Check your internet/email/API key.[/]",
                border_style="red",
            )
        )
        raise RuntimeError("No viral genomes downloaded.")

    # Concatenate to outfile
    console.print(f"\n[bold]Concatenating {len(temp_files)} genomes to {outfile}...[/]")
    with open(outfile, "w") as out:
        for f in temp_files:
            with open(f, "r") as inf:
                out.write(inf.read())
            os.remove(f)  # clean up temporary file

    console.print(
        Panel.fit(
            f"[bold green]✓ DONE: Wrote {success_count} viral genomes to {outfile}[/]",
            border_style="green",
        )
    )


def download_contaminant_genomes(outfile, email, api_key=None, verbose=False):
    """
    Download the curated contaminant genomes.
    Writes final file to outfile, creating parent directory if needed.
    """
    console.print()
    console.rule("[bold magenta]Downloading Contaminant Genomes[/]")
    console.print(f"Target: {len(CONTAMINANT_SPECIES)} species\n")

    # Ensure output directory exists
    outdir = os.path.dirname(outfile)
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    temp_files = []
    success_count = 0

    for species in CONTAMINANT_SPECIES:
        console.print(f"  [dim]→ {species}[/]")
        aid = fetch_assembly_id(species, email, api_key, verbose=verbose)
        if aid:
            fname = download_assembly_fasta(
                aid, outdir, email, api_key, verbose=verbose
            )
            if fname:
                temp_files.append(fname)
                success_count += 1
                console.print(f"    [green]✓ Downloaded ({success_count})[/]")
            else:
                console.print(f"    [red]✗ Failed to download[/]")
        else:
            console.print(f"    [red]✗ No complete assembly found[/]")
        time.sleep(0.5)

    if not temp_files:
        console.print(
            Panel.fit(
                "[red]ERROR: No contaminant genomes could be downloaded.[/]",
                border_style="red",
            )
        )
        raise RuntimeError("No contaminant genomes downloaded.")

    console.print(f"\n[bold]Concatenating {len(temp_files)} genomes to {outfile}...[/]")
    with open(outfile, "w") as out:
        for f in temp_files:
            with open(f, "r") as inf:
                out.write(inf.read())
            os.remove(f)  # clean up temporary file

    console.print(
        Panel.fit(
            f"[bold green]✓ DONE: Wrote {success_count} contaminant genomes to {outfile}[/]",
            border_style="green",
        )
    )


# ============================================================
# 4. Command-line interface
# ============================================================


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
        choices=["viral", "contaminants"],
        help="Which dataset to download",
    )
    parser.add_argument(
        "--genomes-per-category",
        type=int,
        default=5,
        help="Number of genomes per viral category (use 0 to download all)",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable debug-level logging (detailed per-step messages)",
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
            )
        else:  # contaminants
            download_contaminant_genomes(
                outfile=args.outfile,
                email=args.email,
                api_key=args.api_key,
                verbose=args.verbose,
            )
    except Exception as e:
        console.print(Panel.fit(f"[red]FATAL ERROR: {e}[/]", border_style="red"))
        sys.exit(1)


if __name__ == "__main__":
    main()
