#!/usr/bin/env python3
"""
generate_mock.py - Generate a mock FASTA dataset for benchmarking.

Inputs:
  --viral-seq          FASTA with viral genomes (required)
  --prokaryotic-seq    FASTA with prokaryotic genomes (optional; synthetic fallback)
  --eukaryotic-seq     FASTA with eukaryotic genomes (optional; synthetic fallback)

Fractions (must sum to 1.0):
  --virus-frac         fraction of viral contigs (default 0.5)
  --prokaryote-frac    fraction of prokaryotic contigs (default 0.3)
  --eukaryote-frac     fraction of eukaryotic contigs (default 0.2)

Outputs:
  - FASTA file (<name>.fna)
  - Ground truth TSV (<name>.ground_truth.tsv)

Labeling:
  - Viral contigs → source genome accession
  - Strain-mode contigs → parent accession (all strains share label)
  - Prokaryotic/Eukaryotic contigs → source genome accession (or self if synthetic)
  - Duplicates → same label as original

All proportions, mutation rates, duplication, and overlap are adjustable.
"""

import os
import sys
import math
import random
import argparse
from typing import List, Tuple, Optional
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Rich imports
try:
    from rich.console import Console
    from rich.panel import Panel
    from rich.progress import Progress, SpinnerColumn, TextColumn
    from rich import box
except ImportError:
    print("Rich library not installed. Install with: pip install rich", file=sys.stderr)
    sys.exit(1)

console = Console()


def read_fasta_seq(seq_file: str) -> List[SeqRecord]:
    """Return a list of SeqRecord objects from a FASTA file."""
    if not os.path.exists(seq_file) or os.path.getsize(seq_file) == 0:
        return []
    return list(SeqIO.parse(seq_file, "fasta"))


def fragment_sequence(
    seq_record: SeqRecord,
    min_len: int = 500,
    max_len: int = 50000,
    num_fragments: int = 1,
    lognormal_mu: float = 8.5,
    lognormal_sigma: float = 1.2,
) -> List[SeqRecord]:
    """
    Generate num_fragments contigs from a single sequence.
    Fragment lengths are drawn from a truncated lognormal distribution
    (parameters mu, sigma) bounded by [min_len, max_len].
    Overlap is not explicitly enforced but can occur naturally.
    """
    seq_len = len(seq_record.seq)
    if seq_len < min_len:
        return [seq_record]

    effective_max = min(max_len, seq_len)
    if effective_max < min_len:
        effective_max = seq_len

    fragments = []

    for _ in range(num_fragments):
        length = None
        for _ in range(100):  # safety limit
            z = random.gauss(0, 1)
            candidate = int(round(math.exp(lognormal_mu + lognormal_sigma * z)))
            if min_len <= candidate <= effective_max:
                length = candidate
                break
        if length is None:
            length = random.randint(min_len, effective_max)

        max_start = max(0, seq_len - length)
        start = random.randint(0, max_start)
        end = start + length
        fragments.append(seq_record[start:end])

    return fragments if fragments else [seq_record]


def generate_strains(
    species_seq: SeqRecord,
    source_id: str,
    num_strains: int,
    mut_rate_min: float = 0.001,
    mut_rate_max: float = 0.05,
) -> List[Tuple[SeqRecord, str]]:
    """
    Generate mutated strains from a single species sequence.

    Args:
        species_seq (SeqRecord): The parent genome sequence.
        source_id (str): Cluster label for all strains.
        num_strains (int): Number of strains to generate.
        mut_rate_min (float): Minimum mutation rate.
        mut_rate_max (float): Maximum mutation rate.

    Returns:
        List[Tuple[SeqRecord, str]]: List of (strain_record, cluster_label) tuples.
    """
    seq_str = str(species_seq.seq)
    strains = []

    for strain_idx in range(num_strains):
        mut_rate = random.uniform(mut_rate_min, mut_rate_max)
        mutated = list(seq_str)
        for i in range(len(mutated)):
            if random.random() < mut_rate:
                mutated[i] = random.choice("ACGT")
        new_seq = Seq("".join(mutated))
        new_id = f"{source_id}_strain_{strain_idx+1}"
        new_record = SeqRecord(
            new_seq,
            id=new_id,
            description=f"ANI={1-mut_rate:.3f} source={source_id}",
        )
        strains.append((new_record, source_id))

    return strains


def add_contaminant_contigs(
    records: List[SeqRecord],
    ground_truth: List[Tuple[str, str]],
    source_records: List[SeqRecord],
    needed: int,
    label_prefix: str,
    use_synthetic: bool = False,
    min_len: int = 500,
    max_len: int = 50000,
) -> Tuple[List[SeqRecord], List[Tuple[str, str]]]:
    """
    Add contaminant (prokaryotic/eukaryotic) contigs to the dataset.

    If source_records is non‑empty and use_synthetic=False, fragment real genomes.
    Otherwise, generate synthetic random sequences.

    Args:
        records: Current list of SeqRecord objects.
        ground_truth: Current ground truth list.
        source_records: List of SeqRecord objects to fragment.
        needed: Number of contigs to add.
        label_prefix: Prefix for synthetic IDs (e.g., "bacteria", "eukaryote").
        use_synthetic: If True, generate synthetic sequences regardless of source_records.
        min_len, max_len: Fragment length bounds.

    Returns:
        Updated records and ground truth.
    """
    if needed <= 0:
        return records, ground_truth

    console.print(f"[cyan]Adding {needed} {label_prefix} sequences...[/]")

    if not use_synthetic and source_records:
        # Use real genomes
        while len(records) < needed:
            rec = random.choice(source_records)
            num_frags = random.randint(1, 3)
            frags = fragment_sequence(rec, min_len=min_len, max_len=max_len, num_fragments=num_frags)
            for frag in frags:
                new_id = f"{rec.id}_{label_prefix}_{len(records)+1}"
                frag.id = new_id
                frag.description = f"source={rec.id}"
                records.append(frag)
                ground_truth.append((new_id, rec.id))
                if len(records) >= needed:
                    break
            if len(records) >= needed:
                break
    else:
        # Synthetic sequences
        for i in range(needed):
            # Choose sub‑category (e.g., for prokaryotic, maybe different GC content? but simple)
            length = random.randint(min_len, max_len)
            seq = ''.join(random.choices("ACGT", k=length))
            new_id = f"{label_prefix}_synthetic_{i+1}"
            rec = SeqRecord(Seq(seq), id=new_id, description=f"length={length}")
            records.append(rec)
            ground_truth.append((new_id, new_id))  # self as cluster

    return records, ground_truth


def duplicate_records(
    records: List[SeqRecord],
    ground_truth: List[Tuple[str, str]],
    duplication_factor: int,
) -> Tuple[List[SeqRecord], List[Tuple[str, str]]]:
    """
    Duplicate each record `duplication_factor` times.

    Args:
        records (List[SeqRecord]): Original records.
        ground_truth (List[Tuple[str, str]]): Original ground truth (seq_id, cluster).
        duplication_factor (int): Number of copies.

    Returns:
        Tuple[List[SeqRecord], List[Tuple[str, str]]]: Duplicated records and ground truth.
    """
    if duplication_factor <= 1:
        return records, ground_truth

    console.print(f"[yellow]Duplicating fragments {duplication_factor} times...[/]")
    duplicated_records = []
    duplicated_gt = []

    for rec, (seq_id, cluster) in zip(records, ground_truth):
        for copy in range(duplication_factor):
            new_rec = rec[:]
            new_rec.id = f"{seq_id}_dup{copy+1}"
            new_rec.description = rec.description + f" copy={copy+1}"
            duplicated_records.append(new_rec)
            duplicated_gt.append((new_rec.id, cluster))

    return duplicated_records, duplicated_gt


def generate_mock_dataset(
    name: str,
    total_seqs: int,
    outdir: str,
    tmpdir: str,
    viral_seq: str,
    prokaryotic_seq: str,
    eukaryotic_seq: str,
    virus_frac: float = 0.5,
    prokaryote_frac: float = 0.3,
    eukaryote_frac: float = 0.2,
    strain_mode: bool = False,
    mut_rate_min: float = 0.001,
    mut_rate_max: float = 0.05,
    duplication_factor: int = 1,
    overlap_min: int = 0,
    overlap_max: int = 0,
    num_species: int = 10,
    seed: int = 42,
    ground_truth_file: Optional[str] = None,
    no_ground_truth: bool = False,
    force: bool = False,
    verbose: bool = False,
    dry_run: bool = False,
) -> None:
    """
    Main random contig generation function.
    """
    # Validate fractions
    total_frac = virus_frac + prokaryote_frac + eukaryote_frac
    if not abs(total_frac - 1.0) < 1e-9:
        raise ValueError(
            f"The fractions must sum to 1.0. Got virus={virus_frac}, "
            f"prokaryote={prokaryote_frac}, eukaryote={eukaryote_frac} (sum={total_frac:.10f})."
        )

    random.seed(seed)
    outfile = os.path.join(outdir, f"{name}.fna")
    gt_file = ground_truth_file or outfile.replace(".fna", ".ground_truth.tsv")

    if not force and os.path.exists(outfile) and os.path.getsize(outfile) > 1000:
        console.print(f"[yellow][SKIP] {name} already exists. Use --force to regenerate.[/]")
        if not no_ground_truth and not os.path.exists(gt_file):
            console.print("[red][WARN] Ground truth missing; please regenerate with --force.[/]")
        return

    if dry_run:
        console.print("[cyan]DRY-RUN: Would generate dataset with the following parameters:[/]")
        console.print(f"  Name: {name}, total sequences: {total_seqs}")
        console.print(f"  Virus fraction: {virus_frac}, Prokaryote: {prokaryote_frac}, Eukaryote: {eukaryote_frac}")
        console.print(f"  Strain mode: {strain_mode}, duplications: {duplication_factor}")
        console.print(f"  Output: {outfile}, ground truth: {gt_file}")
        return

    console.print(f"[green][GEN] Generating {name} with {total_seqs} sequences...[/]")

    # Load sequences
    viral_records = read_fasta_seq(viral_seq)
    if not viral_records:
        raise RuntimeError("Viral sequences file is empty or missing. Please run download_refseq_viral rule first.")

    prokaryotic_records = read_fasta_seq(prokaryotic_seq) if prokaryotic_seq else []
    eukaryotic_records = read_fasta_seq(eukaryotic_seq) if eukaryotic_seq else []

    # Compute target counts
    n_virus = int(total_seqs * virus_frac)
    n_prok = int(total_seqs * prokaryote_frac)
    n_euk = total_seqs - n_virus - n_prok  # remainder

    records = []
    ground_truth = []

    # ----- Viral sequences -----
    if strain_mode:
        if len(viral_records) < num_species:
            raise ValueError(
                f"Not enough viral genomes for strain mode (need {num_species}, have {len(viral_records)})."
            )
        selected = random.sample(viral_records, num_species)
        per_species = n_virus // num_species
        remainder = n_virus % num_species

        console.print(f"[cyan]Strain mode: {num_species} species, {per_species} strains each (+ {remainder} extra)[/]")

        for species_idx, species_seq in enumerate(selected):
            source_id = species_seq.id
            strains_needed = per_species + (1 if species_idx < remainder else 0)
            if strains_needed <= 0:
                continue

            # Generate full-length strains
            strain_tuples = generate_strains(
                species_seq,
                source_id,
                strains_needed,
                mut_rate_min,
                mut_rate_max,
            )

            # For each strain, fragment into 1-5 contigs
            for strain_rec, cluster_label in strain_tuples:
                num_frags = random.randint(1, 5)
                frags = fragment_sequence(strain_rec, num_fragments=num_frags)
                for frag in frags:
                    new_id = f"{strain_rec.id}_frag_{len(records)+1}"
                    frag.id = new_id
                    frag.description = (
                        f"strain={strain_rec.id} source={source_id} "
                        f"ANI={strain_rec.description.split('ANI=')[-1].split()[0]}"
                    )
                    records.append(frag)
                    ground_truth.append((new_id, cluster_label))
                    if len(records) >= n_virus:
                        break
                if len(records) >= n_virus:
                    break
            if len(records) >= n_virus):
                break

        # Fill remaining viral count with non‑strain fragments if needed
        while len(records) < n_virus:
            rec = random.choice(viral_records)
            frags = fragment_sequence(rec, num_fragments=1)
            for frag in frags:
                new_id = f"{rec.id}_frag_{len(records)+1}"
                frag.id = new_id
                frag.description = f"source={rec.id}"
                records.append(frag)
                ground_truth.append((new_id, rec.id))
                if len(records) >= n_virus:
                    break
    else:
        # Normal mode: fragment viral genomes
        console.print("[cyan]Normal mode: fragmenting viral genomes into contigs[/]")
        while len(records) < n_virus:
            rec = random.choice(viral_records)
            num_frags = random.randint(1, 5)
            frags = fragment_sequence(rec, num_fragments=num_frags)
            for frag in frags:
                new_id = f"{rec.id}_frag_{len(records)+1}"
                frag.id = new_id
                frag.description = f"source={rec.id}"
                records.append(frag)
                ground_truth.append((new_id, rec.id))
                if len(records) >= n_virus:
                    break
            if len(records) >= n_virus:
                break

    # Trim to exactly n_virus (if overshot)
    if len(records) > n_virus:
        records = records[:n_virus]
        ground_truth = ground_truth[:n_virus]

    # ----- Prokaryotic sequences -----
    # Use real records if available, otherwise synthetic
    use_synthetic_prok = not bool(prokaryotic_records)
    records, ground_truth = add_contaminant_contigs(
        records, ground_truth,
        prokaryotic_records,
        n_prok,
        label_prefix="prokaryote",
        use_synthetic=use_synthetic_prok,
    )

    # ----- Eukaryotic sequences -----
    use_synthetic_euk = not bool(eukaryotic_records)
    records, ground_truth = add_contaminant_contigs(
        records, ground_truth,
        eukaryotic_records,
        n_euk,
        label_prefix="eukaryote",
        use_synthetic=use_synthetic_euk,
    )

    # ----- Duplication -----
    records, ground_truth = duplicate_records(records, ground_truth, duplication_factor)

    # Shuffle records (but keep ground truth in sync)
    combined = list(zip(records, ground_truth))
    random.shuffle(combined)
    records, ground_truth = zip(*combined) if combined else ([], [])

    # Write FASTA
    os.makedirs(outdir, exist_ok=True)
    SeqIO.write(records, outfile, "fasta")
    console.print(f"[green][DONE] Wrote {len(records)} sequences to {outfile}[/]")

    # Write ground truth TSV (unless disabled)
    if not no_ground_truth:
        with open(gt_file, "w") as f:
            f.write("sequence_id\ttrue_cluster\n")
            for seq_id, cluster in ground_truth:
                f.write(f"{seq_id}\t{cluster}\n")
        console.print(f"[green][DONE] Wrote ground truth to {gt_file}[/]")


def main():
    parser = argparse.ArgumentParser(
        description="Generate a mock contig FASTA datasets with ground truth inputs (with viral strain generation functionality).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--name", required=True, help="Dataset name (e.g., Mock-10K)")
    parser.add_argument("--num-sequences", type=int, required=True, help="Total number of sequences")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--tmpdir", required=True, help="Temporary directory (unused but kept for consistency)")
    parser.add_argument("--viral-seq", required=True, help="FASTA file with viral genomes")
    parser.add_argument("--prokaryotic-seq", default="", help="FASTA file with prokaryotic genomes")
    parser.add_argument("--eukaryotic-seq", default="", help="FASTA file with eukaryotic genomes")
    parser.add_argument("--virus-frac", type=float, default=0.5, help="Fraction of viral sequences")
    parser.add_argument("--prokaryote-frac", type=float, default=0.3, help="Fraction of prokaryotic sequences")
    parser.add_argument("--eukaryote-frac", type=float, default=0.2, help="Fraction of eukaryotic sequences")
    parser.add_argument("--strain-mode", type=int, default=0, help="Enable strain variation (1) or not (0)")
    parser.add_argument("--mut-rate-min", type=float, default=0.001, help="Minimum mutation rate (0.001 = 99.9% ANI)")
    parser.add_argument("--mut-rate-max", type=float, default=0.05, help="Maximum mutation rate (0.05 = 95% ANI)")
    parser.add_argument("--duplication-factor", type=int, default=1, help="Number of copies of each fragment")
    parser.add_argument("--overlap-min", type=int, default=0, help="Minimum overlap between fragments (bp)")
    parser.add_argument("--overlap-max", type=int, default=0, help="Maximum overlap between fragments (bp)")
    parser.add_argument("--num-species", type=int, default=10, help="Number of species for strain mode")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--ground-truth", help="Output path for ground truth TSV (default: auto)")
    parser.add_argument("--no-ground-truth", action="store_true", help="Do not write ground truth file")
    parser.add_argument("--force", action="store_true", help="Regenerate even if output exists")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
    parser.add_argument("--dry-run", "-n", action="store_true", help="Preview actions without writing files")
    args = parser.parse_args()

    # Validate fractions
    total_frac = args.virus_frac + args.prokaryote_frac + args.eukaryote_frac
    if abs(total_frac - 1.0) > 1e-6:
        raise ValueError(f"Fractions must sum to 1.0 (got {total_frac})")

    # Validate input files (if not dry‑run)
    if not args.dry_run:
        if not os.path.exists(args.viral_seq):
            raise FileNotFoundError(f"Viral sequence file not found: {args.viral_seq}")
        # Prokaryotic and eukaryotic files are optional; if provided, check existence
        if args.prokaryotic_seq and not os.path.exists(args.prokaryotic_seq):
            raise FileNotFoundError(f"Prokaryotic sequence file not found: {args.prokaryotic_seq}")
        if args.eukaryotic_seq and not os.path.exists(args.eukaryotic_seq):
            raise FileNotFoundError(f"Eukaryotic sequence file not found: {args.eukaryotic_seq}")

    generate_mock_dataset(
        name=args.name,
        total_seqs=args.num_sequences,
        outdir=args.outdir,
        tmpdir=args.tmpdir,
        viral_seq=args.viral_seq,
        prokaryotic_seq=args.prokaryotic_seq,
        eukaryotic_seq=args.eukaryotic_seq,
        virus_frac=args.virus_frac,
        prokaryote_frac=args.prokaryote_frac,
        eukaryote_frac=args.eukaryote_frac,
        strain_mode=bool(args.strain_mode),
        mut_rate_min=args.mut_rate_min,
        mut_rate_max=args.mut_rate_max,
        duplication_factor=args.duplication_factor,
        overlap_min=args.overlap_min,
        overlap_max=args.overlap_max,
        num_species=args.num_species,
        seed=args.seed,
        ground_truth_file=args.ground_truth,
        no_ground_truth=args.no_ground_truth,
        force=args.force,
        verbose=args.verbose,
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    main()