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

Fragment length distribution (lognormal):
  --lognormal-mu       mean of log-length (default 8.5, median ~5 kb)
  --lognormal-sigma    spread of log-length (default 1.2)

Outputs:
  - FASTA file (<name>.fna)
  - Ground truth TSV (<name>.ground_truth.tsv) with extensive metadata

All proportions, mutation rates, and duplications are adjustable.
"""

import os
import sys
import math
import random
import argparse
from typing import List, Tuple, Optional, Any
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Rich imports
try:
    from rich.console import Console
    from rich.panel import Panel
    from rich import box
except ImportError:
    print("Rich library not installed. Install with: pip install rich", file=sys.stderr)
    sys.exit(1)

console = Console()


# ----------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------
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
    bounded by [min_len, max_len].
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
) -> List[Tuple[SeqRecord, str, float]]:
    """
    Generate mutated strains from a single species sequence.

    Returns:
        List of (strain_record, source_id, mutation_rate)
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
        strains.append((new_record, source_id, mut_rate))

    return strains


def add_contaminant_contigs(
    records: List[SeqRecord],
    ground_truth: List[Tuple],
    source_records: List[SeqRecord],
    needed: int,
    label_prefix: str,
    category: str,
    use_synthetic: bool = False,
    min_len: int = 500,
    max_len: int = 50000,
    lognormal_mu: float = 8.5,
    lognormal_sigma: float = 1.2,
    verbose: bool = False,
) -> Tuple[List[SeqRecord], List[Tuple]]:
    """
    Add exactly `needed` contaminant (prokaryotic/eukaryotic) contigs.
    """
    if needed <= 0:
        return records, ground_truth

    if verbose:
        console.log(f"[cyan]Adding {needed} {category} sequences...[/]")

    added = 0

    while added < needed:
        if not use_synthetic and source_records:
            # Use real genomes
            rec = random.choice(source_records)
            num_frags = random.randint(1, 3)
            if verbose:
                console.log(
                    f"Fragmenting {rec.id} into {num_frags} fragments for {category}"
                )
            frags = fragment_sequence(
                rec,
                min_len=min_len,
                max_len=max_len,
                num_fragments=num_frags,
                lognormal_mu=lognormal_mu,
                lognormal_sigma=lognormal_sigma,
            )
            for frag in frags:
                new_id = f"{rec.id}_{category}_{len(records)+1}"
                frag.id = new_id
                frag.description = f"source={rec.id}"
                records.append(frag)
                meta = (
                    new_id,
                    rec.id,  # cluster = source genome
                    category,
                    rec.id,
                    len(frag.seq),
                    "no",
                    "NA",
                    1,
                    "no",
                )
                ground_truth.append(meta)
                added += 1
                if verbose:
                    console.log(f"Added contig {new_id} (length {len(frag.seq)})")
                if added >= needed:
                    break
        else:
            # Synthetic sequences
            length = random.randint(min_len, max_len)
            seq = "".join(random.choices("ACGT", k=length))
            new_id = f"{category}_synthetic_{len(records)+1}"
            rec = SeqRecord(Seq(seq), id=new_id, description=f"length={length}")
            records.append(rec)
            meta = (
                new_id,
                new_id,  # cluster = self
                category,
                "synthetic",
                length,
                "no",
                "NA",
                1,
                "yes",
            )
            ground_truth.append(meta)
            added += 1
            if verbose:
                console.log(f"Added synthetic contig {new_id} (length {length})")

    return records, ground_truth


def duplicate_records(
    records: List[SeqRecord],
    ground_truth: List[Tuple],
    duplication_factor: int,
    verbose: bool = False,
) -> Tuple[List[SeqRecord], List[Tuple]]:
    """
    Duplicate each record `duplication_factor` times.
    Updates seq_id and metadata with copy number.
    """
    if duplication_factor <= 1:
        return records, ground_truth

    console.print(f"[yellow]Duplicating fragments {duplication_factor} times...[/]")
    new_records = []
    new_gt = []

    for rec, meta in zip(records, ground_truth):
        (
            seq_id,
            cluster,
            category,
            source_id,
            length,
            is_strain,
            mut_rate,
            dup_copy,
            is_synthetic,
        ) = meta
        for copy in range(1, duplication_factor + 1):
            new_rec = rec[:]
            new_id = f"{seq_id}_dup{copy}"
            new_rec.id = new_id
            new_rec.description = rec.description + f" copy={copy}"
            new_records.append(new_rec)
            new_meta = (
                new_id,
                cluster,
                category,
                source_id,
                length,
                is_strain,
                mut_rate,
                copy,
                is_synthetic,
            )
            new_gt.append(new_meta)
            if verbose:
                console.log(f"Duplicated {seq_id} -> {new_id} (copy {copy})")

    return new_records, new_gt


def write_ground_truth(gt_file: str, ground_truth: List[Tuple]) -> None:
    """Write ground truth TSV with extensive metadata."""
    with open(gt_file, "w") as f:
        f.write(
            "sequence_id\ttrue_cluster\tsource_category\tsource_genome\t"
            "length\tstrain_mutated\tmutation_rate\tduplication_copy\tis_synthetic\n"
        )
        for meta in ground_truth:
            f.write(
                f"{meta[0]}\t{meta[1]}\t{meta[2]}\t{meta[3]}\t"
                f"{meta[4]}\t{meta[5]}\t{meta[6]}\t{meta[7]}\t{meta[8]}\n"
            )


# ----------------------------------------------------------------------
# Main generator
# ----------------------------------------------------------------------
def generate_mock_dataset(
    name: str,
    total_seqs: int,
    outdir: str,
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
    num_species: int = 10,
    seed: int = 42,
    lognormal_mu: float = 8.5,
    lognormal_sigma: float = 1.2,
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
        console.print(
            f"[yellow][SKIP] {name} already exists. Use --force to regenerate.[/]"
        )
        if not no_ground_truth and not os.path.exists(gt_file):
            console.print(
                "[red][WARN] Ground truth missing; please regenerate with --force.[/]"
            )
        return

    if dry_run:
        console.print(
            "[cyan]DRY-RUN: Would generate dataset with the following parameters:[/]"
        )
        console.print(f"Name: {name}, total sequences: {total_seqs}")
        console.print(
            f"Virus fraction: {virus_frac}, Prokaryote: {prokaryote_frac}, Eukaryote: {eukaryote_frac}"
        )
        console.print(f"Strain mode: {strain_mode}, duplications: {duplication_factor}")
        console.print(f"Lognormal mu: {lognormal_mu}, sigma: {lognormal_sigma}")
        console.print(f"Output: {outfile}, ground truth: {gt_file}")
        return

    console.print(f"[green][GEN] Generating {name} with {total_seqs} sequences...[/]")
    if verbose:
        console.log(
            f"Using lognormal_mu={lognormal_mu}, lognormal_sigma={lognormal_sigma}"
        )

    # Load sequences
    viral_records = read_fasta_seq(viral_seq)
    if not viral_records:
        raise RuntimeError("Viral sequences file is empty or missing.")
    prokaryotic_records = read_fasta_seq(prokaryotic_seq) if prokaryotic_seq else []
    eukaryotic_records = read_fasta_seq(eukaryotic_seq) if eukaryotic_seq else []

    if verbose:
        console.log(f"Loaded {len(viral_records)} viral records")
        console.log(f"Loaded {len(prokaryotic_records)} prokaryotic records")
        console.log(f"Loaded {len(eukaryotic_records)} eukaryotic records")

    # Target counts
    n_virus = int(total_seqs * virus_frac)
    n_prok = int(total_seqs * prokaryote_frac)
    n_euk = total_seqs - n_virus - n_prok

    if verbose:
        console.log(f"Target viral sequences: {n_virus}")
        console.log(f"Target prokaryotic sequences: {n_prok}")
        console.log(f"Target eukaryotic sequences: {n_euk}")

    records = []
    ground_truth = []

    # ===== Viral sequences =====
    if strain_mode:
        if len(viral_records) < num_species:
            raise ValueError(
                f"Not enough viral genomes for strain mode (need {num_species}, have {len(viral_records)})."
            )
        selected = random.sample(viral_records, num_species)
        per_species = n_virus // num_species
        remainder = n_virus % num_species

        console.print(
            f"[cyan]Strain mode: {num_species} species, {per_species} strains each (+ {remainder} extra)[/]"
        )
        if verbose:
            console.log(f"Selected {len(selected)} species for strain generation")

        for species_idx, species_seq in enumerate(selected):
            source_id = species_seq.id
            strains_needed = per_species + (1 if species_idx < remainder else 0)
            if strains_needed <= 0:
                continue
            if verbose:
                console.log(f"Species {source_id}: generating {strains_needed} strains")

            strain_tuples = generate_strains(
                species_seq,
                source_id,
                strains_needed,
                mut_rate_min,
                mut_rate_max,
            )

            for strain_rec, cluster_label, mut_rate in strain_tuples:
                num_frags = random.randint(1, 5)
                if verbose:
                    console.log(
                        f"Strain {strain_rec.id}: fragmenting into {num_frags} contigs (mut_rate={mut_rate:.4f})"
                    )
                frags = fragment_sequence(
                    strain_rec,
                    num_fragments=num_frags,
                    lognormal_mu=lognormal_mu,
                    lognormal_sigma=lognormal_sigma,
                )
                for frag in frags:
                    new_id = f"{strain_rec.id}_frag_{len(records)+1}"
                    frag.id = new_id
                    frag.description = (
                        f"strain={strain_rec.id} source={source_id} "
                        f"ANI={strain_rec.description.split('ANI=')[-1].split()[0]}"
                    )
                    records.append(frag)
                    meta = (
                        new_id,
                        cluster_label,
                        "viral",
                        source_id,
                        len(frag.seq),
                        "yes",
                        f"{mut_rate:.6f}",
                        1,
                        "no",
                    )
                    ground_truth.append(meta)
                    if len(records) >= n_virus:
                        break
                if len(records) >= n_virus:
                    break
            if len(records) >= n_virus:
                break

        # Fill remaining viral count with non‑strain fragments
        while len(records) < n_virus:
            rec = random.choice(viral_records)
            if verbose:
                console.log(f"Filling with non-strain fragment from {rec.id}")
            frags = fragment_sequence(
                rec,
                num_fragments=1,
                lognormal_mu=lognormal_mu,
                lognormal_sigma=lognormal_sigma,
            )
            for frag in frags:
                new_id = f"{rec.id}_frag_{len(records)+1}"
                frag.id = new_id
                frag.description = f"source={rec.id}"
                records.append(frag)
                meta = (
                    new_id,
                    rec.id,
                    "viral",
                    rec.id,
                    len(frag.seq),
                    "no",
                    "NA",
                    1,
                    "no",
                )
                ground_truth.append(meta)
                if len(records) >= n_virus:
                    break
    else:
        console.print("[cyan]Normal mode: fragmenting viral genomes into contigs[/]")
        while len(records) < n_virus:
            rec = random.choice(viral_records)
            num_frags = random.randint(1, 5)
            if verbose:
                console.log(f"Fragmenting {rec.id} into {num_frags} contigs")
            frags = fragment_sequence(
                rec,
                num_fragments=num_frags,
                lognormal_mu=lognormal_mu,
                lognormal_sigma=lognormal_sigma,
            )
            for frag in frags:
                new_id = f"{rec.id}_frag_{len(records)+1}"
                frag.id = new_id
                frag.description = f"source={rec.id}"
                records.append(frag)
                meta = (
                    new_id,
                    rec.id,
                    "viral",
                    rec.id,
                    len(frag.seq),
                    "no",
                    "NA",
                    1,
                    "no",
                )
                ground_truth.append(meta)
                if len(records) >= n_virus:
                    break
            if len(records) >= n_virus:
                break

    # Trim to exactly n_virus (if overshot)
    if len(records) > n_virus:
        if verbose:
            console.log(f"Trimming viral records from {len(records)} to {n_virus}")
        records = records[:n_virus]
        ground_truth = ground_truth[:n_virus]

    if verbose:
        console.log(f"After viral generation: {len(records)} sequences")

    # ===== Prokaryotic sequences =====
    use_synthetic_prok = not bool(prokaryotic_records)
    records, ground_truth = add_contaminant_contigs(
        records,
        ground_truth,
        prokaryotic_records,
        n_prok,
        label_prefix="prok",
        category="prokaryotic",
        use_synthetic=use_synthetic_prok,
        min_len=500,
        max_len=50000,
        lognormal_mu=lognormal_mu,
        lognormal_sigma=lognormal_sigma,
        verbose=verbose,
    )

    # ===== Eukaryotic sequences =====
    use_synthetic_euk = not bool(eukaryotic_records)
    records, ground_truth = add_contaminant_contigs(
        records,
        ground_truth,
        eukaryotic_records,
        n_euk,
        label_prefix="euk",
        category="eukaryotic",
        use_synthetic=use_synthetic_euk,
        min_len=500,
        max_len=50000,
        lognormal_mu=lognormal_mu,
        lognormal_sigma=lognormal_sigma,
        verbose=verbose,
    )

    if verbose:
        console.log(f"After adding contaminants: {len(records)} sequences")

    # ===== Duplication =====
    records, ground_truth = duplicate_records(
        records, ground_truth, duplication_factor, verbose=verbose
    )

    # Shuffle records (keep ground truth in sync)
    combined = list(zip(records, ground_truth))
    random.shuffle(combined)
    records, ground_truth = zip(*combined) if combined else ([], [])

    if verbose:
        console.log(f"Final sequence count: {len(records)}")

    # Write FASTA
    os.makedirs(outdir, exist_ok=True)
    SeqIO.write(records, outfile, "fasta")
    console.print(f"[green][DONE] Wrote {len(records)} sequences to {outfile}[/]")

    # Write ground truth TSV
    if not no_ground_truth:
        write_ground_truth(gt_file, ground_truth)
        console.print(f"[green][DONE] Wrote ground truth to {gt_file}[/]")


# ----------------------------------------------------------------------
# Command-line interface
# ----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Generate a mock contig FASTA dataset with ground truth.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--name", required=True, help="Dataset name (e.g., Mock-10K)")
    parser.add_argument(
        "--num-sequences", type=int, required=True, help="Total number of sequences"
    )
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument(
        "--viral-seq", required=True, help="FASTA file with viral genomes"
    )
    parser.add_argument(
        "--prokaryotic-seq",
        default="",
        help="FASTA file with prokaryotic genomes",
    )
    parser.add_argument(
        "--eukaryotic-seq",
        default="",
        help="FASTA file with eukaryotic genomes",
    )
    parser.add_argument(
        "--virus-frac", type=float, default=0.5, help="Fraction of viral sequences"
    )
    parser.add_argument(
        "--prokaryote-frac",
        type=float,
        default=0.3,
        help="Fraction of prokaryotic sequences",
    )
    parser.add_argument(
        "--eukaryote-frac",
        type=float,
        default=0.2,
        help="Fraction of eukaryotic sequences",
    )
    parser.add_argument(
        "--strain-mode",
        type=int,
        default=0,
        help="Enable strain variation (1) or not (0)",
    )
    parser.add_argument(
        "--mut-rate-min",
        type=float,
        default=0.001,
        help="Minimum mutation rate (0.001 = 99.9%% ANI)",
    )
    parser.add_argument(
        "--mut-rate-max",
        type=float,
        default=0.05,
        help="Maximum mutation rate (0.05 = 95%% ANI)",
    )
    parser.add_argument(
        "--duplication-factor",
        type=int,
        default=1,
        help="Number of copies of each fragment",
    )
    parser.add_argument(
        "--num-species", type=int, default=10, help="Number of species for strain mode"
    )
    parser.add_argument(
        "--lognormal-mu",
        type=float,
        default=8.5,
        help="Mean of the log-length distribution (median fragment length = exp(mu))",
    )
    parser.add_argument(
        "--lognormal-sigma",
        type=float,
        default=1.2,
        help="Spread (standard deviation) of the log-length distribution",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument(
        "--ground-truth", help="Output path for ground truth TSV (default: auto)"
    )
    parser.add_argument(
        "--no-ground-truth", action="store_true", help="Do not write ground truth file"
    )
    parser.add_argument(
        "--force", action="store_true", help="Regenerate even if output exists"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Enable verbose logging"
    )
    parser.add_argument(
        "--dry-run",
        "-n",
        action="store_true",
        help="Preview actions without writing files",
    )
    args = parser.parse_args()

    # Validate fractions
    total_frac = args.virus_frac + args.prokaryote_frac + args.eukaryote_frac
    if abs(total_frac - 1.0) > 1e-6:
        raise ValueError(f"Fractions must sum to 1.0 (got {total_frac})")

    # Validate input files
    if not args.dry_run:
        if not os.path.exists(args.viral_seq):
            raise FileNotFoundError(f"Viral sequence file not found: {args.viral_seq}")
        if args.prokaryotic_seq and not os.path.exists(args.prokaryotic_seq):
            raise FileNotFoundError(
                f"Prokaryotic sequence file not found: {args.prokaryotic_seq}"
            )
        if args.eukaryotic_seq and not os.path.exists(args.eukaryotic_seq):
            raise FileNotFoundError(
                f"Eukaryotic sequence file not found: {args.eukaryotic_seq}"
            )

    generate_mock_dataset(
        name=args.name,
        total_seqs=args.num_sequences,
        outdir=args.outdir,
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
        num_species=args.num_species,
        seed=args.seed,
        lognormal_mu=args.lognormal_mu,
        lognormal_sigma=args.lognormal_sigma,
        ground_truth_file=args.ground_truth,
        no_ground_truth=args.no_ground_truth,
        force=args.force,
        verbose=args.verbose,
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    main()
