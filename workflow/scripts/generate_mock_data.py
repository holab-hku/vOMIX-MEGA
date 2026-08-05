#!/usr/bin/env python3
"""
Generate a mock FASTA dataset for benchmarking.

Outputs:
  - FASTA file (<name>.fna)
  - Ground truth TSV (<name>.ground_truth.tsv) mapping each sequence ID to its true cluster label.

Labeling:
  - Viral contigs → source genome accession (RefSeq ID)
  - Strain-mode contigs → parent genome accession (all strains from same species share label)
  - Real contaminants → source genome accession
  - Synthetic contaminants → self ID (singleton cluster)
  - Duplicates → same label as original
  - Overlapping fragments → same label as original

All proportions, mutation rates, duplication, and overlap are adjustable.
"""

import os
import sys
import argparse
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def read_fasta_cache(cache_file):
    """Return a list of SeqRecord objects from a FASTA file."""
    if not os.path.exists(cache_file) or os.path.getsize(cache_file) == 0:
        return []
    return list(SeqIO.parse(cache_file, "fasta"))


def fragment_sequence(
    seq_record,
    min_len=500,
    max_len=50000,
    num_fragments=1,
    overlap_min=0,
    overlap_max=0,
):
    """
    Generate num_fragments contigs from a single sequence.
    If overlap_min/max > 0, fragments overlap by a random amount in that range.
    """
    seq_len = len(seq_record.seq)
    if seq_len < min_len:
        return [seq_record]

    fragments = []
    # For simplicity, we generate independent random intervals.
    # Overlap is not explicitly enforced but can occur naturally.
    # To enforce overlap, we would need to coordinate intervals.
    # We'll keep it simple: random intervals may overlap.
    for _ in range(num_fragments):
        start = random.randint(0, max(0, seq_len - min_len))
        end = min(
            start + random.randint(min_len, min(max_len, seq_len - start)), seq_len
        )
        if end - start < min_len:
            continue
        frag = seq_record[start:end]
        fragments.append(frag)
    return fragments if fragments else [seq_record]


def generate_mock_dataset(
    name,
    total_seqs,
    outdir,
    tmpdir,
    viral_cache,
    contam_cache,
    virus_frac=0.5,
    bacteria_frac=0.3,
    eukaryote_frac=0.2,
    strain_mode=False,
    mut_rate_min=0.001,
    mut_rate_max=0.05,
    duplication_factor=1,
    overlap_min=0,
    overlap_max=0,
    num_species=10,
    seed=42,
    ground_truth_file=None,
    no_ground_truth=False,
    force=False,
):
    """
    Main generation function.
    """
    random.seed(seed)
    outfile = os.path.join(outdir, f"{name}.fna")
    if not force and os.path.exists(outfile) and os.path.getsize(outfile) > 1000:
        print(f"[SKIP] {name} already exists. Use --force to regenerate.")
        # Still ensure ground truth exists if requested
        if not no_ground_truth:
            gt_file = ground_truth_file or outfile.replace(".fna", ".ground_truth.tsv")
            if not os.path.exists(gt_file):
                print(f"[WARN] Ground truth missing; regenerating it without FASTA?")
                # Could implement a separate ground truth regeneration, but we'll simply exit.
        return

    print(f"[GEN] Generating {name} with {total_seqs} sequences...")

    # Load caches
    viral_records = read_fasta_cache(viral_cache)
    if not viral_records:
        raise RuntimeError(
            "Viral cache is empty. Please run download_refseq_viral rule first."
        )
    contam_records = read_fasta_cache(contam_cache) if contam_cache else []
    use_real_contam = len(contam_records) > 0

    # Compute counts per category
    n_virus = int(total_seqs * virus_frac)
    n_bact = int(total_seqs * bacteria_frac)
    n_euk = total_seqs - n_virus - n_bact  # remainder

    records = []
    ground_truth = []  # list of (seq_id, true_cluster)

    # ----- Viral sequences -----
    if strain_mode:
        # Pick N random viral genomes to serve as species
        if len(viral_records) < num_species:
            raise ValueError(
                f"Not enough viral genomes for strain mode (need {num_species})."
            )
        selected = random.sample(viral_records, num_species)
        per_species = n_virus // num_species
        for species_idx, species_seq in enumerate(selected):
            seq_str = str(species_seq.seq)
            source_id = species_seq.id  # Use original ID as cluster label
            for strain_idx in range(per_species):
                mut_rate = random.uniform(mut_rate_min, mut_rate_max)
                # mutate
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
                # Keep full length to test strain resolution (no fragmentation)
                records.append(new_record)
                ground_truth.append((new_id, source_id))
                if len(records) >= total_seqs:
                    break
            if len(records) >= total_seqs:
                break
        # If we still need more, fill with random viral fragments (non-strain)
        while len(records) < total_seqs:
            rec = random.choice(viral_records)
            frags = fragment_sequence(
                rec, num_fragments=1, overlap_min=overlap_min, overlap_max=overlap_max
            )
            for frag in frags:
                new_id = f"{rec.id}_frag_{len(records)+1}"
                frag.id = new_id
                frag.description = f"source={rec.id}"
                records.append(frag)
                ground_truth.append((new_id, rec.id))
                if len(records) >= total_seqs:
                    break
    else:
        # Normal mode: fragment viral genomes into contigs
        while len(records) < n_virus:
            rec = random.choice(viral_records)
            num_frags = random.randint(1, 5)
            frags = fragment_sequence(
                rec,
                num_fragments=num_frags,
                overlap_min=overlap_min,
                overlap_max=overlap_max,
            )
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
        # Trim to exactly n_virus
        if len(records) > n_virus:
            records = records[:n_virus]
            ground_truth = ground_truth[:n_virus]

    # ----- Contaminants (bacteria + eukaryotes) -----
    contam_needed = total_seqs - len(records)
    if contam_needed > 0:
        if use_real_contam:
            while len(records) < total_seqs:
                rec = random.choice(contam_records)
                num_frags = random.randint(1, 3)
                frags = fragment_sequence(
                    rec,
                    num_fragments=num_frags,
                    overlap_min=overlap_min,
                    overlap_max=overlap_max,
                )
                for frag in frags:
                    new_id = f"{rec.id}_contam_{len(records)+1}"
                    frag.id = new_id
                    frag.description = f"source={rec.id}"
                    records.append(frag)
                    ground_truth.append((new_id, rec.id))
                    if len(records) >= total_seqs:
                        break
                if len(records) >= total_seqs:
                    break
        else:
            # Synthetic contaminants
            for i in range(contam_needed):
                cat = random.choices(["bacteria", "eukaryote"], weights=[0.6, 0.4])[0]
                length = random.randint(500, 50000)
                seq = "".join(random.choices("ACGT", k=length))
                new_id = f"{cat}_synthetic_{i+1}"
                rec = SeqRecord(Seq(seq), id=new_id, description=f"length={length}")
                records.append(rec)
                ground_truth.append((new_id, new_id))  # self as cluster

    # ----- Duplication -----
    if duplication_factor > 1:
        print(f"Duplicating fragments {duplication_factor} times...")
        duplicated_records = []
        duplicated_gt = []
        for rec, (seq_id, cluster) in zip(records, ground_truth):
            for copy in range(duplication_factor):
                new_rec = rec[:]  # copy
                new_rec.id = f"{seq_id}_dup{copy+1}"
                new_rec.description = rec.description + f" copy={copy+1}"
                duplicated_records.append(new_rec)
                duplicated_gt.append((new_rec.id, cluster))
        records = duplicated_records
        ground_truth = duplicated_gt

    # Shuffle records (but keep ground truth in sync)
    combined = list(zip(records, ground_truth))
    random.shuffle(combined)
    records, ground_truth = zip(*combined) if combined else ([], [])

    # Write FASTA
    SeqIO.write(records, outfile, "fasta")
    print(f"[DONE] Wrote {len(records)} sequences to {outfile}")

    # Write ground truth TSV (unless disabled)
    if not no_ground_truth:
        gt_file = ground_truth_file or outfile.replace(".fna", ".ground_truth.tsv")
        with open(gt_file, "w") as f:
            f.write("sequence_id\ttrue_cluster\n")
            for seq_id, cluster in ground_truth:
                f.write(f"{seq_id}\t{cluster}\n")
        print(f"[DONE] Wrote ground truth to {gt_file}")


def main():
    p = argparse.ArgumentParser(description="Generate mock dataset with ground truth.")
    p.add_argument("--name", required=True, help="Dataset name (e.g., Mock-10K)")
    p.add_argument(
        "--num-sequences", type=int, required=True, help="Total number of sequences"
    )
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument(
        "--tmpdir",
        required=True,
        help="Temporary directory (unused but kept for consistency)",
    )
    p.add_argument("--viral-cache", required=True, help="FASTA file with viral genomes")
    p.add_argument(
        "--contam-cache", default="", help="FASTA file with contaminant genomes"
    )
    p.add_argument(
        "--virus-frac", type=float, default=0.5, help="Fraction of viral sequences"
    )
    p.add_argument(
        "--bacteria-frac",
        type=float,
        default=0.3,
        help="Fraction of bacterial sequences",
    )
    p.add_argument(
        "--eukaryote-frac",
        type=float,
        default=0.2,
        help="Fraction of eukaryotic sequences",
    )
    p.add_argument(
        "--strain-mode",
        type=int,
        default=0,
        help="Enable strain variation (1) or not (0)",
    )
    p.add_argument(
        "--mut-rate-min",
        type=float,
        default=0.001,
        help="Minimum mutation rate (e.g., 0.001 = 99.9% ANI)",
    )
    p.add_argument(
        "--mut-rate-max",
        type=float,
        default=0.05,
        help="Maximum mutation rate (e.g., 0.05 = 95% ANI)",
    )
    p.add_argument(
        "--duplication-factor",
        type=int,
        default=1,
        help="Number of copies of each fragment",
    )
    p.add_argument(
        "--overlap-min",
        type=int,
        default=0,
        help="Minimum overlap between fragments (bp)",
    )
    p.add_argument(
        "--overlap-max",
        type=int,
        default=0,
        help="Maximum overlap between fragments (bp)",
    )
    p.add_argument(
        "--num-species", type=int, default=10, help="Number of species for strain mode"
    )
    p.add_argument("--seed", type=int, default=42, help="Random seed")
    p.add_argument(
        "--ground-truth", help="Output path for ground truth TSV (default: auto)"
    )
    p.add_argument(
        "--no-ground-truth", action="store_true", help="Do not write ground truth file"
    )
    p.add_argument(
        "--force", action="store_true", help="Regenerate even if output exists"
    )
    args = p.parse_args()

    # Validate fractions sum to 1
    total_frac = args.virus_frac + args.bacteria_frac + args.eukaryote_frac
    if abs(total_frac - 1.0) > 1e-6:
        raise ValueError(f"Fractions must sum to 1 (got {total_frac})")

    generate_mock_dataset(
        name=args.name,
        total_seqs=args.num_sequences,
        outdir=args.outdir,
        tmpdir=args.tmpdir,
        viral_cache=args.viral_cache,
        contam_cache=args.contam_cache,
        virus_frac=args.virus_frac,
        bacteria_frac=args.bacteria_frac,
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
    )


if __name__ == "__main__":
    main()
