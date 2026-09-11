#!/usr/bin/env python3
"""
generate_mock_clust_data.py - Generate a mock FASTA dataset for clustering benchmarking.

Approach:
- For each genome (viral, prokaryotic, eukaryotic), generate random non-overlapping
  fragments with lengths drawn from a log-normal distribution.
- Each fragment becomes a "set" — a cluster that clustering algorithms should recover.
- For each set, generate multiple mutated copies with controlled point mutations
  and indels. The number of copies can follow a uniform or Poisson distribution.
- All categories (viral, prokaryotic, eukaryotic) are mutated in the same way.
- Balanced representation: sequences are distributed across genomes based on genome
  size or a fixed target, ensuring no single genome dominates.
- Exact target counts are enforced at the SET level, so no set is ever reduced to a
  singleton.

Ground truth:
- true_cluster = set_id (i.e., {source_genome}_frag_{fragment_index})
  This is what a sequence-similarity-based clustering algorithm can actually recover.
- Mutation metadata is stored in compact form (counts + truncated position/type lists)
  to avoid extremely long lines that break downstream TSV parsers.

Inputs:
  --viral-seq          FASTA with viral genomes (required)
  --prokaryotic-seq    FASTA with prokaryotic genomes (optional; synthetic fallback)
  --eukaryotic-seq     FASTA with eukaryotic genomes (optional; synthetic fallback)

Fractions (must sum to 1.0):
  --virus-frac         fraction of viral sequences (default 0.5)
  --prokaryote-frac    fraction of prokaryotic sequences (default 0.3)
  --eukaryote-frac     fraction of eukaryotic sequences (default 0.2)

Fragmentation:
  --fragments-per-10kb     number of fragments per 10 kb of genome (default 1)
  --min-fragments-per-genome  minimum fragments per genome (default 10)
  --max-fragments-per-genome  maximum fragments per genome (default 500)

Copies per set:
  --copies-mean        mean number of copies per set (default 3)
  --copies-min         minimum copies per set (default 2)
  --copies-max         maximum copies per set (default 5)
  --copies-distribution  distribution for copies: 'uniform' or 'poisson' (default 'uniform')

Mutation:
  --mut-rate-min       minimum mutation rate (default 0.0001, 99.9% ANI)
  --mut-rate-max       maximum mutation rate (default 0.15, 70% ANI)
  --indel-rate         fraction of mutations that are indels (default 0.1)

Other:
  --total-sequences    total sequences to generate (across all categories)
  --seed               random seed (default 42)
  --outdir             output directory
  --force              overwrite existing files
  --verbose            detailed logging
  --dry-run            preview settings without writing
"""

import os
import sys
import math
import random
import argparse
import time
from collections import defaultdict
from typing import List, Tuple, Optional, Dict, Any
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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
# Helper functions
# ----------------------------------------------------------------------


def poisson(lam: float) -> int:
    """
    Generate a Poisson random variate with mean `lam` using the exponential
    interarrival method. This is a fallback for environments where
    random.poisson is not available (though it usually is).
    """
    if lam <= 0:
        return 0
    n = 0
    total = 0.0
    while total < 1.0:
        total += random.expovariate(lam)
        n += 1
    return n - 1


def read_fasta_seq(seq_file: str) -> List[SeqRecord]:
    """Return a list of SeqRecord objects from a FASTA file."""
    if not os.path.exists(seq_file) or os.path.getsize(seq_file) == 0:
        return []
    return list(SeqIO.parse(seq_file, "fasta"))


def generate_fragment_length(
    lognormal_mu: float = 8.5,
    lognormal_sigma: float = 1.2,
    min_len: int = 500,
    max_len: int = 50000,
) -> int:
    """Generate a fragment length from a truncated log-normal distribution."""
    for _ in range(100):
        z = random.gauss(0, 1)
        candidate = int(round(math.exp(lognormal_mu + lognormal_sigma * z)))
        if min_len <= candidate <= max_len:
            return candidate
    return random.randint(min_len, max_len)


def extract_random_fragments(
    seq_record: SeqRecord,
    n_fragments: int,
    min_len: int = 500,
    max_len: int = 50000,
    lognormal_mu: float = 8.5,
    lognormal_sigma: float = 1.2,
) -> List[SeqRecord]:
    """
    Extract n_fragments non-overlapping random fragments from a sequence.
    Fragments are placed sequentially to avoid overlap.
    """
    seq_len = len(seq_record.seq)
    if seq_len < min_len:
        return [seq_record] if n_fragments > 0 else []

    # Generate lengths
    lengths = []
    for _ in range(n_fragments):
        length = generate_fragment_length(
            lognormal_mu, lognormal_sigma, min_len, max_len
        )
        lengths.append(min(length, seq_len))

    # Scale lengths to fit the genome if total exceeds seq_len
    total_len = sum(lengths)
    if total_len > seq_len:
        scale = seq_len / total_len
        lengths = [max(min_len, int(l * scale)) for l in lengths]
        lengths = [max(min_len, l) for l in lengths]

    # Place fragments sequentially (non-overlapping)
    fragments = []
    current_pos = 0
    for i, length in enumerate(lengths):
        if current_pos >= seq_len:
            break
        end = min(current_pos + length, seq_len)
        if end - current_pos >= min_len:
            fragments.append(seq_record[current_pos:end])
        current_pos = end

    # If we didn't get all fragments, fill with random intervals
    while len(fragments) < n_fragments:
        start = random.randint(0, max(0, seq_len - min_len))
        end = min(start + random.randint(min_len, max_len), seq_len)
        if end - start >= min_len:
            fragments.append(seq_record[start:end])

    return fragments


def mutate_sequence(
    sequence: Seq,
    mutation_rate: float,
    indel_rate: float = 0.1,
) -> Tuple[Seq, List[Tuple[int, str, Optional[str]]]]:
    """
    Introduce mutations into a sequence.

    Returns:
        (mutated_sequence, mutations)
    mutations: list of (position, mutation_type, detail)
    mutation_type: 'substitution', 'insertion', 'deletion'
    """
    seq_str = list(str(sequence))
    seq_len = len(seq_str)
    mutations = []

    # Number of mutations: Poisson(mutation_rate * seq_len)
    mean_muts = mutation_rate * seq_len
    if mean_muts > 0:
        n_mutations = max(1, poisson(mean_muts))
    else:
        n_mutations = 0

    for _ in range(n_mutations):
        if len(seq_str) <= 1:
            break
        pos = random.randint(0, len(seq_str) - 1)

        is_indel = random.random() < indel_rate

        if is_indel:
            if len(seq_str) > 1 and random.random() < 0.5:
                # Deletion
                removed = seq_str.pop(pos)
                mutations.append((pos, "deletion", removed))
            else:
                # Insertion
                new_base = random.choice("ACGT")
                seq_str.insert(pos, new_base)
                mutations.append((pos, "insertion", new_base))
        else:
            # Substitution
            old_base = seq_str[pos]
            new_base = random.choice([b for b in "ACGT" if b != old_base])
            seq_str[pos] = new_base
            mutations.append((pos, "substitution", f"{old_base}->{new_base}"))

    return Seq("".join(seq_str)), mutations


def summarize_mutations(
    mutations: List[Tuple[int, str, Optional[str]]],
    max_positions: int = 20,
) -> Tuple[str, str, str]:
    """
    Produce compact summaries of a mutation list.

    Returns:
        (mutation_summary, mutation_positions_short, mutation_types_short)

    - mutation_summary:      "S:12,I:3,D:2" (counts per type)
    - mutation_positions_short: first `max_positions` positions, comma-joined
    - mutation_types_short:  first `max_positions` types, comma-joined
    """
    counts = {"substitution": 0, "insertion": 0, "deletion": 0}
    for _, mtype, _ in mutations:
        if mtype in counts:
            counts[mtype] += 1

    summary = (
        f"S:{counts['substitution']},I:{counts['insertion']},D:{counts['deletion']}"
    )

    positions = [str(m[0]) for m in mutations[:max_positions]]
    types = [m[1] for m in mutations[:max_positions]]

    positions_short = ",".join(positions)
    types_short = ",".join(types)

    if len(mutations) > max_positions:
        positions_short += f",...(+{len(mutations) - max_positions} more)"
        types_short += f",...(+{len(mutations) - max_positions} more)"

    return summary, positions_short, types_short


def create_mutation_set(
    source_fragment: SeqRecord,
    source_genome: str,
    fragment_idx: int,
    category: str,
    num_copies: int,
    mut_rate_min: float = 0.001,
    mut_rate_max: float = 0.05,
    indel_rate: float = 0.1,
    verbose: bool = False,
) -> List[Tuple[SeqRecord, Dict[str, Any]]]:
    """
    Create a set of mutated copies from a source fragment.
    set_id = {source_genome}_frag_{fragment_idx}
    Mutation metadata is stored compactly (counts + truncated lists).
    """
    results = []
    set_id = f"{source_genome}_frag_{fragment_idx}"

    # Original copy (no mutations)
    orig_rec = source_fragment[:]
    orig_id = f"{source_genome}_frag_{fragment_idx}_copy_1"
    orig_rec.id = orig_id
    orig_rec.description = (
        f"source={source_genome} frag={fragment_idx} set={set_id} copy=1 mut_rate=0.000"
    )
    results.append(
        (
            orig_rec,
            {
                "set_id": set_id,
                "copy_number": 1,
                "mutation_rate": 0.0,
                "mutation_summary": "S:0,I:0,D:0",
                "mutation_positions": "",
                "mutation_types": "",
                "is_duplicate": False,
            },
        )
    )

    if verbose:
        console.log(f"      set {set_id} → {num_copies} copies")

    # Mutated copies (2..num_copies)
    for copy_idx in range(2, num_copies + 1):
        mut_rate = random.uniform(mut_rate_min, mut_rate_max)
        mutated_seq, mutations = mutate_sequence(
            source_fragment.seq, mut_rate, indel_rate
        )

        summary, pos_short, types_short = summarize_mutations(mutations)

        mut_rec = SeqRecord(
            mutated_seq,
            id=f"{source_genome}_frag_{fragment_idx}_copy_{copy_idx}",
            description=(
                f"source={source_genome} frag={fragment_idx} set={set_id} "
                f"copy={copy_idx} mut_rate={mut_rate:.4f}"
            ),
        )
        results.append(
            (
                mut_rec,
                {
                    "set_id": set_id,
                    "copy_number": copy_idx,
                    "mutation_rate": mut_rate,
                    "mutation_summary": summary,
                    "mutation_positions": pos_short,
                    "mutation_types": types_short,
                    "is_duplicate": False,
                },
            )
        )

    return results


def generate_balanced_category(
    source_records: List[SeqRecord],
    category: str,
    target_sequences: int,
    fragments_per_10kb: float = 1.0,
    min_fragments_per_genome: int = 10,
    max_fragments_per_genome: int = 500,
    copies_mean: float = 3.0,
    copies_min: int = 2,
    copies_max: int = 5,
    copies_distribution: str = "uniform",
    mut_rate_min: float = 0.001,
    mut_rate_max: float = 0.05,
    indel_rate: float = 0.1,
    lognormal_mu: float = 8.5,
    lognormal_sigma: float = 1.2,
    min_len: int = 500,
    max_len: int = 50000,
    use_synthetic: bool = False,
    verbose: bool = False,
) -> Tuple[List[SeqRecord], List[Dict]]:
    """
    Generate a balanced set of fragments from a list of source genomes.

    Returns:
        (records, ground_truth)
    """
    if target_sequences <= 0:
        return [], []

    # ============================================================
    # Synthetic mode: singleton sets (each is its own cluster)
    # ============================================================
    if use_synthetic or not source_records:
        if verbose:
            console.log(
                f"[cyan]Generating {target_sequences} synthetic {category} sequences (singletons)[/]"
            )
        records = []
        gt = []
        for i in range(target_sequences):
            length = random.randint(min_len, max_len)
            seq = "".join(random.choices("ACGT", k=length))
            seq_id = f"{category}_synthetic_{i+1}"
            rec = SeqRecord(Seq(seq), id=seq_id, description=f"length={length}")
            records.append(rec)
            gt.append(
                {
                    "sequence_id": seq_id,
                    "true_cluster": seq_id,
                    "set_id": seq_id,
                    "source_genome": "synthetic",
                    "source_fragment": "1",
                    "length": length,
                    "mutation_rate": 0.0,
                    "mutation_summary": "S:0,I:0,D:0",
                    "mutation_positions": "",
                    "mutation_types": "",
                    "copy_number": 1,
                    "is_duplicate": False,
                    "source_category": category,
                    "is_synthetic": "yes",
                }
            )
            if verbose and (i + 1) % 1000 == 0:
                console.log(f"  Synthetic {i+1}/{target_sequences}")
        return records, gt

    # ============================================================
    # Real genomes mode
    # ============================================================
    n_genomes = len(source_records)
    if n_genomes == 0:
        return [], []

    if verbose:
        console.log(
            f"[cyan]Generating {target_sequences} sequences from {n_genomes} {category} genomes[/]"
        )

    # Determine average copies per set
    if copies_distribution == "poisson":
        avg_copies = copies_mean
    else:
        avg_copies = (copies_min + copies_max) / 2
    if avg_copies <= 0:
        avg_copies = 1

    total_fragments_needed = (
        int(target_sequences / avg_copies) if avg_copies > 0 else target_sequences
    )
    total_fragments_needed = max(
        total_fragments_needed, n_genomes * min_fragments_per_genome
    )

    # Compute fragments per genome based on genome length
    total_genome_length = sum(len(g.seq) for g in source_records)
    if total_genome_length == 0:
        total_genome_length = 1

    fragments_per_genome = {}
    for genome in source_records:
        len_based = max(1, int(len(genome.seq) / 10000 * fragments_per_10kb))
        len_based = max(
            min_fragments_per_genome, min(len_based, max_fragments_per_genome)
        )
        fragments_per_genome[genome.id] = len_based

    # Scale to match total_fragments_needed
    current_total = sum(fragments_per_genome.values())
    if current_total < total_fragments_needed:
        scale = total_fragments_needed / current_total
        for genome in source_records:
            new_val = int(fragments_per_genome[genome.id] * scale)
            new_val = max(
                min_fragments_per_genome, min(new_val, max_fragments_per_genome)
            )
            fragments_per_genome[genome.id] = new_val

    # Distribute remaining
    remaining = total_fragments_needed - sum(fragments_per_genome.values())
    if remaining > 0:
        sorted_ids = sorted(
            fragments_per_genome.keys(),
            key=lambda x: fragments_per_genome[x],
            reverse=True,
        )
        for i in range(remaining):
            genome_id = sorted_ids[i % len(sorted_ids)]
            if fragments_per_genome[genome_id] < max_fragments_per_genome:
                fragments_per_genome[genome_id] += 1

    if verbose:
        console.log(
            f"  Total fragments needed: {total_fragments_needed}, "
            f"actual: {sum(fragments_per_genome.values())}"
        )

    # ============================================================
    # Generate fragments and mutation sets
    # ============================================================
    records = []
    ground_truth = []

    for idx, genome in enumerate(source_records):
        n_frags = fragments_per_genome[genome.id]
        if verbose:
            console.log(
                f"  Processing genome {idx+1}/{n_genomes}: {genome.id} "
                f"(length {len(genome.seq)}) → {n_frags} fragments"
            )

        frags = extract_random_fragments(
            genome,
            n_frags,
            min_len=min_len,
            max_len=max_len,
            lognormal_mu=lognormal_mu,
            lognormal_sigma=lognormal_sigma,
        )

        for frag_idx, frag in enumerate(frags):
            if copies_distribution == "poisson":
                num_copies = max(1, int(round(poisson(copies_mean))))
                num_copies = max(copies_min, min(num_copies, copies_max))
            else:
                num_copies = random.randint(copies_min, copies_max)

            mutation_set = create_mutation_set(
                frag,
                genome.id,
                frag_idx + 1,
                category,
                num_copies,
                mut_rate_min=mut_rate_min,
                mut_rate_max=mut_rate_max,
                indel_rate=indel_rate,
                verbose=verbose,
            )

            for mut_record, meta in mutation_set:
                records.append(mut_record)
                gt_entry = {
                    "sequence_id": mut_record.id,
                    "true_cluster": meta["set_id"],  # fragment-level (attainable)
                    "set_id": meta["set_id"],
                    "source_genome": genome.id,
                    "source_fragment": f"{frag_idx+1}",
                    "length": len(mut_record.seq),
                    "mutation_rate": meta["mutation_rate"],
                    "mutation_summary": meta["mutation_summary"],
                    "mutation_positions": meta["mutation_positions"],
                    "mutation_types": meta["mutation_types"],
                    "copy_number": meta["copy_number"],
                    "is_duplicate": meta["is_duplicate"],
                    "source_category": category,
                    "is_synthetic": "no",
                }
                ground_truth.append(gt_entry)

            if verbose and (len(records) % 100 == 0):
                console.log(f"    Generated {len(records)} sequences so far")

    if verbose:
        console.log(
            f"  Generated {len(records)} sequences for {category} "
            f"(target: {target_sequences})"
        )

    # ============================================================
    # Exact-count enforcement at the SET level (never breaks sets)
    # ============================================================
    def _group_by_set(recs, gts):
        groups = defaultdict(list)
        for i, g in enumerate(gts):
            groups[g["set_id"]].append(i)
        return groups

    if len(records) != target_sequences:
        if verbose:
            console.log(
                f"[yellow]Adjusting {category} sequences from {len(records)} "
                f"to {target_sequences} (set-aware)[/]"
            )

        if len(records) > target_sequences:
            # Subsample whole sets (never break a set apart)
            groups = _group_by_set(records, ground_truth)
            set_ids = list(groups.keys())
            random.shuffle(set_ids)

            keep_indices: List[int] = []
            for sid in set_ids:
                idx_list = groups[sid]
                if len(keep_indices) + len(idx_list) <= target_sequences:
                    keep_indices.extend(idx_list)
                if len(keep_indices) >= target_sequences:
                    break

            keep_indices.sort()
            records = [records[i] for i in keep_indices]
            ground_truth = [ground_truth[i] for i in keep_indices]

            if verbose:
                kept_sets = len(set(g["set_id"] for g in ground_truth))
                console.log(
                    f"[yellow]Kept {len(records)} sequences across {kept_sets} sets[/]"
                )
        else:
            # Need more: duplicate whole sets (with new set_id)
            extra_needed = target_sequences - len(records)
            if verbose:
                console.log(
                    f"[yellow]Adding {extra_needed} extra sequences by "
                    f"duplicating whole sets[/]"
                )

            groups = _group_by_set(records, ground_truth)
            set_ids = list(groups.keys())
            i = 0
            while len(records) < target_sequences and set_ids:
                sid = set_ids[i % len(set_ids)]
                idx_list = groups[sid]
                size = len(idx_list)
                if len(records) + size <= target_sequences:
                    for j in idx_list:
                        rec = records[j]
                        gt = ground_truth[j]
                        new_rec = rec[:]
                        new_id = f"{rec.id}_extraset_{i+1}"
                        new_rec.id = new_id
                        new_rec.description = rec.description + f" extraset={i+1}"
                        records.append(new_rec)
                        new_gt = gt.copy()
                        new_gt["sequence_id"] = new_id
                        new_gt["set_id"] = f"{gt['set_id']}_dup{i+1}"
                        new_gt["true_cluster"] = new_gt["set_id"]
                        ground_truth.append(new_gt)
                i += 1

    # Safety check: no singleton sets
    groups = _group_by_set(records, ground_truth)
    for sid, idx_list in groups.items():
        if len(idx_list) < copies_min:
            raise RuntimeError(
                f"Internal error: set {sid} has {len(idx_list)} records "
                f"(< copies_min={copies_min}). This should never happen."
            )

    if verbose:
        console.log(f"  Final {category} sequences: {len(records)}")

    return records, ground_truth


def write_ground_truth(gt_file: str, ground_truth: List[Dict]) -> None:
    """Write ground truth TSV with compact mutation metadata."""
    if not ground_truth:
        return
    with open(gt_file, "w") as f:
        f.write(
            "sequence_id\ttrue_cluster\tset_id\tsource_genome\tsource_fragment\t"
            "length\tmutation_rate\tmutation_summary\tmutation_positions_short\t"
            "mutation_types_short\tcopy_number\tis_duplicate\tsource_category\t"
            "is_synthetic\n"
        )
        for entry in ground_truth:
            f.write(
                f"{entry['sequence_id']}\t"
                f"{entry['true_cluster']}\t"
                f"{entry['set_id']}\t"
                f"{entry['source_genome']}\t"
                f"{entry['source_fragment']}\t"
                f"{entry['length']}\t"
                f"{entry['mutation_rate']}\t"
                f"{entry['mutation_summary']}\t"
                f"{entry['mutation_positions']}\t"
                f"{entry['mutation_types']}\t"
                f"{entry['copy_number']}\t"
                f"{str(entry['is_duplicate']).lower()}\t"
                f"{entry['source_category']}\t"
                f"{entry['is_synthetic']}\n"
            )


def print_summary_table(
    name: str,
    total_sequences: int,
    viral_stats: Dict,
    prok_stats: Dict,
    euk_stats: Dict,
) -> None:
    """Print a summary table of generated data."""
    table = Table(title=f"Summary for {name}", box=box.ROUNDED)
    table.add_column("Category", style="cyan")
    table.add_column("Genomes", justify="right")
    table.add_column("Fragments", justify="right")
    table.add_column("Sets", justify="right")
    table.add_column("Sequences", justify="right", style="green")

    for cat, stats in [
        ("Viral", viral_stats),
        ("Prokaryotic", prok_stats),
        ("Eukaryotic", euk_stats),
    ]:
        if stats["sequences"] > 0:
            table.add_row(
                cat,
                str(stats["genomes"]),
                str(stats["fragments"]),
                str(stats["sets"]),
                str(stats["sequences"]),
            )

    table.add_row("", "", "", "", "")
    table.add_row("[bold]Total[/]", "", "", "", f"[bold]{total_sequences}[/]")
    console.print(table)


# ----------------------------------------------------------------------
# Main generator
# ----------------------------------------------------------------------
def generate_mock_dataset(
    name: str,
    total_sequences: int,
    outdir: str,
    viral_seq: str,
    prokaryotic_seq: str,
    eukaryotic_seq: str,
    virus_frac: float = 0.5,
    prokaryote_frac: float = 0.3,
    eukaryote_frac: float = 0.2,
    fragments_per_10kb: float = 1.0,
    min_fragments_per_genome: int = 10,
    max_fragments_per_genome: int = 500,
    copies_mean: float = 3.0,
    copies_min: int = 2,
    copies_max: int = 5,
    copies_distribution: str = "uniform",
    mut_rate_min: float = 0.001,
    mut_rate_max: float = 0.05,
    indel_rate: float = 0.1,
    seed: int = 42,
    lognormal_mu: float = 8.5,
    lognormal_sigma: float = 1.2,
    min_len: int = 500,
    max_len: int = 50000,
    ground_truth_file: Optional[str] = None,
    no_ground_truth: bool = False,
    force: bool = False,
    verbose: bool = False,
    dry_run: bool = False,
) -> None:
    """Main generation function."""
    start_time = time.time()

    # Validate fractions
    total_frac = virus_frac + prokaryote_frac + eukaryote_frac
    if not abs(total_frac - 1.0) < 1e-9:
        raise ValueError(
            f"The fractions must sum to 1.0. Got virus={virus_frac}, "
            f"prokaryote={prokaryote_frac}, eukaryote={eukaryote_frac} "
            f"(sum={total_frac:.10f})."
        )

    # Set the seed for all random operations
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
        console.print(f"  Name: {name}, total sequences: {total_sequences}")
        console.print(
            f"  Fractions: virus={virus_frac}, prokaryote={prokaryote_frac}, "
            f"eukaryote={eukaryote_frac}"
        )
        console.print(
            f"  Fragments per 10kb: {fragments_per_10kb} "
            f"(min={min_fragments_per_genome}, max={max_fragments_per_genome})"
        )
        console.print(
            f"  Copies: distribution={copies_distribution}, mean={copies_mean}, "
            f"min={copies_min}, max={copies_max}"
        )
        console.print(
            f"  Mutation: rate={mut_rate_min}–{mut_rate_max}, indel={indel_rate}"
        )
        console.print(f"  Seed: {seed}")
        console.print(f"  Output: {outfile}, ground truth: {gt_file}")
        return

    console.log(f"[green][START] Generating {name} with {total_sequences} sequences[/]")
    console.log(f"Seed: {seed}")

    # Load sequences
    console.log("Loading viral genomes...")
    viral_records = read_fasta_seq(viral_seq)
    if not viral_records:
        raise RuntimeError("Viral sequences file is empty or missing.")
    console.log(f"  Loaded {len(viral_records)} viral records")

    console.log("Loading prokaryotic genomes...")
    prokaryotic_records = read_fasta_seq(prokaryotic_seq) if prokaryotic_seq else []
    console.log(f"  Loaded {len(prokaryotic_records)} prokaryotic records")

    console.log("Loading eukaryotic genomes...")
    eukaryotic_records = read_fasta_seq(eukaryotic_seq) if eukaryotic_seq else []
    console.log(f"  Loaded {len(eukaryotic_records)} eukaryotic records")

    # Target counts
    n_viral = int(total_sequences * virus_frac)
    n_prok = int(total_sequences * prokaryote_frac)
    n_euk = total_sequences - n_viral - n_prok

    console.log(
        f"Target sequences: viral={n_viral}, prokaryotic={n_prok}, "
        f"eukaryotic={n_euk}"
    )

    all_records = []
    all_ground_truth = []
    stats = {}

    # Viral
    if n_viral > 0:
        console.log("[cyan]Generating viral sets...[/]")
        recs, gt = generate_balanced_category(
            viral_records,
            category="viral",
            target_sequences=n_viral,
            fragments_per_10kb=fragments_per_10kb,
            min_fragments_per_genome=min_fragments_per_genome,
            max_fragments_per_genome=max_fragments_per_genome,
            copies_mean=copies_mean,
            copies_min=copies_min,
            copies_max=copies_max,
            copies_distribution=copies_distribution,
            mut_rate_min=mut_rate_min,
            mut_rate_max=mut_rate_max,
            indel_rate=indel_rate,
            lognormal_mu=lognormal_mu,
            lognormal_sigma=lognormal_sigma,
            min_len=min_len,
            max_len=max_len,
            use_synthetic=False,
            verbose=verbose,
        )
        all_records.extend(recs)
        all_ground_truth.extend(gt)
        stats["viral"] = {
            "genomes": len(viral_records),
            "fragments": len(set(g["source_fragment"] for g in gt)) if gt else 0,
            "sets": len(set(g["set_id"] for g in gt)) if gt else 0,
            "sequences": len(recs),
        }
        console.log(f"  Generated {len(recs)} viral sequences")

    # Prokaryotic
    use_synthetic_prok = not bool(prokaryotic_records)
    if n_prok > 0:
        console.log("[magenta]Generating prokaryotic sets...[/]")
        recs, gt = generate_balanced_category(
            prokaryotic_records,
            category="prokaryotic",
            target_sequences=n_prok,
            fragments_per_10kb=fragments_per_10kb,
            min_fragments_per_genome=min_fragments_per_genome,
            max_fragments_per_genome=max_fragments_per_genome,
            copies_mean=copies_mean,
            copies_min=copies_min,
            copies_max=copies_max,
            copies_distribution=copies_distribution,
            mut_rate_min=mut_rate_min,
            mut_rate_max=mut_rate_max,
            indel_rate=indel_rate,
            lognormal_mu=lognormal_mu,
            lognormal_sigma=lognormal_sigma,
            min_len=min_len,
            max_len=max_len,
            use_synthetic=use_synthetic_prok,
            verbose=verbose,
        )
        all_records.extend(recs)
        all_ground_truth.extend(gt)
        stats["prokaryotic"] = {
            "genomes": len(prokaryotic_records) if not use_synthetic_prok else 0,
            "fragments": len(set(g["source_fragment"] for g in gt)) if gt else 0,
            "sets": len(set(g["set_id"] for g in gt)) if gt else 0,
            "sequences": len(recs),
        }
        console.log(f"  Generated {len(recs)} prokaryotic sequences")

    # Eukaryotic
    use_synthetic_euk = not bool(eukaryotic_records)
    if n_euk > 0:
        console.log("[yellow]Generating eukaryotic sets...[/]")
        recs, gt = generate_balanced_category(
            eukaryotic_records,
            category="eukaryotic",
            target_sequences=n_euk,
            fragments_per_10kb=fragments_per_10kb,
            min_fragments_per_genome=min_fragments_per_genome,
            max_fragments_per_genome=max_fragments_per_genome,
            copies_mean=copies_mean,
            copies_min=copies_min,
            copies_max=copies_max,
            copies_distribution=copies_distribution,
            mut_rate_min=mut_rate_min,
            mut_rate_max=mut_rate_max,
            indel_rate=indel_rate,
            lognormal_mu=lognormal_mu,
            lognormal_sigma=lognormal_sigma,
            min_len=min_len,
            max_len=max_len,
            use_synthetic=use_synthetic_euk,
            verbose=verbose,
        )
        all_records.extend(recs)
        all_ground_truth.extend(gt)
        stats["eukaryotic"] = {
            "genomes": len(eukaryotic_records) if not use_synthetic_euk else 0,
            "fragments": len(set(g["source_fragment"] for g in gt)) if gt else 0,
            "sets": len(set(g["set_id"] for g in gt)) if gt else 0,
            "sequences": len(recs),
        }
        console.log(f"  Generated {len(recs)} eukaryotic sequences")

    # Shuffle records
    if verbose:
        console.log("Shuffling records...")
    combined = list(zip(all_records, all_ground_truth))
    random.shuffle(combined)
    all_records, all_ground_truth = zip(*combined) if combined else ([], [])

    if verbose:
        console.log(f"Final sequence count: {len(all_records)}")

    # Write FASTA
    console.log(f"Writing FASTA to {outfile}...")
    os.makedirs(outdir, exist_ok=True)
    SeqIO.write(all_records, outfile, "fasta")
    console.log(f"[green]Wrote {len(all_records)} sequences to {outfile}[/]")

    # Write ground truth TSV
    if not no_ground_truth:
        console.log(f"Writing ground truth to {gt_file}...")
        write_ground_truth(gt_file, all_ground_truth)
        console.log(f"[green]Wrote ground truth to {gt_file}[/]")

    # Print summary
    elapsed = time.time() - start_time
    print_summary_table(
        name,
        len(all_records),
        stats.get("viral", {"genomes": 0, "fragments": 0, "sets": 0, "sequences": 0}),
        stats.get(
            "prokaryotic", {"genomes": 0, "fragments": 0, "sets": 0, "sequences": 0}
        ),
        stats.get(
            "eukaryotic", {"genomes": 0, "fragments": 0, "sets": 0, "sequences": 0}
        ),
    )
    console.log(f"[green]Done in {elapsed:.1f} seconds[/]")


# ----------------------------------------------------------------------
# Command-line interface
# ----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Generate a mock contig FASTA dataset for clustering benchmarks.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--name", required=True, help="Dataset name (e.g., Mock-300K)")
    parser.add_argument(
        "--total-sequences",
        type=int,
        required=True,
        help="Total number of sequences to generate (across all categories)",
    )
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument(
        "--viral-seq", required=True, help="FASTA file with viral genomes"
    )
    parser.add_argument(
        "--prokaryotic-seq",
        default="",
        help="FASTA file with prokaryotic genomes (optional)",
    )
    parser.add_argument(
        "--eukaryotic-seq",
        default="",
        help="FASTA file with eukaryotic genomes (optional)",
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

    # Fragmentation
    parser.add_argument(
        "--fragments-per-10kb",
        type=float,
        default=1.0,
        help="Number of fragments per 10 kb of genome length",
    )
    parser.add_argument(
        "--min-fragments-per-genome",
        type=int,
        default=10,
        help="Minimum number of fragments per genome",
    )
    parser.add_argument(
        "--max-fragments-per-genome",
        type=int,
        default=500,
        help="Maximum number of fragments per genome",
    )

    # Copies per set
    parser.add_argument(
        "--copies-mean",
        type=float,
        default=3.0,
        help="Mean number of copies per set (used for Poisson distribution)",
    )
    parser.add_argument(
        "--copies-min", type=int, default=2, help="Minimum copies per set"
    )
    parser.add_argument(
        "--copies-max", type=int, default=5, help="Maximum copies per set"
    )
    parser.add_argument(
        "--copies-distribution",
        choices=["uniform", "poisson"],
        default="uniform",
        help="Distribution for number of copies per set ('uniform' or 'poisson')",
    )

    # Mutation
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
        "--indel-rate",
        type=float,
        default=0.1,
        help="Fraction of mutations that are insertions or deletions",
    )

    # Fragment length distribution
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
    parser.add_argument(
        "--min-len", type=int, default=500, help="Minimum fragment length (bp)"
    )
    parser.add_argument(
        "--max-len", type=int, default=50000, help="Maximum fragment length (bp)"
    )

    parser.add_argument(
        "--seed", type=int, default=42, help="Random seed for reproducibility"
    )
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
        total_sequences=args.total_sequences,
        outdir=args.outdir,
        viral_seq=args.viral_seq,
        prokaryotic_seq=args.prokaryotic_seq,
        eukaryotic_seq=args.eukaryotic_seq,
        virus_frac=args.virus_frac,
        prokaryote_frac=args.prokaryote_frac,
        eukaryote_frac=args.eukaryote_frac,
        fragments_per_10kb=args.fragments_per_10kb,
        min_fragments_per_genome=args.min_fragments_per_genome,
        max_fragments_per_genome=args.max_fragments_per_genome,
        copies_mean=args.copies_mean,
        copies_min=args.copies_min,
        copies_max=args.copies_max,
        copies_distribution=args.copies_distribution,
        mut_rate_min=args.mut_rate_min,
        mut_rate_max=args.mut_rate_max,
        indel_rate=args.indel_rate,
        seed=args.seed,
        lognormal_mu=args.lognormal_mu,
        lognormal_sigma=args.lognormal_sigma,
        min_len=args.min_len,
        max_len=args.max_len,
        ground_truth_file=args.ground_truth,
        no_ground_truth=args.no_ground_truth,
        force=args.force,
        verbose=args.verbose,
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    main()
