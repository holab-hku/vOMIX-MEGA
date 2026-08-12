#!/usr/bin/env python
"""
Aggregate assembly statistics from multiple FASTA files.

Usage:
    assembly_stats.py -i file1.fa file2.fa ... --size-dist-file dist.tsv --stat-file stats.tsv

Input files are expected to be in the format:
    assembly/{assembler}/samples/{assembly_id}/output/final.contigs.fa
The assembly_id is extracted from the path.
"""

import os
import sys
import pandas as pd
from Bio.SeqIO import parse
from argparse import ArgumentParser


def store_lengths(files):
    """
    Read each FASTA file and store contig lengths in a dictionary.
    Returns dict: {assembly_id: {contig_id: length}}
    """
    r = {}
    for f in files:
        # Extract assembly_id from path: .../samples/{assembly_id}/output/...
        # We go up two levels from the file: output/ -> assembly_id/
        # So the assembly_id is the directory name of the parent of the parent.
        parent = os.path.dirname(f)  # .../samples/{assembly_id}/output
        grandparent = os.path.dirname(parent)  # .../samples/{assembly_id}
        assembly_id = os.path.basename(grandparent)  # assembly_id

        r[assembly_id] = {}
        try:
            for record in parse(f, "fasta"):
                r[assembly_id][record.id] = len(record.seq)
        except Exception as e:
            sys.stderr.write(f"Warning: Could not parse {f}: {e}\n")
            continue

        # If the assembly has no contigs, add an empty dict
        if not r[assembly_id]:
            sys.stderr.write(f"Warning: No contigs found in {f}\n")
    return r


def size_distribute(r, thresholds=None):
    """
    Compute size distribution for each assembly.
    Returns a DataFrame with columns: assembly, min_length, num_contigs, total_length, %
    """
    if thresholds is None:
        thresholds = [
            0,
            100,
            250,
            500,
            1000,
            2500,
            5000,
            10000,
            15000,
            20000,
            25000,
            30000,
            35000,
            40000,
            45000,
            50000,
            75000,
            100000,
            125000,
            150000,
            200000,
            250000,
            500000,
        ]

    rows = []
    for name, contig_lengths in r.items():
        if not contig_lengths:
            continue
        df = pd.DataFrame.from_dict(contig_lengths, orient="index", columns=["length"])
        total_bp = df["length"].sum()

        for threshold in thresholds:
            subset = df[df["length"] >= threshold]
            if subset.empty:
                break
            n = len(subset)
            s = subset["length"].sum()
            pct = (s / total_bp) * 100 if total_bp > 0 else 0
            rows.append(
                {
                    "assembly": name,
                    "min_length": threshold,
                    "num_contigs": n,
                    "total_length": s,
                    "%": pct,
                }
            )

    return pd.DataFrame(rows)


def calculate_n_stats(df):
    """Calculate N50 and N90 from a sorted DataFrame (ascending)."""
    df_sorted = df.sort_values("length", ascending=True)
    total = df_sorted["length"].sum()
    if total == 0:
        return 0, 0
    half = total * 0.5
    tenth = total * 0.1
    cumsum = 0
    n50 = n90 = 0
    for length in df_sorted["length"]:
        cumsum += length
        if n50 == 0 and cumsum >= half:
            n50 = length
        if n90 == 0 and cumsum >= tenth:
            n90 = length
        if n50 and n90:
            break
    return n50, n90


def calculate_length_stats(contig_lengths):
    """Compute summary statistics for one assembly."""
    if not contig_lengths:
        return {
            "contigs": 0,
            "total_size_bp": 0,
            "min_length": 0,
            "max_length": 0,
            "avg_length": 0,
            "median_length": 0,
            "N50_length": 0,
            "N90_length": 0,
        }
    df = pd.DataFrame.from_dict(contig_lengths, orient="index", columns=["length"])
    total = df["length"].sum()
    n50, n90 = calculate_n_stats(df)
    return {
        "contigs": len(df),
        "total_size_bp": int(total),
        "min_length": int(df["length"].min()),
        "max_length": int(df["length"].max()),
        "avg_length": float(df["length"].mean()),
        "median_length": float(df["length"].median()),
        "N50_length": n50,
        "N90_length": n90,
    }


def generate_stat_df(r):
    """Build a DataFrame of summary statistics for all assemblies."""
    stats = {}
    for name, contig_lengths in r.items():
        stats[name] = calculate_length_stats(contig_lengths)
    index_order = [
        "contigs",
        "total_size_bp",
        "min_length",
        "max_length",
        "avg_length",
        "median_length",
        "N50_length",
        "N90_length",
    ]
    return pd.DataFrame(stats).T[index_order]


def main():
    parser = ArgumentParser(description="Aggregate assembly statistics.")
    parser.add_argument(
        "-i",
        "--infile",
        type=str,
        nargs="+",
        help="FASTA files of contigs",
        required=True,
    )
    parser.add_argument(
        "--size-dist-file", type=str, help="Write size distribution table to this file"
    )
    parser.add_argument(
        "--stat-file", type=str, help="Write summary statistics to this file"
    )
    args = parser.parse_args()

    # Store lengths
    lengths_dict = store_lengths(args.infile)

    # Generate stats
    stat_df = generate_stat_df(lengths_dict)

    if args.stat_file:
        stat_df.to_csv(args.stat_file, sep="\t", index=True)
    else:
        sys.stdout.write(stat_df.to_csv(sep="\t", index=True))

    # Size distribution
    if args.size_dist_file:
        size_df = size_distribute(lengths_dict)
        if not size_df.empty:
            size_df.to_csv(args.size_dist_file, sep="\t", index=False)
        else:
            # Write empty file with headers
            open(args.size_dist_file, "w").close()


if __name__ == "__main__":
    main()
