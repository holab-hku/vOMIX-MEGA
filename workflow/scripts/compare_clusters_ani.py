#!/usr/bin/env python3
"""
Cluster ANI evaluation script.

Takes pairwise ANI values (from FastANI or skani) and computes:
- Within-cluster ANI
- Below-threshold fraction (ANI < 0.95)
- Silhouette score
- Singleton fraction
- Davies-Bouldin index

Also reports per-sequence expected ANI (from mutation_rate) and predicted ANI.

Usage:
    python compare_clusters_ani.py \
        --ground-truth ground_truth.tsv \
        --ani-file fastani_output.tsv \
        --cluster-files cluster1.csv cluster2.csv ... \
        --sample-names exp1 exp2 ... \
        --output-dir results/ \
        --ani-tool fastani
"""

import argparse
import sys
import os
from itertools import combinations
from collections import defaultdict
import pandas as pd
import numpy as np
from sklearn.metrics import silhouette_score, davies_bouldin_score
from rich.console import Console
from rich.table import Table
from rich.progress import Progress, BarColumn, TextColumn, TimeElapsedColumn
from rich import box

console = Console()

# ------------------------------------------------------------
# Parsing functions
# ------------------------------------------------------------


def parse_cluster_file(filepath):
    """Same as before: first column = representative, second = members list."""
    df = pd.read_csv(filepath, header=None, quotechar='"', skipinitialspace=True)
    rep_keywords = ["rep", "representative", "cluster", "id"]
    set_keywords = ["set", "members", "cluster_set", "sequences", "members"]
    first_vals = df.iloc[0].astype(str).str.lower()
    is_header = False
    if len(first_vals) >= 2:
        col0 = first_vals[0]
        col1 = first_vals[1]
        if any(kw in col0 for kw in rep_keywords) and any(
            kw in col1 for kw in set_keywords
        ):
            is_header = True
    if is_header:
        df = df.iloc[1:].reset_index(drop=True)
    if df.shape[1] < 2:
        raise ValueError(f"File {filepath} has fewer than 2 columns.")
    df.columns = ["cluster_id", "cluster_set"] + [
        f"extra_{i}" for i in range(2, df.shape[1])
    ]
    records = []
    for _, row in df.iterrows():
        cluster = str(row["cluster_id"]).strip()
        if not cluster:
            continue
        members = str(row["cluster_set"]).split(",")
        for seq in members:
            seq = seq.strip()
            if seq:
                records.append({"sequence_id": seq, "cluster_id": cluster})
    return pd.DataFrame(records)


def parse_ground_truth(filepath):
    df = pd.read_csv(filepath, sep="\t")
    required = ["sequence_id", "true_cluster", "mutation_rate"]
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in ground truth.")
    # Convert mutation_rate to float, NA -> NaN
    df["mutation_rate"] = pd.to_numeric(df["mutation_rate"], errors="coerce")
    return df


def parse_ani_file(filepath, tool="fastani"):
    """
    Parse pairwise ANI TSV.
    Expected columns:
      - fastani: query, ref, ani, fragments, total_fragments, ...
      - skani: query, ref, ani, ...
    Returns a dict mapping (seq1, seq2) -> ani (float), symmetric.
    """
    df = pd.read_csv(filepath, sep="\t", header=None)
    # Determine column indices based on tool
    if tool.lower() == "fastani":
        # Typical fastani output: query, ref, ani, ...
        # We'll assume first column = query, second = ref, third = ani
        query_col, ref_col, ani_col = 0, 1, 2
    elif tool.lower() == "skani":
        # skani output: query, ref, ani, ...
        query_col, ref_col, ani_col = 0, 1, 2
    else:
        raise ValueError(f"Unknown ANI tool: {tool}. Use 'fastani' or 'skani'.")

    ani_dict = {}
    for _, row in df.iterrows():
        q = str(row[query_col]).strip()
        r = str(row[ref_col]).strip()
        ani = float(row[ani_col])
        if q == r:
            continue
        # Store both directions
        ani_dict[(q, r)] = ani
        ani_dict[(r, q)] = ani
    return ani_dict


# ------------------------------------------------------------
# Metric functions
# ------------------------------------------------------------


def within_cluster_ani(seqs, ani_dict):
    """
    Compute mean pairwise ANI among all sequences in a cluster.
    Returns mean ANI, and list of pairwise ANI values.
    """
    pairs = []
    for s1, s2 in combinations(seqs, 2):
        ani = ani_dict.get((s1, s2))
        if ani is None:
            continue
        pairs.append(ani)
    if not pairs:
        return np.nan, []
    return np.mean(pairs), pairs


def below_threshold_fraction(pair_anis, threshold=0.95):
    if not pair_anis:
        return np.nan
    below = sum(1 for a in pair_anis if a < threshold)
    return below / len(pair_anis)


def singleton_fraction(cluster_df):
    sizes = cluster_df.groupby("cluster_id").size()
    singletons = (sizes == 1).sum()
    return singletons / len(sizes) if len(sizes) > 0 else np.nan


def compute_silhouette(seqs, ani_dict):
    """
    Compute silhouette score for a set of sequences.
    Requires pairwise ANI for all pairs.
    Returns silhouette score or np.nan if missing pairs.
    """
    n = len(seqs)
    if n < 2:
        return np.nan
    # Build distance matrix (1 - ANI)
    dist_matrix = np.zeros((n, n))
    for i, s1 in enumerate(seqs):
        for j, s2 in enumerate(seqs):
            if i == j:
                dist_matrix[i, j] = 0.0
            else:
                ani = ani_dict.get((s1, s2))
                if ani is None:
                    return np.nan  # incomplete
                dist_matrix[i, j] = 1.0 - ani
    # Need cluster labels for each sequence
    # We'll use cluster assignments; but we are computing silhouette for a single cluster? No, silhouette is for all clusters together.
    # Better to compute silhouette for the entire dataset using all clusters.
    return None  # Placeholder; we'll compute globally


def compute_davies_bouldin(seqs, ani_dict):
    """Similar, requires full distance matrix."""
    # We'll implement globally as well.
    pass


# ------------------------------------------------------------
# Main analysis function
# ------------------------------------------------------------


def analyze_clusters_ani(cluster_df, ground_truth, ani_dict, sample_id=None):
    """
    For one clustering result, compute ANI-based metrics.
    Returns a dict of overall metrics and a per-sequence DataFrame.
    """
    # Merge with ground truth to get mutation_rate and expected_ANI
    merged = pd.merge(cluster_df, ground_truth, on="sequence_id", how="inner")
    if merged.empty:
        console.print(f"[yellow]No overlap with ground truth for {sample_id}[/]")
        return None, None

    # Add expected ANI: 1 - mutation_rate (if not NaN)
    merged["expected_ANI"] = 1.0 - merged["mutation_rate"]

    # Get cluster assignments
    labels = merged[["sequence_id", "cluster_id"]]

    # For each cluster, find representative (cluster_id itself)
    # Build mapping from cluster_id to representative sequence ID (assumed equal)
    # But the representative might not be a sequence in the cluster? In our CSV, rep is the cluster_id string.
    # We'll just use the cluster_id as the representative's sequence ID.
    # However, some sequences might not have the exact same ID; we'll assume cluster_id is in the sequence list.
    # If not, we'll pick the first sequence in cluster as representative.
    def get_rep(cluster_id, group):
        if cluster_id in group["sequence_id"].values:
            return cluster_id
        else:
            return group.iloc[0]["sequence_id"]

    # Compute per-sequence predicted ANI (to representative)
    merged["predicted_ANI"] = np.nan
    merged["ani_available"] = False

    for cluster_id, group in merged.groupby("cluster_id"):
        rep = get_rep(cluster_id, group)
        for idx, row in group.iterrows():
            seq = row["sequence_id"]
            if seq == rep:
                merged.loc[idx, "predicted_ANI"] = 1.0  # self
                merged.loc[idx, "ani_available"] = True
            else:
                ani = ani_dict.get((rep, seq))
                if ani is not None:
                    merged.loc[idx, "predicted_ANI"] = ani
                    merged.loc[idx, "ani_available"] = True

    # Overall metrics across all clusters
    all_seqs = merged["sequence_id"].tolist()
    # Within-cluster ANI: mean of all pair ANI within each cluster
    cluster_pair_anis = {}
    cluster_anis = []
    for cluster_id, group in merged.groupby("cluster_id"):
        seqs = group["sequence_id"].tolist()
        mean_ani, pair_list = within_cluster_ani(seqs, ani_dict)
        cluster_pair_anis[cluster_id] = pair_list
        if not np.isnan(mean_ani):
            cluster_anis.append(mean_ani)

    overall_within_ani = np.mean(cluster_anis) if cluster_anis else np.nan

    # Below-threshold fraction: count pairs <0.95 / total pairs
    all_pairs = []
    for pair_list in cluster_pair_anis.values():
        all_pairs.extend(pair_list)
    btf = below_threshold_fraction(all_pairs)

    # Singleton fraction
    sf = singleton_fraction(merged)

    # Silhouette score: need complete distance matrix for all sequences
    silhouette = np.nan
    # Build distance matrix if all pairs available
    n = len(all_seqs)
    if n >= 2:
        dist_matrix = np.zeros((n, n))
        complete = True
        for i, s1 in enumerate(all_seqs):
            for j, s2 in enumerate(all_seqs):
                if i == j:
                    dist_matrix[i, j] = 0.0
                else:
                    ani = ani_dict.get((s1, s2))
                    if ani is None:
                        complete = False
                        break
                    dist_matrix[i, j] = 1.0 - ani
            if not complete:
                break
        if complete:
            try:
                # Need cluster labels in same order as all_seqs
                label_map = {
                    seq: row["cluster_id"]
                    for seq, row in merged.set_index("sequence_id").iterrows()
                }
                labels = [label_map[seq] for seq in all_seqs]
                silhouette = silhouette_score(dist_matrix, labels, metric="precomputed")
            except Exception as e:
                if sample_id:
                    console.log(f"[yellow]Silhouette failed for {sample_id}: {e}[/]")
                silhouette = np.nan

    # Davies-Bouldin
    davies = np.nan
    if complete and n >= 2:
        try:
            label_map = {
                seq: row["cluster_id"]
                for seq, row in merged.set_index("sequence_id").iterrows()
            }
            labels = [label_map[seq] for seq in all_seqs]
            davies = davies_bouldin_score(dist_matrix, labels)
        except:
            davies = np.nan

    metrics = {
        "sample_id": sample_id,
        "n_sequences": len(merged),
        "n_clusters": len(merged["cluster_id"].unique()),
        "within_cluster_ANI_mean": overall_within_ani,
        "below_threshold_fraction": btf,
        "singleton_fraction": sf,
        "silhouette": silhouette,
        "davies_bouldin": davies,
    }

    # Per-sequence table with expected and predicted ANI
    seq_table = merged[
        ["sequence_id", "cluster_id", "expected_ANI", "predicted_ANI"]
    ].copy()
    seq_table["ani_available"] = merged["ani_available"]
    if sample_id:
        seq_table["sample_id"] = sample_id

    return metrics, seq_table


# ------------------------------------------------------------
# Main CLI
# ------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description="Evaluate clustering using ANI.")
    parser.add_argument(
        "--ground-truth", required=True, help="Ground truth TSV with mutation_rate."
    )
    parser.add_argument(
        "--ani-file", required=True, help="Pairwise ANI TSV (fastani/skani output)."
    )
    parser.add_argument(
        "--ani-tool",
        default="fastani",
        choices=["fastani", "skani"],
        help="Tool used to generate ANI file (default: fastani).",
    )
    parser.add_argument(
        "--cluster-files",
        nargs="+",
        required=True,
        help="List of clustering CSV files.",
    )
    parser.add_argument(
        "--sample-names",
        nargs="*",
        default=None,
        help="Sample IDs for each cluster file.",
    )
    parser.add_argument(
        "--output-dir", default=".", help="Directory to write output files."
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output.")
    parser.add_argument(
        "-n", "--dry-run", action="store_true", help="Dry run (no files written)."
    )
    args = parser.parse_args()

    # Validate input files
    for f in [args.ground_truth, args.ani_file] + args.cluster_files:
        if not os.path.isfile(f):
            console.print(f"[red]ERROR:[/] File not found: {f}")
            sys.exit(1)

    if not args.dry_run:
        os.makedirs(args.output_dir, exist_ok=True)

    # Parse ground truth
    console.print("[cyan]Parsing ground truth...[/]")
    gt = parse_ground_truth(args.ground_truth)

    # Parse ANI file
    console.print(f"[cyan]Parsing ANI file ({args.ani_tool})...[/]")
    ani_dict = parse_ani_file(args.ani_file, args.ani_tool)
    console.log(f"[green]✓[/] Loaded {len(ani_dict)} ANI pairs.")

    # Parse cluster files
    cluster_dfs = []
    sample_ids = []
    for i, f in enumerate(args.cluster_files):
        try:
            df = parse_cluster_file(f)
            cluster_dfs.append(df)
            sid = (
                args.sample_names[i]
                if args.sample_names and i < len(args.sample_names)
                else os.path.splitext(os.path.basename(f))[0]
            )
            sample_ids.append(sid)
        except Exception as e:
            console.print(f"[red]Error parsing {f}:[/] {e}")
            if args.verbose:
                console.print_exception()
            continue

    if not cluster_dfs:
        console.print("[red]No valid cluster files. Exiting.[/]")
        sys.exit(1)

    # Process each clustering
    all_metrics = []
    all_seq_tables = []
    for sid, df in zip(sample_ids, cluster_dfs):
        console.print(f"[cyan]Processing {sid}...[/]")
        metrics, seq_table = analyze_clusters_ani(df, gt, ani_dict, sample_id=sid)
        if metrics:
            all_metrics.append(metrics)
        if seq_table is not None and not seq_table.empty:
            all_seq_tables.append(seq_table)
            if not args.dry_run:
                out_seq = os.path.join(args.output_dir, f"{sid}_ani_per_sequence.tsv")
                seq_table.to_csv(out_seq, sep="\t", index=False)
                console.log(f"[green]✓[/] Per-sequence table written to {out_seq}")

    # Summary metrics
    if all_metrics:
        metrics_df = pd.DataFrame(all_metrics)
        if not args.dry_run:
            out_metrics = os.path.join(
                args.output_dir, "ani_cluster_metrics_summary.tsv"
            )
            metrics_df.to_csv(out_metrics, sep="\t", index=False)
            console.log(f"[green]✓[/] Cluster metrics summary written to {out_metrics}")
        # Print a nice table
        table = Table(title="ANI-based Cluster Metrics", box=box.ROUNDED)
        table.add_column("Sample")
        table.add_column("N seqs")
        table.add_column("N clusters")
        table.add_column("Within-ANI")
        table.add_column("BTF")
        table.add_column("Singleton")
        table.add_column("Silhouette")
        table.add_column("Davies-Bouldin")
        for _, row in metrics_df.iterrows():
            table.add_row(
                row["sample_id"],
                str(row["n_sequences"]),
                str(row["n_clusters"]),
                (
                    f"{row['within_cluster_ANI_mean']:.4f}"
                    if not pd.isna(row["within_cluster_ANI_mean"])
                    else "NA"
                ),
                (
                    f"{row['below_threshold_fraction']:.4f}"
                    if not pd.isna(row["below_threshold_fraction"])
                    else "NA"
                ),
                (
                    f"{row['singleton_fraction']:.4f}"
                    if not pd.isna(row["singleton_fraction"])
                    else "NA"
                ),
                f"{row['silhouette']:.4f}" if not pd.isna(row["silhouette"]) else "NA",
                (
                    f"{row['davies_bouldin']:.4f}"
                    if not pd.isna(row["davies_bouldin"])
                    else "NA"
                ),
            )
        console.print(table)
    else:
        console.print("[yellow]No metrics computed.[/]")

    # Combined per-sequence table
    if all_seq_tables:
        combined_seq = pd.concat(all_seq_tables, ignore_index=True)
        if not args.dry_run:
            out_comb = os.path.join(args.output_dir, "all_ani_per_sequence.tsv")
            combined_seq.to_csv(out_comb, sep="\t", index=False)
            console.log(
                f"[green]✓[/] Combined per-sequence table written to {out_comb}"
            )

    console.print("[bold green]Done.[/]")


if __name__ == "__main__":
    main()
