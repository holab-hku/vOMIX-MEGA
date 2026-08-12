#!/usr/bin/env python3
"""
Cluster comparison and evaluation script.

This script computes external metrics (ARI, precision, recall, F1, taxonomic purity)
and pairwise agreement metrics (CCS, CCR, RSS) for one or more clustering results,
comparing them against a ground‑truth clustering.

Usage:
    python compare_clusters.py \
        --ground-truth ground_truth.tsv \
        --cluster-files file1.csv file2.csv ... \
        --sample-names exp1 exp2 ... \
        --output-dir results/ \
        [--verbose] [--dry-run]

Dependencies:
    - pandas, numpy, scikit-learn, rich
    (see environment.yml or install via pip)
"""

import argparse
import sys
import os
from itertools import combinations
from collections import defaultdict, Counter
from pathlib import Path

import pandas as pd
import numpy as np
from sklearn.metrics import adjusted_rand_score
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.progress import Progress, BarColumn, TextColumn, TimeElapsedColumn
from rich import box

# ============================================================
#  Global console for rich output
# ============================================================

console = Console()

# ============================================================
#  Helper functions
# ============================================================


def parse_cluster_file(filepath, verbose=False):
    """
    Read a clustering CSV robustly.
    The first column is always the cluster representative ID.
    The second column is always the comma‑separated list of member sequence IDs.
    A header row is automatically detected and skipped if present.
    """
    if verbose:
        console.log(f"[cyan]Parsing cluster file:[/] {filepath}")

    # Read with header=None so we can inspect the first row
    df = pd.read_csv(filepath, header=None, quotechar='"', skipinitialspace=True)

    # Keywords to detect a header row
    rep_keywords = ["rep", "representative", "cluster", "id"]
    set_keywords = ["set", "members", "cluster_set", "sequences", "members"]

    # Check first row strings (lowercase)
    first_vals = df.iloc[0].astype(str).str.lower()
    is_header = False
    if len(first_vals) >= 2:
        col0 = first_vals[0]
        col1 = first_vals[1]
        if any(kw in col0 for kw in rep_keywords) and any(
            kw in col1 for kw in set_keywords
        ):
            is_header = True
            if verbose:
                console.log("[cyan]Detected header row, skipping.[/]")

    if is_header:
        df = df.iloc[1:].reset_index(drop=True)

    # Ensure we have at least two columns
    if df.shape[1] < 2:
        raise ValueError(f"File {filepath} has fewer than 2 columns.")

    # Assign column names by position
    df.columns = ["cluster_id", "cluster_set"] + [
        f"extra_{i}" for i in range(2, df.shape[1])
    ]

    # Expand members
    records = []
    for _, row in df.iterrows():
        cluster = str(row["cluster_id"]).strip()
        # In case cluster_id is empty, skip
        if not cluster:
            continue
        members = str(row["cluster_set"]).split(",")
        for seq in members:
            seq = seq.strip()
            if seq:
                records.append({"sequence_id": seq, "cluster_id": cluster})

    if verbose:
        console.log(f"[green]✓[/] Parsed {len(records)} sequences from {filepath}")

    return pd.DataFrame(records)


def parse_ground_truth(filepath, verbose=False):
    """Read ground truth TSV. Must have 'sequence_id' and 'true_cluster'."""
    if verbose:
        console.log(f"[cyan]Parsing ground truth:[/] {filepath}")
    df = pd.read_csv(filepath, sep="\t")
    for col in ["sequence_id", "true_cluster"]:
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in ground truth file.")
    if verbose:
        console.log(f"[green]✓[/] Loaded {len(df)} ground‑truth sequences")
    return df


def pair_counts(pred, true):
    """Return TP, FP, FN, TN pair counts between two label arrays."""
    n = len(pred)
    if n == 0:
        return 0, 0, 0, 0
    pred_counts = Counter(pred)
    true_counts = Counter(true)
    joint = defaultdict(Counter)
    for p, t in zip(pred, true):
        joint[p][t] += 1

    TP = 0
    for p, true_counter in joint.items():
        for t, cnt in true_counter.items():
            if cnt >= 2:
                TP += cnt * (cnt - 1) // 2

    FP = 0
    for p, cnt in pred_counts.items():
        if cnt >= 2:
            FP += cnt * (cnt - 1) // 2
    FP -= TP

    FN = 0
    for t, cnt in true_counts.items():
        if cnt >= 2:
            FN += cnt * (cnt - 1) // 2
    FN -= TP

    total_pairs = n * (n - 1) // 2
    TN = total_pairs - TP - FP - FN
    return TP, FP, FN, TN


# ============================================================
#  Individual metric functions (external)
# ============================================================


def adjusted_rand_index(true_labels, pred_labels):
    """Compute Adjusted Rand Index."""
    try:
        return adjusted_rand_score(true_labels, pred_labels)
    except Exception:
        return np.nan


def pair_precision_recall_f1(pred_labels, true_labels):
    """Compute precision, recall, and F1 based on pair counting."""
    TP, FP, FN, _ = pair_counts(pred_labels, true_labels)
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0.0
    f1 = (
        2 * precision * recall / (precision + recall)
        if (precision + recall) > 0
        else 0.0
    )
    return precision, recall, f1


def taxonomic_agreement(merged_df):
    """
    Compute average purity of each cluster with respect to 'source_category'.
    Returns float or np.nan if column missing.
    """
    if "source_category" not in merged_df.columns:
        return np.nan
    purity_scores = []
    for _, group in merged_df.groupby("cluster_id")["source_category"]:
        if len(group) > 0:
            mode = group.mode()
            if len(mode) > 0:
                majority = mode.iloc[0]
                purity = (group == majority).sum() / len(group)
                purity_scores.append(purity)
    return np.mean(purity_scores) if purity_scores else np.nan


def cluster_count(labels):
    """Return number of unique clusters."""
    return len(set(labels))


# ============================================================
#  Individual metric functions (pairwise agreement)
# ============================================================


def ccs(labels1, labels2):
    """Cluster Correspondence Score (CCS)."""
    n = len(labels1)
    if n < 2:
        return np.nan
    total_pairs = n * (n - 1) // 2
    TP, _, _, _ = pair_counts(labels1, labels2)
    return TP / total_pairs if total_pairs > 0 else 0.0


def ccr(labels1, labels2):
    """Cluster Count Ratio (CCR)."""
    n1 = len(set(labels1))
    n2 = len(set(labels2))
    return n2 / n1 if n1 > 0 else np.nan


def rss(cluster_df1, cluster_df2):
    """
    Representative Sequence Stability (RSS).
    Uses first alphabetically as representative.
    """

    def cluster_to_seq(df):
        d = defaultdict(list)
        for _, row in df.iterrows():
            d[row["cluster_id"]].append(row["sequence_id"])
        return d

    map1 = cluster_to_seq(cluster_df1)
    map2 = cluster_to_seq(cluster_df2)
    rep1 = {c: sorted(seqs)[0] for c, seqs in map1.items()}
    rep2 = {c: sorted(seqs)[0] for c, seqs in map2.items()}

    stable = 0
    total = 0
    for c1, seqs1 in map1.items():
        set1 = set(seqs1)
        best_c2 = None
        best_overlap = 0
        for c2, seqs2 in map2.items():
            overlap = len(set1 & set(seqs2))
            if overlap > best_overlap:
                best_overlap = overlap
                best_c2 = c2
        if best_c2 is not None:
            total += 1
            overlap_frac = best_overlap / len(seqs1) if len(seqs1) > 0 else 0.0
            if overlap_frac >= 0.5 and rep1[c1] == rep2[best_c2]:
                stable += 1
    return stable / total if total > 0 else 0.0


# ============================================================
#  Aggregated metric functions
# ============================================================


def cluster_external_metrics(
    cluster_df, ground_truth_df, sample_id=None, file_path=None
):
    """Compute all external metrics; returns (metrics_dict, merged_df)."""
    merged = pd.merge(cluster_df, ground_truth_df, on="sequence_id", how="inner")
    if merged.empty:
        return {}, pd.DataFrame()

    pred_labels = merged["cluster_id"].values
    true_labels = merged["true_cluster"].values

    metrics = {
        "ARI": adjusted_rand_index(true_labels, pred_labels),
        "n_clusters": cluster_count(pred_labels),
        "precision": pair_precision_recall_f1(pred_labels, true_labels)[0],
        "recall": pair_precision_recall_f1(pred_labels, true_labels)[1],
        "f1": pair_precision_recall_f1(pred_labels, true_labels)[2],
        "taxonomic_agreement": taxonomic_agreement(merged),
    }

    merged["match"] = merged["cluster_id"] == merged["true_cluster"]
    if sample_id:
        merged["sample_id"] = sample_id
    if file_path:
        merged["file_path"] = file_path

    return metrics, merged


def cluster_agreement_metrics(cluster_df1, cluster_df2):
    """Compute pairwise agreement metrics; returns dict with CCS, CCR, RSS."""
    merged = pd.merge(cluster_df1, cluster_df2, on="sequence_id", suffixes=("_1", "_2"))
    if merged.empty or len(merged) < 2:
        return {"CCS": np.nan, "CCR": np.nan, "RSS": np.nan}

    labels1 = merged["cluster_id_1"].values
    labels2 = merged["cluster_id_2"].values

    return {
        "CCS": ccs(labels1, labels2),
        "CCR": ccr(labels1, labels2),
        "RSS": rss(cluster_df1, cluster_df2),
    }


# ============================================================
#  Main CLI with rich logging and progress
# ============================================================


def show_parameters(args):
    """Display a summary table of input parameters."""
    table = Table(title="Configuration", box=box.ROUNDED)
    table.add_column("Parameter", style="cyan", no_wrap=True)
    table.add_column("Value", style="white")
    table.add_row("Ground truth", args.ground_truth)
    table.add_row("Number of cluster files", str(len(args.cluster_files)))
    table.add_row(
        "Sample names",
        ", ".join(args.sample_names) if args.sample_names else "auto‑detected",
    )
    table.add_row("Output directory", args.output_dir)
    table.add_row("Verbose", str(args.verbose))
    table.add_row("Dry‑run", str(args.dry_run))
    console.print(table)


def main():
    parser = argparse.ArgumentParser(
        description="Compare clustering outputs against ground truth.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--ground-truth",
        required=True,
        help="Path to ground truth TSV file (must have 'sequence_id' and 'true_cluster').",
    )
    parser.add_argument(
        "--cluster-files",
        nargs="+",
        required=True,
        help="List of clustering CSV files to evaluate.",
    )
    parser.add_argument(
        "--sample-names",
        nargs="*",
        default=None,
        help="Optional list of sample IDs (same order as cluster files). "
        "If not provided, basenames are used.",
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory to write output files (default: current dir).",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print detailed progress and debug information.",
    )
    parser.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="Simulate execution without writing any files.",
    )
    args = parser.parse_args()

    # Early validation
    if not os.path.isfile(args.ground_truth):
        console.print(
            f"[bold red]ERROR:[/] Ground truth file not found: {args.ground_truth}"
        )
        sys.exit(1)

    for f in args.cluster_files:
        if not os.path.isfile(f):
            console.print(f"[bold red]ERROR:[/] Cluster file not found: {f}")
            sys.exit(1)

    # Show configuration
    console.rule("[bold green]Cluster Comparison Tool")
    show_parameters(args)

    if args.dry_run:
        console.print("[bold yellow]DRY RUN – no files will be written.[/]")

    # Prepare output directory
    if not args.dry_run:
        os.makedirs(args.output_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # Parse files with progress bar
    # ------------------------------------------------------------------
    console.print("\n[bold]Loading data...[/]")
    gt = None
    try:
        gt = parse_ground_truth(args.ground_truth, verbose=args.verbose)
    except Exception as e:
        console.print(f"[bold red]Failed to parse ground truth:[/] {e}")
        sys.exit(1)

    cluster_dfs = []
    sample_ids = []
    file_paths = []

    with Progress(
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TimeElapsedColumn(),
        console=console,
    ) as progress:
        task = progress.add_task(
            "[cyan]Parsing cluster files...", total=len(args.cluster_files)
        )
        for i, f in enumerate(args.cluster_files):
            try:
                df = parse_cluster_file(f, verbose=args.verbose)
                cluster_dfs.append(df)
                file_paths.append(os.path.abspath(f))
                sid = (
                    args.sample_names[i]
                    if args.sample_names and i < len(args.sample_names)
                    else os.path.splitext(os.path.basename(f))[0]
                )
                sample_ids.append(sid)
                progress.advance(task)
            except Exception as e:
                console.print(f"[bold red]Error parsing {f}:[/] {e}")
                if args.verbose:
                    console.print_exception()
                # Continue with other files

    if len(cluster_dfs) == 0:
        console.print("[bold red]No valid cluster files loaded. Exiting.[/]")
        sys.exit(1)

    # ------------------------------------------------------------------
    # 1. External metrics per file
    # ------------------------------------------------------------------
    console.print("\n[bold]Computing external metrics...[/]")
    external_summary = []
    merged_tables = []

    with Progress(
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TimeElapsedColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("[cyan]External metrics", total=len(cluster_dfs))
        for sid, df, fpath in zip(sample_ids, cluster_dfs, file_paths):
            metrics, merged = cluster_external_metrics(
                df, gt, sample_id=sid, file_path=fpath
            )
            if merged.empty:
                console.print(
                    f"[yellow]Warning: No overlap for sample {sid}. Skipping metrics.[/]"
                )
            else:
                metrics["sample_id"] = sid
                external_summary.append(metrics)
                merged_tables.append(merged)
                if not args.dry_run:
                    out_merged = os.path.join(
                        args.output_dir, f"{sid}_merged_with_gt.tsv"
                    )
                    merged.to_csv(out_merged, sep="\t", index=False)
            progress.advance(task)

    if external_summary:
        summary_df = pd.DataFrame(external_summary)
        summary_df = summary_df[
            [
                "sample_id",
                "ARI",
                "n_clusters",
                "precision",
                "recall",
                "f1",
                "taxonomic_agreement",
            ]
        ]
        if not args.dry_run:
            summary_path = os.path.join(args.output_dir, "external_metrics_summary.tsv")
            summary_df.to_csv(summary_path, sep="\t", index=False)
            console.log(
                f"[green]✓[/] External metrics summary written to [cyan]{summary_path}[/]"
            )
        else:
            console.log("[yellow]Dry‑run: would write external_metrics_summary.tsv[/]")
    else:
        console.print("[yellow]No external metrics computed.[/]")

    # ------------------------------------------------------------------
    # 2. Pairwise agreement (all pairs)
    # ------------------------------------------------------------------
    if len(cluster_dfs) >= 2:
        console.print("\n[bold]Computing pairwise agreement metrics...[/]")
        pairwise_results = []
        total_pairs = len(list(combinations(range(len(cluster_dfs)), 2)))
        with Progress(
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeElapsedColumn(),
            console=console,
        ) as progress:
            task = progress.add_task("[cyan]Pairwise comparisons", total=total_pairs)
            for (i, df1, sid1), (j, df2, sid2) in combinations(
                zip(range(len(cluster_dfs)), cluster_dfs, sample_ids), 2
            ):
                ag = cluster_agreement_metrics(df1, df2)
                if ag:
                    record = {
                        "sample1": sid1,
                        "sample2": sid2,
                        "CCS": ag["CCS"],
                        "CCR": ag["CCR"],
                        "RSS": ag["RSS"],
                    }
                    pairwise_results.append(record)
                progress.advance(task)

        if pairwise_results:
            pair_df = pd.DataFrame(pairwise_results)
            pair_df = pair_df[["sample1", "sample2", "CCS", "CCR", "RSS"]]
            if not args.dry_run:
                pair_path = os.path.join(args.output_dir, "pairwise_agreement.tsv")
                pair_df.to_csv(pair_path, sep="\t", index=False)
                console.log(
                    f"[green]✓[/] Pairwise agreement written to [cyan]{pair_path}[/]"
                )
            else:
                console.log("[yellow]Dry‑run: would write pairwise_agreement.tsv[/]")
        else:
            console.print("[yellow]No pairwise metrics computed.[/]")
    else:
        console.print("[yellow]Only one cluster file – skipping pairwise metrics.[/]")

    # ------------------------------------------------------------------
    # 3. Combined merged table (all samples)
    # ------------------------------------------------------------------
    if merged_tables:
        combined = pd.concat(merged_tables, ignore_index=True)
        if not args.dry_run:
            combined_path = os.path.join(args.output_dir, "all_merged_with_gt.tsv")
            combined.to_csv(combined_path, sep="\t", index=False)
            console.log(
                f"[green]✓[/] Combined merged table written to [cyan]{combined_path}[/]"
            )
        else:
            console.log("[yellow]Dry‑run: would write all_merged_with_gt.tsv[/]")

    # ------------------------------------------------------------------
    # Final summary
    # ------------------------------------------------------------------
    console.rule("[bold green]Summary")
    success_count = len(external_summary)
    failure_count = len(cluster_dfs) - success_count
    table = Table(title="Results Summary", box=box.SIMPLE)
    table.add_column("Metric", style="cyan")
    table.add_column("Value", style="white")
    table.add_row("Samples processed", str(len(cluster_dfs)))
    table.add_row(
        "Successful external metrics", f"{success_count} / {len(cluster_dfs)}"
    )
    table.add_row(
        "Failures", str(failure_count) if failure_count > 0 else "[green]None[/]"
    )
    table.add_row(
        "Output directory",
        (
            args.output_dir
            if not args.dry_run
            else "[yellow]Dry‑run (no files written)[/]"
        ),
    )
    console.print(table)

    if failure_count > 0:
        console.print("[yellow]⚠ Some samples failed – see messages above.[/]")
    else:
        console.print("[green]✓ All samples processed successfully.[/]")


if __name__ == "__main__":
    main()
