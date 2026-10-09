# ============================================================
# compare_clusters.py
# ============================================================
# Cluster comparison and evaluation script.
#
# Computes external metrics (ARI, Rand Index, precision, recall, F1,
# taxonomic purity, per-set metrics, adjusted/macro precision) and
# pairwise agreement metrics (CCS, CCR, RSS) for one or more
# clustering results.
#
# Dependencies: pandas, numpy, scikit-learn, rich
#
# ============================================================
# METRIC REFERENCE
# ============================================================
#
# --- Family A: Pooled pair-level metrics (classic clustering metrics) ---
#   rand_index      Raw Rand Index = (TP_pairs + TN_pairs) / total_pairs.
#   ARI             Adjusted Rand Index (chance-corrected).
#   precision       TP_pairs / (TP_pairs + FP_pairs).
#   recall          TP_pairs / (TP_pairs + FN_pairs).
#   f1              Harmonic mean of pooled precision and recall.
#
# --- Family B: Adjusted (macro) pair-level metrics ---
#   adjusted_precision      Mean of per-cluster pair precision.
#   adjusted_recall         Mean of per-true-set pair recall.
#   adjusted_f1             Harmonic mean of adjusted_precision & recall.
#   adjusted_rand_index_macro Macro analogue of the Rand Index.
#
# --- Family C: Per-set macro metrics (set recovery) ---
#   per_set_recall          Macro mean of recall_T over all sets.
#   per_set_precision       Macro mean of precision_T over all sets.
#   per_set_f1              Macro mean of f1_T over all sets.
#   median_best_cluster_size Median |best_cluster_T| across all sets.
#   mean_best_cluster_size   Mean   |best_cluster_T| across all sets.
#
# --- Family C-bis: Recovery rate buckets (per-set) ---
#   Counts and fractions of ground-truth sets retrieved at >= 100%, 90%,
#   80%, 70%, 60%, 50%, and <50%.
#
# --- Family D: Sequence-level assignment accuracy ---
#   seq_accuracy_mapped     Fraction of sequences with match=True.
#
# --- Family E: Singleton statistics ---
#   gt_singleton_count / rate   True singletons in dataset.
#   pred_singleton_count / rate Predicted singletons (computed pre-merge).
#
# --- Family F: Cluster-size diagnostics ---
#   Per predicted cluster purity, true set overlap, and pair precision.
#
# --- Family G: Input counts ---
#   Sequence and cluster counts across input, ground truth, and merged sets.
#
# --- Family H: Pairwise agreement between cluster files ---
#   CCS     Cluster Correspondence Score (includes TN pairs).
#   CCR     Cluster Count Ratio (n_clusters_2 / n_clusters_1).
#   RSS     Representative Sequence Stability (vectorized O(N log N)).
#
# --- Family I: Taxonomic purity ---
#   taxonomic_agreement     Majority category consistency across clusters.
# ============================================================

import argparse
import sys
import os
import re
from itertools import combinations
from collections import defaultdict, Counter

import pandas as pd
import numpy as np
from sklearn.metrics import adjusted_rand_score
from rich.console import Console
from rich.table import Table
from rich.progress import Progress, BarColumn, TextColumn, TimeElapsedColumn
from rich import box

console = Console()


# ============================================================
#  Parsing helpers
# ============================================================


def parse_cluster_file(filepath, sep="auto", verbose=False):
    if verbose:
        console.log(f"[cyan]Parsing cluster file:[/] {filepath}")

    if sep == "auto":
        first_line = ""
        with open(filepath, "r") as fh:
            for line in fh:
                if line.strip():
                    first_line = line
                    break
        if not first_line:
            raise ValueError(f"File {filepath} appears to be empty.")
        if "\t" in first_line:
            sep_char = "\t"
            if verbose:
                console.log("[dim]  Detected TAB separator[/]")
        elif "," in first_line:
            sep_char = ","
            if verbose:
                console.log("[dim]  Detected COMMA separator[/]")
        else:
            raise ValueError(
                f"Cannot detect separator in {filepath}. "
                "Use --sep tab|comma to specify explicitly."
            )
    elif sep == "tab":
        sep_char = "\t"
    elif sep == "comma":
        sep_char = ","
    else:
        sep_char = sep

    try:
        df = pd.read_csv(
            filepath,
            header=None,
            sep=sep_char,
            quotechar='"',
            skipinitialspace=True,
            engine="python",
            on_bad_lines="warn",
        )
    except pd.errors.ParserError as e:
        raise ValueError(f"Failed to parse {filepath} as '{sep_char}'-separated: {e}")

    if df.shape[1] < 2:
        raise ValueError(
            f"File {filepath} has {df.shape[1]} column(s) when parsed with "
            f"separator '{sep_char}'. Expected at least 2."
        )

    rep_keywords = {"rep", "representative", "cluster_id", "cluster", "id"}
    set_keywords = {"set", "members", "cluster_set", "sequences", "seqs"}
    first_vals = df.iloc[0].astype(str).str.lower().str.strip()
    is_header = False
    if len(first_vals) >= 2:
        col0 = first_vals.iloc[0]
        col1 = first_vals.iloc[1]
        if (col0 in rep_keywords or any(kw in col0 for kw in rep_keywords)) and (
            col1 in set_keywords or any(kw in col1 for kw in set_keywords)
        ):
            is_header = True
            if verbose:
                console.log("[cyan]Detected header row, skipping.[/]")
    if is_header:
        df = df.iloc[1:].reset_index(drop=True)

    df.columns = ["cluster_id", "cluster_set"] + [
        f"extra_{i}" for i in range(2, df.shape[1])
    ]

    records = []
    skipped = 0
    for row in df.itertuples(index=False):
        cluster = str(row.cluster_id).strip()
        if not cluster or cluster.lower() == "nan":
            skipped += 1
            continue
        members = str(row.cluster_set).split(",")
        for seq in members:
            seq = seq.strip()
            if seq:
                records.append({"sequence_id": seq, "cluster_id": cluster})

    if not records:
        raise ValueError(
            f"No valid sequences found in {filepath}. "
            "Check that the file has two columns and is tab- or comma-separated."
        )
    if verbose:
        console.log(
            f"[green]OK[/] Parsed {len(records)} sequences "
            f"({df.shape[0]} rows, {skipped} skipped) from {filepath}"
        )
    return pd.DataFrame(records)


def parse_ground_truth(filepath, verbose=False):
    if verbose:
        console.log(f"[cyan]Parsing ground truth:[/] {filepath}")
    df = pd.read_csv(filepath, sep="\t")
    for col in ["sequence_id", "true_cluster"]:
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in ground truth file.")
    if verbose:
        console.log(f"[green]OK[/] Loaded {len(df)} ground-truth sequences")
    return df


# ============================================================
#  Pair-counting helpers
# ============================================================


def pair_counts(pred, true):
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
#  Individual metrics
# ============================================================


def rand_index(pred_labels, true_labels):
    """Raw Rand Index (TP + TN) / total_pairs."""
    TP, FP, FN, TN = pair_counts(pred_labels, true_labels)
    total = TP + FP + FN + TN
    return (TP + TN) / total if total > 0 else np.nan


def adjusted_rand_index(true_labels, pred_labels):
    """sklearn Adjusted Rand Index (chance-corrected)."""
    try:
        return adjusted_rand_score(true_labels, pred_labels)
    except Exception:
        return np.nan


def pair_precision_recall_f1(pred_labels, true_labels):
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
    return len(set(labels))


# ============================================================
#  Per-predicted-cluster summary
# ============================================================


def per_cluster_summary(merged_df):
    if merged_df.empty:
        return pd.DataFrame()
    rows = []
    for cid, group in merged_df.groupby("cluster_id"):
        n = len(group)
        if n == 0:
            continue
        counts = group["true_cluster"].value_counts()
        n_true = int(len(counts))
        dominant = counts.idxmax()
        dom_count = int(counts.iloc[0])
        purity = dom_count / n if n > 0 else 0.0
        n_pairs = n * (n - 1) // 2
        tp_pairs = int(sum(int(c) * (int(c) - 1) // 2 for c in counts.values))
        fp_pairs = n_pairs - tp_pairs
        pair_prec = (tp_pairs / n_pairs) if n_pairs > 0 else np.nan
        rows.append(
            {
                "cluster_id": cid,
                "n_members": n,
                "n_true_clusters": n_true,
                "dominant_true_cluster": dominant,
                "dominant_count": dom_count,
                "purity": round(purity, 6),
                "n_pairs": n_pairs,
                "n_tp_pairs": tp_pairs,
                "n_fp_pairs": fp_pairs,
                "pair_precision": (
                    round(pair_prec, 6) if not pd.isna(pair_prec) else np.nan
                ),
            }
        )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df = df.sort_values("n_members", ascending=False).reset_index(drop=True)
    return df


# ============================================================
#  Adjusted (macro) precision / recall
# ============================================================


def adjusted_precision_from_clusters(per_cluster_df):
    if per_cluster_df.empty:
        return np.nan
    non_singleton = per_cluster_df[per_cluster_df["n_members"] >= 2]
    if non_singleton.empty:
        return np.nan
    return float(non_singleton["pair_precision"].mean())


def adjusted_recall_from_merged(merged_df):
    if merged_df.empty:
        return np.nan
    recalls = []
    for _, group in merged_df.groupby("true_cluster"):
        n = len(group)
        if n < 2:
            continue
        total_pairs = n * (n - 1) // 2
        counts = group.groupby("cluster_id").size()
        tp_pairs = int(sum(int(c) * (int(c) - 1) // 2 for c in counts.values))
        if total_pairs > 0:
            recalls.append(tp_pairs / total_pairs)
    return float(np.mean(recalls)) if recalls else np.nan


def harmonic_mean(a, b):
    if a is None or b is None or pd.isna(a) or pd.isna(b):
        return np.nan
    if a + b == 0:
        return 0.0
    return 2 * a * b / (a + b)


# ============================================================
#  Per-set metrics (macro) with recovery buckets
# ============================================================


def per_set_metrics(merged_df):
    empty_result = {
        "per_set_recall": np.nan,
        "per_set_precision": np.nan,
        "per_set_f1": np.nan,
        "median_best_cluster_size": np.nan,
        "mean_best_cluster_size": np.nan,
        "n_sets_total": 0,
        "n_sets_fully_retrieved": 0,
        "n_sets_at_least_90": 0,
        "n_sets_at_least_80": 0,
        "n_sets_at_least_70": 0,
        "n_sets_at_least_60": 0,
        "n_sets_at_least_50": 0,
        "n_sets_below_50": 0,
        "pct_sets_fully_retrieved": np.nan,
        "pct_sets_at_least_90": np.nan,
        "pct_sets_at_least_80": np.nan,
        "pct_sets_at_least_70": np.nan,
        "pct_sets_at_least_60": np.nan,
        "pct_sets_at_least_50": np.nan,
        "pct_sets_below_50": np.nan,
    }
    if merged_df.empty:
        return empty_result

    recalls = []
    precisions = []
    f1s = []
    best_sizes = []
    cluster_sizes = merged_df.groupby("cluster_id").size().to_dict()

    for true_cluster, group in merged_df.groupby("true_cluster"):
        overlap = group.groupby("cluster_id").size()
        if overlap.empty:
            continue
        best_pred = overlap.idxmax()
        tp = int(overlap.loc[best_pred])
        n_true = int(len(group))
        n_pred = int(cluster_sizes.get(best_pred, tp))
        r = tp / n_true if n_true > 0 else 0.0
        p = tp / n_pred if n_pred > 0 else 0.0
        f = (2 * p * r / (p + r)) if (p + r) > 0 else 0.0
        recalls.append(r)
        precisions.append(p)
        f1s.append(f)
        best_sizes.append(n_pred)

    if not recalls:
        return empty_result

    recalls_arr = np.array(recalls)
    n_sets = int(len(recalls_arr))
    n_fully = int((recalls_arr >= 1.0 - 1e-12).sum())
    n_ge_90 = int((recalls_arr >= 0.90).sum())
    n_ge_80 = int((recalls_arr >= 0.80).sum())
    n_ge_70 = int((recalls_arr >= 0.70).sum())
    n_ge_60 = int((recalls_arr >= 0.60).sum())
    n_ge_50 = int((recalls_arr >= 0.50).sum())
    n_lt_50 = int((recalls_arr < 0.50).sum())

    def _pct(n):
        return float(n / n_sets) if n_sets > 0 else np.nan

    return {
        "per_set_recall": float(np.mean(recalls_arr)),
        "per_set_precision": float(np.mean(precisions)),
        "per_set_f1": float(np.mean(f1s)),
        "median_best_cluster_size": float(np.median(best_sizes)),
        "mean_best_cluster_size": float(np.mean(best_sizes)),
        "n_sets_total": n_sets,
        "n_sets_fully_retrieved": n_fully,
        "n_sets_at_least_90": n_ge_90,
        "n_sets_at_least_80": n_ge_80,
        "n_sets_at_least_70": n_ge_70,
        "n_sets_at_least_60": n_ge_60,
        "n_sets_at_least_50": n_ge_50,
        "n_sets_below_50": n_lt_50,
        "pct_sets_fully_retrieved": _pct(n_fully),
        "pct_sets_at_least_90": _pct(n_ge_90),
        "pct_sets_at_least_80": _pct(n_ge_80),
        "pct_sets_at_least_70": _pct(n_ge_70),
        "pct_sets_at_least_60": _pct(n_ge_60),
        "pct_sets_at_least_50": _pct(n_ge_50),
        "pct_sets_below_50": _pct(n_lt_50),
    }


# ============================================================
#  Cluster label matching
# ============================================================


def map_clusters_to_truth(merged_df):
    """
    Maps predicted clusters to ground truth. Strips '_copy_\\d+' suffixes first.
    If the base ID exists in the ground truth, it maps directly.
    Otherwise, falls back to majority voting.
    """
    mapping = {}
    if merged_df.empty:
        return mapping

    known_true = set(merged_df["true_cluster"].unique())

    for cid, group in merged_df.groupby("cluster_id"):
        base_cid = re.sub(r"_copy_\d+$", "", str(cid))

        if base_cid in known_true:
            mapping[cid] = base_cid
        else:
            counts = group["true_cluster"].value_counts()
            if len(counts) > 0:
                mapping[cid] = counts.idxmax()
            else:
                mapping[cid] = cid
    return mapping


# ============================================================
#  Pairwise agreement
# ============================================================


def ccs(labels1, labels2):
    """
    Cluster Correspondence Score.
    Fraction of sequence pairs that both runs agree on (TP + TN) / total_pairs.
    """
    n = len(labels1)
    if n < 2:
        return np.nan
    TP, FP, FN, TN = pair_counts(labels1, labels2)
    total_pairs = TP + FP + FN + TN
    return (TP + TN) / total_pairs if total_pairs > 0 else 0.0


def ccr(labels1, labels2):
    """Cluster Count Ratio: n_clusters_2 / n_clusters_1."""
    n1 = len(set(labels1))
    n2 = len(set(labels2))
    return n2 / n1 if n1 > 0 else np.nan


def rss(cluster_df1, cluster_df2):
    """
    Representative Sequence Stability.
    Vectorized implementation O(N log N) comparing representative sequence preservation.
    """
    rep1 = cluster_df1.groupby("cluster_id")["sequence_id"].first().to_dict()
    rep2 = cluster_df2.groupby("cluster_id")["sequence_id"].first().to_dict()

    merged = pd.merge(cluster_df1, cluster_df2, on="sequence_id", suffixes=("_1", "_2"))
    if merged.empty:
        return 0.0

    overlap_counts = (
        merged.groupby(["cluster_id_1", "cluster_id_2"])
        .size()
        .reset_index(name="overlap")
    )
    sizes1 = cluster_df1.groupby("cluster_id").size().to_dict()

    best_matches = overlap_counts.sort_values(
        "overlap", ascending=False
    ).drop_duplicates("cluster_id_1")

    stable = 0
    total = len(sizes1)

    for row in best_matches.itertuples():
        c1, c2, ov = row.cluster_id_1, row.cluster_id_2, row.overlap
        size1 = sizes1.get(c1, 0)

        if (ov / size1) >= 0.5 and rep1.get(c1) == rep2.get(c2):
            stable += 1

    return stable / total if total > 0 else 0.0


# ============================================================
#  Aggregated metrics
# ============================================================


def cluster_external_metrics(
    cluster_df, ground_truth_df, sample_id=None, file_path=None
):
    n_sequences_input = int(cluster_df["sequence_id"].nunique())
    n_clusters_total = int(cluster_df["cluster_id"].nunique())
    n_sequences_gt = int(ground_truth_df["sequence_id"].nunique())

    # Pre-calculate sizes and singletons directly on input datasets before inner merge
    pred_sizes_input = cluster_df.groupby("cluster_id").size().to_dict()
    pred_singletons_input = {cid for cid, sz in pred_sizes_input.items() if sz == 1}

    gt_sizes_input = ground_truth_df.groupby("true_cluster").size().to_dict()
    gt_singletons_input = {tid for tid, sz in gt_sizes_input.items() if sz == 1}

    merged = pd.merge(cluster_df, ground_truth_df, on="sequence_id", how="inner")
    n_sequences_merged = int(merged["sequence_id"].nunique()) if not merged.empty else 0

    if merged.empty:
        return (
            {
                "n_sequences_input": n_sequences_input,
                "n_sequences_gt": n_sequences_gt,
                "n_sequences_merged": n_sequences_merged,
                "n_clusters_input": n_clusters_total,
            },
            pd.DataFrame(),
            pd.DataFrame(),
        )

    pred_labels = merged["cluster_id"].values
    true_labels = merged["true_cluster"].values

    ri = rand_index(pred_labels, true_labels)
    ari = adjusted_rand_index(true_labels, pred_labels)
    prec, rec, f1 = pair_precision_recall_f1(pred_labels, true_labels)
    n_clusters = cluster_count(pred_labels)
    tax_agr = taxonomic_agreement(merged)

    per_cluster_df = per_cluster_summary(merged)

    adj_prec = adjusted_precision_from_clusters(per_cluster_df)
    adj_rec = adjusted_recall_from_merged(merged)
    adj_f1 = harmonic_mean(adj_prec, adj_rec)
    adj_ri_macro = harmonic_mean(adj_prec, adj_rec)

    pset = per_set_metrics(merged)

    mapping = map_clusters_to_truth(merged)

    gt_set_sizes = None
    if "set_id" in ground_truth_df.columns:
        gt_set_sizes = ground_truth_df.groupby("set_id").size()

    cluster_purity_map = {}
    cluster_n_true_map = {}
    for cid, group in merged.groupby("cluster_id"):
        counts = group["true_cluster"].value_counts()
        cluster_purity_map[cid] = float(counts.iloc[0] / counts.sum())
        cluster_n_true_map[cid] = int(len(counts))

    merged["cluster_assigned"] = merged["cluster_id"]
    merged["cluster_assigned_mapped"] = merged["cluster_id"].map(mapping)
    merged["match"] = merged["cluster_assigned_mapped"] == merged["true_cluster"]

    # Accurate singleton assignment using pre-merge inputs
    merged["gt_is_singleton"] = merged["true_cluster"].isin(gt_singletons_input)
    merged["pred_is_singleton"] = merged["cluster_id"].isin(pred_singletons_input)
    merged["pred_cluster_size"] = (
        merged["cluster_id"].map(pred_sizes_input).fillna(0).astype(int)
    )
    merged["pred_cluster_purity"] = (
        merged["cluster_id"].map(cluster_purity_map).fillna(0.0)
    )
    merged["pred_cluster_n_true"] = (
        merged["cluster_id"].map(cluster_n_true_map).fillna(0).astype(int)
    )

    if "set_id" in merged.columns:
        if gt_set_sizes is not None:
            merged["set_size_gt"] = (
                merged["set_id"].map(gt_set_sizes).fillna(0).astype(int)
            )
        else:
            merged["set_size_gt"] = merged.groupby("set_id")["sequence_id"].transform(
                "size"
            )
        correct_by_set = merged.loc[merged["match"]].groupby("set_id").size()
        merged["set_size_correctly_mapped"] = (
            merged["set_id"].map(correct_by_set).fillna(0).astype(int)
        )

    seq_accuracy = float(merged["match"].mean()) if len(merged) > 0 else np.nan

    gt_singleton_rate = float(merged["gt_is_singleton"].mean())
    pred_singleton_rate = float(merged["pred_is_singleton"].mean())
    gt_singleton_count = int(merged["gt_is_singleton"].sum())
    pred_singleton_count = int(merged["pred_is_singleton"].sum())

    metrics = {
        "n_sequences_input": n_sequences_input,
        "n_sequences_gt": n_sequences_gt,
        "n_sequences_merged": n_sequences_merged,
        "n_clusters_input": n_clusters_total,
        "seq_accuracy_mapped": seq_accuracy,
        "rand_index": ri,
        "ARI": ari,
        "adjusted_rand_index_macro": adj_ri_macro,
        "adjusted_precision": adj_prec,
        "adjusted_recall": adj_rec,
        "adjusted_f1": adj_f1,
        "per_set_recall": pset["per_set_recall"],
        "per_set_precision": pset["per_set_precision"],
        "per_set_f1": pset["per_set_f1"],
        "median_best_cluster_size": pset["median_best_cluster_size"],
        "mean_best_cluster_size": pset["mean_best_cluster_size"],
        "n_sets_total": pset["n_sets_total"],
        "n_sets_fully_retrieved": pset["n_sets_fully_retrieved"],
        "n_sets_at_least_90": pset["n_sets_at_least_90"],
        "n_sets_at_least_80": pset["n_sets_at_least_80"],
        "n_sets_at_least_70": pset["n_sets_at_least_70"],
        "n_sets_at_least_60": pset["n_sets_at_least_60"],
        "n_sets_at_least_50": pset["n_sets_at_least_50"],
        "n_sets_below_50": pset["n_sets_below_50"],
        "pct_sets_fully_retrieved": pset["pct_sets_fully_retrieved"],
        "pct_sets_at_least_90": pset["pct_sets_at_least_90"],
        "pct_sets_at_least_80": pset["pct_sets_at_least_80"],
        "pct_sets_at_least_70": pset["pct_sets_at_least_70"],
        "pct_sets_at_least_60": pset["pct_sets_at_least_60"],
        "pct_sets_at_least_50": pset["pct_sets_at_least_50"],
        "pct_sets_below_50": pset["pct_sets_below_50"],
        "n_clusters": n_clusters,
        "precision": prec,
        "recall": rec,
        "f1": f1,
        "taxonomic_agreement": tax_agr,
        "gt_singleton_count": gt_singleton_count,
        "gt_singleton_rate": gt_singleton_rate,
        "pred_singleton_count": pred_singleton_count,
        "pred_singleton_rate": pred_singleton_rate,
    }

    if sample_id:
        merged["sample_id"] = sample_id
    if file_path:
        metrics["file_path"] = file_path

    preferred_order = [
        "sequence_id",
        "true_cluster",
        "cluster_assigned",
        "cluster_assigned_mapped",
        "match",
        "gt_is_singleton",
        "pred_is_singleton",
        "set_id",
        "set_size_gt",
        "set_size_correctly_mapped",
        "pred_cluster_size",
        "pred_cluster_purity",
        "pred_cluster_n_true",
        "source_genome",
        "source_fragment",
        "length",
        "mutation_rate",
        "mutation_summary",
        "mutation_positions_short",
        "mutation_types_short",
        "copy_number",
        "is_duplicate",
        "source_category",
        "is_synthetic",
        "sample_id",
    ]
    ordered = [c for c in preferred_order if c in merged.columns]
    remaining = [c for c in merged.columns if c not in ordered]
    merged = merged[ordered + remaining]

    return metrics, merged, per_cluster_df


def cluster_agreement_metrics(cluster_df1, cluster_df2):
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
#  CLI helpers
# ============================================================


def show_parameters(args):
    table = Table(title="Configuration", box=box.ROUNDED)
    table.add_column("Parameter", style="cyan", no_wrap=True)
    table.add_column("Value", style="white")
    table.add_row("Ground truth", args.ground_truth)
    table.add_row("Number of cluster files", str(len(args.cluster_files)))
    table.add_row(
        "Sample names",
        ", ".join(args.sample_names) if args.sample_names else "auto-detected",
    )
    table.add_row("Separator", args.sep)
    table.add_row("Output directory", args.output_dir)
    table.add_row("Verbose", str(args.verbose))
    table.add_row("Dry-run", str(args.dry_run))
    console.print(table)


# ============================================================
#  Main
# ============================================================


def main():
    parser = argparse.ArgumentParser(
        description="Compare clustering results against ground truth and evaluate pairwise agreement."
    )
    parser.add_argument(
        "-g", "--ground-truth", required=True, help="Path to ground truth TSV file."
    )
    parser.add_argument(
        "-c",
        "--cluster-files",
        nargs="+",
        required=True,
        help="One or more cluster files.",
    )
    parser.add_argument(
        "-s",
        "--sample-names",
        nargs="+",
        default=None,
        help="Optional sample labels for each cluster file.",
    )
    parser.add_argument(
        "--sep",
        choices=["auto", "tab", "comma"],
        default="auto",
        help="Separator format (default: auto).",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default="cluster_comparison_results",
        help="Directory to save output CSVs.",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Print verbose execution log."
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Parse files without writing output CSVs.",
    )

    args = parser.parse_args()
    show_parameters(args)

    if args.sample_names and len(args.sample_names) != len(args.cluster_files):
        console.print(
            "[red]Error: --sample-names length must match --cluster-files length.[/]"
        )
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)
    gt_df = parse_ground_truth(args.ground_truth, verbose=args.verbose)

    metrics_list = []
    cluster_dfs = []

    with Progress(
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TimeElapsedColumn(),
        console=console,
    ) as progress:
        task = progress.add_task(
            "Evaluating cluster files...", total=len(args.cluster_files)
        )

        for i, cfile in enumerate(args.cluster_files):
            sample_id = (
                args.sample_names[i]
                if args.sample_names
                else os.path.basename(cfile).split(".")[0]
            )
            c_df = parse_cluster_file(cfile, sep=args.sep, verbose=args.verbose)
            cluster_dfs.append((sample_id, c_df))

            metrics, merged_df, per_cluster_df = cluster_external_metrics(
                c_df, gt_df, sample_id=sample_id, file_path=cfile
            )
            metrics_list.append(metrics)

            if not args.dry_run and not merged_df.empty:
                merged_df.to_csv(
                    os.path.join(args.output_dir, f"{sample_id}_sequence_matches.csv"),
                    index=False,
                )
                per_cluster_df.to_csv(
                    os.path.join(
                        args.output_dir, f"{sample_id}_per_cluster_summary.csv"
                    ),
                    index=False,
                )
            progress.update(task, advance=1)

    # Save aggregated external metrics
    metrics_summary_df = pd.DataFrame(metrics_list)
    if not args.dry_run:
        metrics_summary_df.to_csv(
            os.path.join(args.output_dir, "external_metrics_summary.csv"), index=False
        )

    # Display summary table in CLI
    summary_table = Table(title="Clustering Metric Summary", box=box.ROUNDED)
    summary_table.add_column("Sample", style="cyan")
    summary_table.add_column("Mapped Acc", style="magenta")
    summary_table.add_column("ARI", style="green")
    summary_table.add_column("Adj Precision", style="yellow")
    summary_table.add_column("Adj Recall", style="yellow")
    summary_table.add_column("Adj F1", style="bold yellow")
    summary_table.add_column("Per-Set F1", style="blue")

    for m in metrics_list:
        summary_table.add_row(
            str(m.get("file_path", "Sample")),
            f"{m.get('seq_accuracy_mapped', 0):.4f}",
            f"{m.get('ARI', 0):.4f}",
            f"{m.get('adjusted_precision', 0):.4f}",
            f"{m.get('adjusted_recall', 0):.4f}",
            f"{m.get('adjusted_f1', 0):.4f}",
            f"{m.get('per_set_f1', 0):.4f}",
        )
    console.print(summary_table)

    # Pairwise agreement metrics if >= 2 runs
    if len(cluster_dfs) >= 2:
        console.print("\n[cyan]Computing pairwise agreement metrics...[/]")
        pairwise_rows = []
        for (id1, cdf1), (id2, cdf2) in combinations(cluster_dfs, 2):
            agree = cluster_agreement_metrics(cdf1, cdf2)
            pairwise_rows.append({"Run_1": id1, "Run_2": id2, **agree})

        pairwise_df = pd.DataFrame(pairwise_rows)
        if not args.dry_run:
            pairwise_df.to_csv(
                os.path.join(args.output_dir, "pairwise_agreement.csv"), index=False
            )

        pair_table = Table(title="Pairwise Run Agreement", box=box.ROUNDED)
        pair_table.add_column("Run 1", style="cyan")
        pair_table.add_column("Run 2", style="cyan")
        pair_table.add_column("CCS", style="green")
        pair_table.add_column("CCR", style="yellow")
        pair_table.add_column("RSS", style="magenta")

        for row in pairwise_rows:
            pair_table.add_row(
                row["Run_1"],
                row["Run_2"],
                f"{row['CCS']:.4f}",
                f"{row['CCR']:.4f}",
                f"{row['RSS']:.4f}",
            )
        console.print(pair_table)

    console.print(f"[green]Done![/] Output saved to: [bold]{args.output_dir}[/]")


if __name__ == "__main__":
    main()
