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
# There are several "families" of metrics below. They can disagree
# strongly when the clustering has a long-tailed distribution of
# predicted cluster sizes (a few giant clusters + many small ones).
#
# --- Family A: Pooled pair-level metrics (classic clustering metrics) ---
#
#   rand_index      Raw Rand Index = (TP_pairs + TN_pairs) / total_pairs.
#                   Range [0, 1]. Not chance-corrected.
#   ARI             Adjusted Rand Index (chance-corrected).
#                   sklearn.metrics.adjusted_rand_score.
#                   Range roughly [-1, 1]; 0 = random, 1 = perfect.
#   precision       TP_pairs / (TP_pairs + FP_pairs).
#   recall          TP_pairs / (TP_pairs + FN_pairs).
#   f1              Harmonic mean of pooled precision and recall.
#
#   A "pair" is an unordered pair of sequences (i, j).
#   - TP pair: predicted same cluster AND truly same cluster.
#   - FP pair: predicted same cluster BUT truly different clusters.
#   - FN pair: predicted different clusters BUT truly same cluster.
#   - TN pair: predicted different clusters AND truly different.
#
#   WARNING: Each cluster of size s contributes C(s, 2) = s*(s-1)/2
#   pairs to the denominator. A single cluster of size 86 alone casts
#   3,655 pair votes -- as many as 610 clean size-4 clusters. These
#   metrics are therefore dominated by the largest predicted clusters
#   and are NOT representative of typical behaviour.
#
# --- Family B: Adjusted (macro) pair-level metrics ---
#
#   adjusted_precision      Mean of per-cluster pair precision over
#                           non-singleton clusters. Each cluster gets
#                           ONE vote regardless of size.
#                           Formula per cluster: TP_pairs / n_pairs.
#
#   adjusted_recall         Mean of per-true-set pair recall over
#                           non-singleton true sets. Each true set
#                           gets ONE vote regardless of size.
#                           Formula per true set: TP_pairs / n_pairs.
#
#   adjusted_f1             Harmonic mean of adjusted_precision and
#                           adjusted_recall.
#
#   adjusted_rand_index_macro
#                           Harmonic mean of adjusted_precision and
#                           adjusted_recall; a macro analogue of the
#                           Rand Index that gives each cluster and each
#                           true set an equal vote. Range [0, 1].
#                           (Same value as adjusted_f1; kept as a
#                           separately-named column for clarity.)
#
#   These are the "fair" versions of precision/recall when the
#   algorithm has a handful of over-merged clusters. Use these for
#   any headline claim about overall clustering quality.
#
# --- Family C: Per-set macro metrics (set recovery) ---
#
#   For each ground-truth set T:
#     best_cluster_T      Predicted cluster with the largest overlap with T.
#     recall_T            |T ∩ best_cluster_T| / |T|
#     precision_T         |T ∩ best_cluster_T| / |best_cluster_T|
#     f1_T                Harmonic mean of the above.
#
#   per_set_recall             Macro mean of recall_T over all sets.
#   per_set_precision          Macro mean of precision_T over all sets.
#   per_set_f1                 Macro mean of f1_T over all sets.
#   median_best_cluster_size   Median |best_cluster_T| across all sets.
#   mean_best_cluster_size     Mean   |best_cluster_T| across all sets.
#
#   These answer: "for each source set, how much of it survived into
#   one predicted cluster, and how clean was that cluster?".
#
# --- Family C-bis: Recovery rate buckets (per-set) ---
#
#   n_sets_total            Total number of distinct ground-truth sets.
#   n_sets_fully_retrieved  Sets with recall_T == 1.0 (100% grouped).
#   n_sets_at_least_90      Sets with recall_T >= 0.90.
#   n_sets_at_least_80      Sets with recall_T >= 0.80.
#   n_sets_at_least_70      Sets with recall_T >= 0.70.
#   n_sets_at_least_60      Sets with recall_T >= 0.60.
#   n_sets_at_least_50      Sets with recall_T >= 0.50.
#   n_sets_below_50         Sets with recall_T <  0.50 (failures).
#   pct_sets_*              Same counts expressed as fractions.
#
# --- Family D: Sequence-level assignment accuracy ---
#
#   For every predicted cluster, its members are mapped to the
#   majority true_cluster ("cluster_assigned_mapped"). For each
#   sequence, match=True iff its mapped label equals its own
#   true_cluster.
#
#   seq_accuracy_mapped     Fraction of sequences with match=True.
#
#   Lenient -- a 60/40 mixed cluster still shows 60% True.
#
# --- Family E: Singleton statistics ---
#
#   gt_singleton_count      # sequences whose true set has size 1.
#   gt_singleton_rate       gt_singleton_count / n_sequences_merged.
#   pred_singleton_count    # sequences in a size-1 predicted cluster.
#   pred_singleton_rate     pred_singleton_count / n_sequences_merged.
#
# --- Family F: Cluster-size diagnostics (in cluster-level) ---
#
#   Per predicted cluster:
#     n_members             Number of sequences in this predicted cluster.
#     n_true_clusters       Number of distinct true sets represented here.
#     dominant_true_cluster The true set with the largest count.
#     dominant_count        Count of members in the dominant true set.
#     purity                dominant_count / n_members.
#     n_pairs               C(n_members, 2), total pairs in the cluster.
#     n_tp_pairs            Sum over true sets s of C(count_s, 2).
#     n_fp_pairs            n_pairs - n_tp_pairs.
#     pair_precision        n_tp_pairs / n_pairs.
#
#   Sort by n_fp_pairs descending to find exactly which clusters
#   are causing low pooled precision.
#
# --- Family G: Input counts ---
#
#   n_sequences_input       Unique sequence_ids in the cluster file.
#   n_sequences_gt          Unique sequence_ids in the ground truth.
#   n_sequences_merged      Sequences present in both (intersection).
#   n_clusters_input        Number of distinct predicted cluster IDs.
#   n_clusters              Number of predicted clusters in merged set.
#
# --- Family H: Pairwise agreement between two cluster files (optional) ---
#
#   CCS     Cluster Correspondence Score: fraction of all pairs that
#           two clustering runs agree on.
#   CCR     Cluster Count Ratio: n_clusters_2 / n_clusters_1.
#   RSS     Representative Sequence Stability: fraction of clusters
#           whose representative sequence is preserved in the other
#           run when overlap >= 50%.
#
# --- Family I: Taxonomic purity ---
#
#   taxonomic_agreement     Mean over predicted clusters of the
#                           fraction of members whose source_category
#                           (viral/prokaryotic/eukaryotic) matches the
#                           cluster's majority category.
#
# ============================================================

import argparse
import sys
import os
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
            filepath, header=None, sep=sep_char, quotechar='"',
            skipinitialspace=True, engine="python", on_bad_lines="warn",
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
        if (col0 in rep_keywords or any(kw in col0 for kw in rep_keywords)) and \
           (col1 in set_keywords or any(kw in col1 for kw in set_keywords)):
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
    for _, row in df.iterrows():
        cluster = str(row["cluster_id"]).strip()
        if not cluster or cluster.lower() == "nan":
            skipped += 1
            continue
        members = str(row["cluster_set"]).split(",")
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
    f1 = (2 * precision * recall / (precision + recall)
          if (precision + recall) > 0 else 0.0)
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
        rows.append({
            "cluster_id": cid,
            "n_members": n,
            "n_true_clusters": n_true,
            "dominant_true_cluster": dominant,
            "dominant_count": dom_count,
            "purity": round(purity, 6),
            "n_pairs": n_pairs,
            "n_tp_pairs": tp_pairs,
            "n_fp_pairs": fp_pairs,
            "pair_precision": (round(pair_prec, 6) if not pd.isna(pair_prec) else np.nan),
        })
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
    if merged_df.empty:
        return {
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
        return {
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

    recalls_arr = np.array(recalls)
    n_sets = int(len(recalls_arr))
    n_fully = int((recalls_arr >= 1.0 - 1e-12).sum())
    n_ge_90 = int((recalls_arr >= 0.90).sum())
    n_ge_80 = int((recalls_arr >= 0.80).sum())
    n_ge_70 = int((recalls_arr >= 0.70).sum())
    n_ge_60 = int((recalls_arr >= 0.60).sum())
    n_ge_50 = int((recalls_arr >= 0.50).sum())
    n_lt_50 = int((recalls_arr <  0.50).sum())

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
    mapping = {}
    if merged_df.empty:
        return mapping
    for cid, group in merged_df.groupby("cluster_id"):
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
    n = len(labels1)
    if n < 2:
        return np.nan
    total_pairs = n * (n - 1) // 2
    TP, _, _, _ = pair_counts(labels1, labels2)
    return TP / total_pairs if total_pairs > 0 else 0.0


def ccr(labels1, labels2):
    n1 = len(set(labels1))
    n2 = len(set(labels2))
    return n2 / n1 if n1 > 0 else np.nan


def rss(cluster_df1, cluster_df2):
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
#  Aggregated metrics
# ============================================================

def cluster_external_metrics(
    cluster_df, ground_truth_df, sample_id=None, file_path=None
):
    n_sequences_input = int(cluster_df["sequence_id"].nunique())
    n_clusters_total = int(cluster_df["cluster_id"].nunique())
    n_sequences_gt = int(ground_truth_df["sequence_id"].nunique())

    merged = pd.merge(cluster_df, ground_truth_df, on="sequence_id", how="inner")
    n_sequences_merged = int(merged["sequence_id"].nunique()) if not merged.empty else 0

    if merged.empty:
        return {
            "n_sequences_input": n_sequences_input,
            "n_sequences_gt": n_sequences_gt,
            "n_sequences_merged": n_sequences_merged,
            "n_clusters_input": n_clusters_total,
        }, pd.DataFrame(), pd.DataFrame()

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
    # Macro analogue of the Rand Index: harmonic mean of macro precision and
    # macro recall. Same numeric value as adjusted_f1 but reported separately
    # as a named column to make its Rand-Index semantics explicit.
    adj_ri_macro = harmonic_mean(adj_prec, adj_rec)

    pset = per_set_metrics(merged)

    mapping = map_clusters_to_truth(merged)

    gt_sizes = ground_truth_df.groupby("true_cluster").size()
    gt_singletons = set(gt_sizes[gt_sizes == 1].index)
    pred_sizes = merged.groupby("cluster_id").size()
    pred_singletons = set(pred_sizes[pred_sizes == 1].index)

    gt_set_sizes = None
    if "set_id" in ground_truth_df.columns:
        gt_set_sizes = ground_truth_df.groupby("set_id").size()

    cluster_size_map = merged.groupby("cluster_id").size().to_dict()
    cluster_purity_map = {}
    cluster_n_true_map = {}
    for cid, group in merged.groupby("cluster_id"):
        counts = group["true_cluster"].value_counts()
        cluster_purity_map[cid] = float(counts.iloc[0] / counts.sum())
        cluster_n_true_map[cid] = int(len(counts))

    merged["cluster_assigned"] = merged["cluster_id"]
    merged["cluster_assigned_mapped"] = merged["cluster_id"].map(mapping)
    merged["match"] = merged["cluster_assigned_mapped"] == merged["true_cluster"]
    merged["gt_is_singleton"] = merged["true_cluster"].isin(gt_singletons)
    merged["pred_is_singleton"] = merged["cluster_id"].isin(pred_singletons)
    merged["pred_cluster_size"] = merged["cluster_id"].map(cluster_size_map).fillna(0).astype(int)
    merged["pred_cluster_purity"] = merged["cluster_id"].map(cluster_purity_map).fillna(0.0)
    merged["pred_cluster_n_true"] = merged["cluster_id"].map(cluster_n_true_map).fillna(0).astype(int)

    if "set_id" in merged.columns:
        if gt_set_sizes is not None:
            merged["set_size_gt"] = (
                merged["set_id"].map(gt_set_sizes).fillna(0).astype(int)
            )
        else:
            merged["set_size_gt"] = (
                merged.groupby("set_id")["sequence_id"].transform("size")
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
        description="Compare clustering outputs against ground truth.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--ground-truth", required=True)
    parser.add_argument("--cluster-files", nargs="+", required=True)
    parser.add_argument("--sample-names", nargs="*", default=None)
    parser.add_argument("--output-dir", default=".")
    parser.add_argument("--sep", choices=["auto", "tab", "comma"], default="auto")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-n", "--dry-run", action="store_true")
    args = parser.parse_args()

    if not os.path.isfile(args.ground_truth):
        console.print(f"[bold red]ERROR:[/] Ground truth not found: {args.ground_truth}")
        sys.exit(1)
    for f in args.cluster_files:
        if not os.path.isfile(f):
            console.print(f"[bold red]ERROR:[/] Cluster file not found: {f}")
            sys.exit(1)

    console.rule("[bold green]Cluster Comparison Tool")
    show_parameters(args)
    if args.dry_run:
        console.print("[bold yellow]DRY RUN - no files will be written.[/]")

    # ---- Set up output directories ----
    merged_dir = None
    clusters_dir = None
    if not args.dry_run:
        os.makedirs(args.output_dir, exist_ok=True)
        merged_dir = os.path.join(args.output_dir, "sequence-level")
        clusters_dir = os.path.join(args.output_dir, "cluster-level")
        os.makedirs(merged_dir, exist_ok=True)
        os.makedirs(clusters_dir, exist_ok=True)
        console.log(f"[cyan]Per-sample merged GT files   -> {merged_dir}[/]")
        console.log(f"[cyan]Per-sample cluster summaries -> {clusters_dir}[/]")

    console.print("\n[bold]Loading data...[/]")
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
                df = parse_cluster_file(f, sep=args.sep, verbose=args.verbose)
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

    if len(cluster_dfs) == 0:
        console.print("[bold red]No valid cluster files loaded. Exiting.[/]")
        sys.exit(1)
    if len(cluster_dfs) == 1:
        console.print(
            "[cyan]Single cluster file detected - computing external metrics only "
            "(pairwise agreement is skipped).[/]"
        )

    console.print("\n[bold]Computing external metrics...[/]")
    external_summary = []
    merged_tables = []
    per_cluster_tables = {}

    with Progress(
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TimeElapsedColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("[cyan]External metrics", total=len(cluster_dfs))
        for sid, df, fpath in zip(sample_ids, cluster_dfs, file_paths):
            metrics, merged, per_cluster_df = cluster_external_metrics(
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
                per_cluster_tables[sid] = per_cluster_df

                if not args.dry_run:
                    # Per-sample merged GT file -> subfolder
                    out_merged = os.path.join(
                        merged_dir, f"{sid}_merged_with_gt.tsv"
                    )
                    merged.to_csv(out_merged, sep="\t", index=False)

                    # Per-sample predicted-cluster summary -> subfolder
                    out_clusters = os.path.join(
                        clusters_dir, f"{sid}_predicted_clusters_summary.tsv"
                    )
                    per_cluster_df.to_csv(out_clusters, sep="\t", index=False)
            progress.advance(task)

    if external_summary:
        summary_df = pd.DataFrame(external_summary)
        cols = [
            "sample_id",
            "file_path",
            "n_sequences_input",
            "n_sequences_gt",
            "n_sequences_merged",
            "n_clusters_input",
            "seq_accuracy_mapped",
            "rand_index",
            "ARI",
            "adjusted_rand_index_macro",
            "adjusted_precision",
            "adjusted_recall",
            "adjusted_f1",
            "per_set_recall",
            "per_set_precision",
            "per_set_f1",
            "median_best_cluster_size",
            "mean_best_cluster_size",
            "n_sets_total",
            "n_sets_fully_retrieved",
            "n_sets_at_least_90",
            "n_sets_at_least_80",
            "n_sets_at_least_70",
            "n_sets_at_least_60",
            "n_sets_at_least_50",
            "n_sets_below_50",
            "pct_sets_fully_retrieved",
            "pct_sets_at_least_90",
            "pct_sets_at_least_80",
            "pct_sets_at_least_70",
            "pct_sets_at_least_60",
            "pct_sets_at_least_50",
            "pct_sets_below_50",
            "n_clusters",
            "precision",
            "recall",
            "f1",
            "taxonomic_agreement",
            "gt_singleton_count",
            "gt_singleton_rate",
            "pred_singleton_count",
            "pred_singleton_rate",
        ]
        cols = [c for c in cols if c in summary_df.columns]
        summary_df = summary_df[cols]

        if not args.dry_run:
            summary_path = os.path.join(args.output_dir, "external_metrics_summary.tsv")
            summary_df.to_csv(summary_path, sep="\t", index=False)
            console.log(f"[green]OK[/] Summary written to [cyan]{summary_path}[/]")
        else:
            console.log("[yellow]Dry-run: would write external_metrics_summary.tsv[/]")
    else:
        console.print("[yellow]No external metrics computed.[/]")

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
            task = progress.add_task("[cyan]Pairwise", total=total_pairs)
            for (i, df1, sid1), (j, df2, sid2) in combinations(
                zip(range(len(cluster_dfs)), cluster_dfs, sample_ids), 2
            ):
                ag = cluster_agreement_metrics(df1, df2)
                if ag:
                    pairwise_results.append({
                        "sample1": sid1, "sample2": sid2,
                        "CCS": ag["CCS"], "CCR": ag["CCR"], "RSS": ag["RSS"],
                    })
                progress.advance(task)

        if pairwise_results:
            pair_df = pd.DataFrame(pairwise_results)[["sample1", "sample2", "CCS", "CCR", "RSS"]]
            if not args.dry_run:
                pair_path = os.path.join(args.output_dir, "pairwise_agreement.tsv")
                pair_df.to_csv(pair_path, sep="\t", index=False)
                console.log(f"[green]OK[/] Pairwise written to [cyan]{pair_path}[/]")
            else:
                console.log("[yellow]Dry-run: would write pairwise_agreement.tsv[/]")
        else:
            console.print("[yellow]No pairwise metrics computed.[/]")
    else:
        console.print("[yellow]Only one cluster file - skipping pairwise metrics.[/]")

    if merged_tables:
        combined = pd.concat(merged_tables, ignore_index=True)
        if not args.dry_run:
            combined_path = os.path.join(args.output_dir, "all_merged_with_gt.tsv")
            combined.to_csv(combined_path, sep="\t", index=False)
            console.log(f"[green]OK[/] Combined merged written to [cyan]{combined_path}[/]")

    if per_cluster_tables:
        combined_clusters = []
        for sid, pcdf in per_cluster_tables.items():
            if pcdf.empty:
                continue
            tmp = pcdf.copy()
            tmp.insert(0, "sample_id", sid)
            combined_clusters.append(tmp)
        if combined_clusters:
            combined_clusters_df = pd.concat(combined_clusters, ignore_index=True)
            if not args.dry_run:
                cc_path = os.path.join(
                    args.output_dir, "all_predicted_clusters_summary.tsv"
                )
                combined_clusters_df.to_csv(cc_path, sep="\t", index=False)
                console.log(
                    f"[green]OK[/] Combined predicted-cluster summary written to "
                    f"[cyan]{cc_path}[/]"
                )

    console.rule("[bold green]Summary")
    success_count = len(external_summary)
    failure_count = len(cluster_dfs) - success_count
    table = Table(title="Results Summary", box=box.SIMPLE)
    table.add_column("Metric", style="cyan")
    table.add_column("Value", style="white")
    table.add_row("Samples processed", str(len(cluster_dfs)))
    table.add_row("Successful external metrics", f"{success_count} / {len(cluster_dfs)}")
    table.add_row(
        "Failures", str(failure_count) if failure_count > 0 else "[green]None[/]"
    )
    table.add_row(
        "Output directory",
        args.output_dir if not args.dry_run else "[yellow]Dry-run (no files written)[/]",
    )
    if not args.dry_run:
        table.add_row("Per-sample merged files", merged_dir)
        table.add_row("Per-sample cluster summaries", clusters_dir)
    console.print(table)

    if failure_count > 0:
        console.print("[yellow]Some samples failed - see messages above.[/]")
    else:
        console.print("[green]All samples processed successfully.[/]")


if __name__ == "__main__":
    main()