#!/usr/bin/env python
"""
Parse fastp and fastplong JSON reports to aggregate statistics.

Usage:
    fastp_parse.py --names names.txt --jsons jsons.txt > stats.csv

Where names.txt and jsons.txt contain space‑separated sample IDs and JSON paths.
"""

import argparse
import sys
import pandas as pd
import json


def parse_json_files(names, csvs):
    """
    Read each JSON, extract before/after filtering stats,
    and produce a merged CSV.
    """
    df = pd.DataFrame(
        columns=[
            "total_reads_before_filtering",
            "total_reads_after_filtering",
            "total_bases_before_filtering",
            "total_bases_after_filtering",
            "gc_content_before_filtering",
            "gc_content_after_filtering",
        ],
        index=names,
    )

    for sample_id, csv in zip(names, csvs):
        try:
            with open(csv, "r") as jsonf:
                data = json.load(jsonf)
        except Exception as e:
            print(f"Warning: Could not read {csv}: {e}", file=sys.stderr)
            continue

        # Skip dummy long‑read placeholders
        if data.get("longread"):
            continue

        # Parse standard fastp/fastplong structure
        try:
            before = data["summary"]["before_filtering"]
            after = data["summary"]["after_filtering"]
        except KeyError:
            print(
                f"Warning: JSON missing 'summary' key in {csv}, skipping.",
                file=sys.stderr,
            )
            continue

        df.at[sample_id, "total_reads_before_filtering"] = before.get("total_reads", 0)
        df.at[sample_id, "total_bases_before_filtering"] = before.get("total_bases", 0)
        df.at[sample_id, "gc_content_before_filtering"] = before.get("gc_content", 0)
        df.at[sample_id, "total_reads_after_filtering"] = after.get("total_reads", 0)
        df.at[sample_id, "total_bases_after_filtering"] = after.get("total_bases", 0)
        df.at[sample_id, "gc_content_after_filtering"] = after.get("gc_content", 0)

    # Output to STDOUT as CSV
    print(df.to_csv(index=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge fastp/fastplong JSON reports")
    parser.add_argument(
        "--names", type=str, help="File containing space‑separated sample IDs"
    )
    parser.add_argument(
        "--jsons", type=str, help="File containing space‑separated JSON file paths"
    )
    args = parser.parse_args()

    with open(args.names, "r") as f:
        names = f.read().split()
    with open(args.jsons, "r") as f:
        json_files = f.read().split()

    parse_json_files(names, json_files)
