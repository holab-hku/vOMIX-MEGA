import csv
import sys
import argparse
from pathlib import Path

import pandas as pd


def detect_sep(path):
    with open(path, "r", encoding="utf-8", newline="") as f:
        sample = f.read(2048)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=[",", "\t", ";", "|"])
        return dialect.delimiter
    except csv.Error:
        suffix = Path(path).suffix.lower()
        if suffix == ".tsv":
            return "\t"
        if suffix == ".csv":
            return ","
        return "\t"


def main():
    parser = argparse.ArgumentParser(description="Column-bind different tables.")

    # Define our expected arguments
    parser.add_argument(
        "-i", "--inputs", nargs="+", required=True, help="List of input files"
    )
    parser.add_argument("-o", "--output", required=True, help="Output file path")
    parser.add_argument(
        "--column-prefix",
        nargs="+",
        default=[],
        help="List of prefixes (must match number of inputs)",
    )
    parser.add_argument(
        "--column-suffix",
        nargs="+",
        default=[],
        help="List of suffixes (must match number of inputs)",
    )

    args = parser.parse_args()

    # Validation: Ensure prefixes/suffixes match the number of input files
    if args.column_prefix and len(args.column_prefix) != len(args.inputs):
        sys.exit(
            f"Error: Provided {len(args.column_prefix)} prefixes for {len(args.inputs)} input files."
        )
    if args.column_suffix and len(args.column_suffix) != len(args.inputs):
        sys.exit(
            f"Error: Provided {len(args.column_suffix)} suffixes for {len(args.inputs)} input files."
        )

    dfs = []
    row_counts = []

    for idx, input_file in enumerate(args.inputs):
        # 1. Recognize the delimitation
        sep = detect_sep(input_file)

        try:
            df = pd.read_csv(input_file, sep=sep)
        except Exception as e:
            print(f"Error reading {input_file}: {e}", file=sys.stderr)
            sys.exit(1)

        # Track row counts for our warning system later
        row_counts.append((Path(input_file).name, len(df)))

        # 2. Apply chronological prefixes and suffixes if provided
        prefix = args.column_prefix[idx] if args.column_prefix else ""
        suffix = args.column_suffix[idx] if args.column_suffix else ""

        if prefix or suffix:
            # Rename columns dynamically
            df.columns = [f"{prefix}{col}{suffix}" for col in df.columns]

        dfs.append(df)

    # 3. Warn about row mismatches
    unique_counts = set([count for name, count in row_counts])
    if len(unique_counts) > 1:
        warning_msg = "\n[!] Warning: Row Count Mismatch Detected!\n"
        for name, count in row_counts:
            warning_msg += f"  -> {name}: {count} rows\n"
        warning_msg += "Pandas will safely column-bind, and missing cells will be filled with empty values (NaN).\n"

        print(warning_msg, file=sys.stderr)

    # Column bind them (axis=1 means columns, as opposed to axis=0 for rows)
    merged_df = pd.concat(dfs, axis=1)

    # Write to output safely
    merged_df.to_csv(args.output, sep="\t", index=False)
    print(f"Successfully column-bound {len(args.inputs)} files into {args.output}")


if __name__ == "__main__":
    main()
