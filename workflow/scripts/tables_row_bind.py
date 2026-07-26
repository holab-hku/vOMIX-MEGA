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
    parser = argparse.ArgumentParser(
        description="Row-bind multiple tables with identical columns."
    )

    # Define our argparse arguments
    parser.add_argument(
        "-i", "--inputs", nargs="+", required=True, help="List of input files"
    )
    parser.add_argument("-o", "--output", required=True, help="Output file path")
    parser.add_argument(
        "-s",
        "--sample-ids",
        nargs="*",
        help="Optional list of sample IDs corresponding to the input files",
    )

    args = parser.parse_args()

    output_file = args.output
    input_files = args.inputs
    sample_ids = args.sample_ids

    # Validate that if sample IDs are provided, they match the number of inputs
    if sample_ids and len(sample_ids) != len(input_files):
        print(
            f"Error: Number of sample IDs ({len(sample_ids)}) does not match number of input files ({len(input_files)}).",
            file=sys.stderr,
        )
        sys.exit(1)

    dfs = []
    reference_columns = None
    reference_file = None

    for idx, input_file in enumerate(input_files):
        # Recognize the delimitation of the file
        sep = detect_sep(input_file)

        try:
            df = pd.read_csv(input_file, sep=sep)
        except Exception as e:
            print(f"Error reading {input_file}: {e}", file=sys.stderr)
            sys.exit(1)

        # Check for column consistency against the first file loaded
        if reference_columns is None:
            reference_columns = set(df.columns)
            reference_file = Path(input_file).name
        else:
            current_columns = set(df.columns)
            if current_columns != reference_columns:
                # Generate a clean, descriptive error message
                missing_cols = reference_columns - current_columns
                extra_cols = current_columns - reference_columns

                error_msg = "\n[!] Column Mismatch Error!\n"
                error_msg += f"The file '{Path(input_file).name}' does not match the columns of the first file ('{reference_file}').\n"

                if missing_cols:
                    error_msg += f"  -> Missing columns: {', '.join(missing_cols)}\n"
                if extra_cols:
                    error_msg += (
                        f"  -> Unexpected extra columns: {', '.join(extra_cols)}\n"
                    )

                print(error_msg, file=sys.stderr)
                sys.exit(1)

        # Add a column to note which file it came from (base name + extension)
        df["source_file"] = Path(input_file).name

        # If sample IDs were provided, add the sample_id column
        if sample_ids:
            # We use insert at index 0 so it appears as the very first column
            df.insert(0, "sample_id", sample_ids[idx])

        dfs.append(df)

    if not dfs:
        print("No data to merge.", file=sys.stderr)
        sys.exit(1)

    # Concatenate them safely
    merged_df = pd.concat(dfs, ignore_index=True)

    # Automatically choose output separator based on file extension
    out_sep = "," if Path(output_file).suffix.lower() == ".csv" else "\t"

    # Write to output
    merged_df.to_csv(output_file, sep=out_sep, index=False)


if __name__ == "__main__":
    main()
