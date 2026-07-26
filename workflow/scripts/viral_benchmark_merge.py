import argparse
import pandas as pd
import os
from functools import reduce


def clean_contig_id(cid):
    """
    Strips spaces and pipes to isolate the raw contig ID.
    Converts: 'k141_12863 flag=0 multi=9' -> 'k141_12863'
    Converts: 'k141_11585||full'          -> 'k141_11585'
    """
    return str(cid).split(" ")[0].split("||")[0]


def load_and_prep(file_path, id_column, prefix):
    """
    Loads a TSV/TXT file, cleans the ID, adds prefixes, and sets the index.
    """
    if not file_path or not os.path.exists(file_path):
        return None

    try:
        # All tools output tab-delimited files, even DVF's .txt
        df = pd.read_csv(file_path, sep="\t")
    except pd.errors.EmptyDataError:
        return None

    if id_column not in df.columns:
        return None

    # Drop VirFinder's unnamed artifact column if it exists
    if "Unnamed: 0" in df.columns:
        df = df.drop(columns=["Unnamed: 0"])

    # Clean IDs to ensure seamless row matching
    df[id_column] = df[id_column].apply(clean_contig_id)

    # Deduplicate in case any tool outputs multiple lines for the same contig
    df = df.drop_duplicates(subset=[id_column])

    # Set the cleaned ID as the index
    df = df.set_index(id_column)

    # Add software prefix to the remaining columns
    df.columns = [f"{prefix}_{col}" for col in df.columns]

    # Standardize index name
    df.index.name = "contig_id"
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Merge outputs from various viral identification tools."
    )
    parser.add_argument("--genomad")
    parser.add_argument("--dvf")
    parser.add_argument("--phamer")
    parser.add_argument("--phavip")
    parser.add_argument("--virsorter2")
    parser.add_argument("--virfinder")
    parser.add_argument("--vibrant")
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    dataframes = []

    # Map the tools to their specific sequence ID column headers
    tool_mappings = [
        (args.genomad, "seq_name", "genomad"),
        (args.dvf, "name", "dvf"),
        (args.phamer, "Accession", "phamer"),
        (args.phavip, "Accession", "phavip"),
        (args.virsorter2, "seqname", "virsorter2"),
        (args.virfinder, "name", "virfinder"),
        (args.vibrant, "scaffold", "vibrant"),
    ]

    for file_path, id_col, prefix in tool_mappings:
        df = load_and_prep(file_path, id_col, prefix)
        if df is not None:
            dataframes.append(df)

    if not dataframes:
        print("No valid input files found.")
        # Create an empty file to satisfy Snakemake
        pd.DataFrame().to_csv(args.output, sep="\t")
        return

    # Outer join all dataframes on the index ('contig_id')
    merged_df = reduce(
        lambda left, right: pd.merge(
            left, right, left_index=True, right_index=True, how="outer"
        ),
        dataframes,
    )

    merged_df.to_csv(args.output, sep="\t")


if __name__ == "__main__":
    main()
