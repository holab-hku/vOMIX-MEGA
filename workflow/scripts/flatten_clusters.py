#!/usr/bin/env python3
"""
Flatten hierarchical cluster files into a single TSV with all members.

Inputs:
  --final-clstr  : Path to the final .clstr file (at layer cluster_iter, chunk 0)
  --layer-dir    : Base directory containing all layers
  --cluster-iter : Number of layers (final layer number)
  --output       : Output TSV file (representative, member_list)
  --format       : Format of .clstr files (cdhit, tsv, vsearch). Default: cdhit.
  --verbose, -v  : Enable verbose logging

Output:
  TSV with columns: representative, member_list (comma‑separated)
"""

import os
import sys
import argparse
import glob
import re
from collections import defaultdict

# Rich imports
try:
    from rich.console import Console
    from rich.panel import Panel
except ImportError:
    # Fallback if rich is not installed (should not happen in pipeline)
    class Console:
        def log(self, *args, **kwargs):
            print(*args)

        def print(self, *args, **kwargs):
            print(*args)

    console = Console()
else:
    console = Console()

# ----------------------------------------------------------------------
# Parsers for different cluster file formats
# ----------------------------------------------------------------------


def parse_cdhit_clstr(filepath, verbose=False):
    """
    Parse CD-HIT .clstr file into dict: rep -> list of members.
    """
    cluster_map = {}
    with open(filepath, "r") as f:
        current_rep = None
        current_members = []
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                if current_rep is not None:
                    cluster_map[current_rep] = current_members
                current_rep = None
                current_members = []
            elif line:
                match = re.search(r">(\S+)", line)
                if match:
                    seq = match.group(1)
                    if "*" in line:
                        current_rep = seq
                    current_members.append(seq)
        if current_rep is not None:
            cluster_map[current_rep] = current_members
    if verbose:
        console.log(f"Parsed {len(cluster_map)} clusters from {filepath}")
    return cluster_map


def parse_tsv_clstr(filepath, verbose=False):
    """
    Parse TSV cluster file (rep, member_list) into dict: rep -> list of members.
    """
    cluster_map = {}
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                rep = parts[0].strip()
                members = parts[1].strip()
                if members.startswith('"') and members.endswith('"'):
                    members = members[1:-1]
                member_list = [m.strip() for m in members.split(",") if m.strip()]
                cluster_map[rep] = member_list
            else:
                # Single sequence cluster
                cluster_map[parts[0].strip()] = [parts[0].strip()]
    if verbose:
        console.log(f"Parsed {len(cluster_map)} clusters from {filepath}")
    return cluster_map


def parse_vsearch_uc(filepath, verbose=False):
    """
    Parse vsearch .uc file into dict: rep -> list of members.
    """
    cluster_map = {}
    current_rep = None
    with open(filepath, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if not parts:
                continue
            record_type = parts[0]
            if record_type == "C":  # centroid / representative
                seq = parts[8]
                current_rep = seq
                cluster_map[current_rep] = [seq]
            elif record_type == "S" and current_rep is not None:
                seq = parts[8]
                if seq != current_rep:
                    cluster_map[current_rep].append(seq)
    if verbose:
        console.log(f"Parsed {len(cluster_map)} clusters from {filepath}")
    return cluster_map


def parse_clstr(filepath, format, verbose=False):
    """
    Dispatch to the appropriate parser.
    """
    if format == "cdhit":
        return parse_cdhit_clstr(filepath, verbose)
    elif format == "tsv":
        return parse_tsv_clstr(filepath, verbose)
    elif format == "vsearch":
        return parse_vsearch_uc(filepath, verbose)
    else:
        raise ValueError(f"Unsupported format: {format}")


# ----------------------------------------------------------------------
# Recursive flattening
# ----------------------------------------------------------------------


def flatten_clusters(layer_maps, final_reps, verbose=False):
    """
    Recursively expand representatives to all original sequences.
    """
    max_layer = max(layer_maps.keys()) if layer_maps else 0
    if verbose:
        console.log(f"Max layer: {max_layer}")

    def expand(seq, layer):
        if layer == 0:
            return [seq]
        if seq not in layer_maps[layer]:
            return [seq]
        members = layer_maps[layer][seq]
        expanded = []
        for m in members:
            if m == seq:
                continue
            expanded.extend(expand(m, layer - 1))
        return expanded

    flat_map = {}
    for rep in final_reps:
        all_members = expand(rep, max_layer)
        # Deduplicate while preserving order
        seen = set()
        unique_members = []
        for m in all_members:
            if m not in seen:
                seen.add(m)
                unique_members.append(m)
        if rep not in unique_members:
            unique_members.insert(0, rep)
        flat_map[rep] = unique_members
    return flat_map


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description="Flatten hierarchical cluster files.")
    parser.add_argument(
        "--final-clstr",
        required=True,
        help="Final cluster file (layer cluster_iter, chunk 0)",
    )
    parser.add_argument(
        "--layer-dir", required=True, help="Base directory containing all layers"
    )
    parser.add_argument(
        "--cluster-iter",
        type=int,
        required=True,
        help="Number of layers (final layer number)",
    )
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument(
        "--format",
        default="cdhit",
        choices=["cdhit", "tsv", "vsearch"],
        help="Format of .clstr files (default: cdhit)",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose logging"
    )
    args = parser.parse_args()

    console.log(f"[cyan]Flattening clusters with format {args.format}[/]")
    if args.verbose:
        console.log(f"Final cluster file: {args.final_clstr}")
        console.log(f"Layer directory: {args.layer_dir}")
        console.log(f"Cluster iterations: {args.cluster_iter}")
        console.log(f"Output file: {args.output}")

    # Collect all cluster files per layer
    layer_maps = {}
    for layer in range(1, args.cluster_iter + 1):
        pattern = os.path.join(args.layer_dir, f"layer_{layer}", "chunk_*.fa.clstr")
        files = glob.glob(pattern)
        if args.verbose:
            console.log(f"Layer {layer}: found {len(files)} files")
        if not files:
            console.print(
                f"[yellow]Warning: No cluster files found for layer {layer} in {args.layer_dir}[/]"
            )
            continue
        combined = {}
        for f in files:
            try:
                cmap = parse_clstr(f, args.format, verbose=args.verbose)
                combined.update(cmap)
            except Exception as e:
                console.print(f"[red]Error parsing {f}: {e}[/]")
                if args.verbose:
                    console.print_exception()
        if combined:
            layer_maps[layer] = combined
            if args.verbose:
                console.log(f"Layer {layer}: {len(combined)} unique clusters")

    if not layer_maps:
        console.print("[red]Error: No cluster files found in any layer.[/]")
        sys.exit(1)

    # Parse final cluster file to get final representatives
    try:
        final_rep_map = parse_clstr(args.final_clstr, args.format, verbose=args.verbose)
    except Exception as e:
        console.print(
            f"[red]Error parsing final cluster file {args.final_clstr}: {e}[/]"
        )
        if args.verbose:
            console.print_exception()
        sys.exit(1)

    if not final_rep_map:
        console.print("[red]Error: Final cluster file is empty or invalid.[/]")
        sys.exit(1)

    final_reps = list(final_rep_map.keys())
    if args.verbose:
        console.log(f"Found {len(final_reps)} final representatives")

    # Flatten
    flat = flatten_clusters(layer_maps, final_reps, verbose=args.verbose)

    # Write output
    with open(args.output, "w") as out:
        for rep, members in flat.items():
            out.write(f"{rep}\t{','.join(members)}\n")

    console.print(f"[green]✓ Flattened {len(flat)} clusters to {args.output}[/]")


if __name__ == "__main__":
    main()
