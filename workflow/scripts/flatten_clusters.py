#!/usr/bin/env python3
"""
Flatten hierarchical cluster files into a single TSV with all members.
"""

import os
import sys
import argparse
import glob
import re
from collections import defaultdict
from rich.console import Console
from rich.table import Table
from rich.panel import Panel

console = Console()


def parse_tsv_clstr(filepath, verbose=False):
    """
    Robust parser for TSV cluster files.
    Handles headers, quoted fields, comma/space separated members.
    """
    cluster_map = {}
    with open(filepath, "r") as f:
        lines = f.readlines()
        if not lines:
            return cluster_map

        # Detect header
        first = lines[0].strip().split("\t")
        is_header = False
        if len(first) >= 2:
            header_lower = [col.lower() for col in first]
            if any(
                "rep" in col or "representative" in col for col in header_lower
            ) or any("member" in col or "cluster" in col for col in header_lower):
                is_header = True
                if verbose:
                    console.log(
                        f"  Detected header in {os.path.basename(filepath)}: {first}"
                    )

        start = 1 if is_header else 0
        if verbose:
            console.log(
                f"  Parsing {len(lines)-start} data lines from {os.path.basename(filepath)}"
            )

        for idx, line in enumerate(lines[start:], start=1):
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                # If only one column, treat as singleton
                rep = parts[0].strip()
                cluster_map[rep] = [rep]
                if verbose:
                    console.log(f"    Line {idx}: singleton (only one column) -> {rep}")
                continue

            rep = parts[0].strip()
            # Combine remaining columns into a single string, then strip quotes
            rest = "\t".join(parts[1:]).strip()
            # Remove outer quotes if present
            if rest.startswith('"') and rest.endswith('"'):
                rest = rest[1:-1]
            elif rest.startswith("'") and rest.endswith("'"):
                rest = rest[1:-1]

            if verbose and idx <= 3:
                console.log(f"    Line {idx}: rep={rep}, rest={rest}")

            if not rest:
                cluster_map[rep] = [rep]
                continue

            # Split members
            # First try comma
            if "," in rest:
                members = [m.strip() for m in rest.split(",") if m.strip()]
            else:
                # Try whitespace (tabs or spaces)
                members = re.split(r"[\s,]+", rest)
                members = [m for m in members if m and m != rep]
                # If still no members, take the whole rest as one member
                if not members and rest:
                    members = [rest]

            if not members:
                members = [rep]
            # Ensure rep is first (if not already)
            if rep not in members:
                members.insert(0, rep)
            # Remove duplicate rep if present more than once
            if members.count(rep) > 1:
                members = [rep] + [m for m in members if m != rep]
            cluster_map[rep] = members

    if verbose:
        console.log(
            f"  Parsed {len(cluster_map)} clusters from {os.path.basename(filepath)}"
        )
    return cluster_map


def parse_cdhit_clstr(filepath, verbose=False):
    """Parse CD-HIT .clstr file."""
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
    return cluster_map


def parse_vsearch_uc(filepath, verbose=False):
    """Parse vsearch .uc file."""
    cluster_map = {}
    current_rep = None
    with open(filepath, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if not parts:
                continue
            record_type = parts[0]
            if record_type == "C":
                seq = parts[8]
                current_rep = seq
                cluster_map[current_rep] = [seq]
            elif record_type == "S" and current_rep is not None:
                seq = parts[8]
                if seq != current_rep:
                    cluster_map[current_rep].append(seq)
    return cluster_map


def parse_clstr(filepath, format, verbose=False):
    if format == "cdhit":
        return parse_cdhit_clstr(filepath, verbose)
    elif format == "tsv":
        return parse_tsv_clstr(filepath, verbose)
    elif format == "vsearch":
        return parse_vsearch_uc(filepath, verbose)
    else:
        raise ValueError(f"Unsupported format: {format}")


def flatten_clusters(layer_maps, final_reps, verbose=False):
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
        seen = set()
        unique = []
        for m in all_members:
            if m not in seen:
                seen.add(m)
                unique.append(m)
        if rep not in unique:
            unique.insert(0, rep)
        flat_map[rep] = unique
    return flat_map


def print_layer_summary(layer, cluster_map, verbose=True):
    if not verbose or not cluster_map:
        return
    total_seqs = sum(len(members) for members in cluster_map.values())
    num_clusters = len(cluster_map)
    avg_size = total_seqs / num_clusters if num_clusters > 0 else 0
    singletons = sum(1 for members in cluster_map.values() if len(members) == 1)
    singleton_pct = singletons / num_clusters * 100 if num_clusters > 0 else 0
    console.log(
        f"    Layer {layer}: {num_clusters} clusters, {total_seqs} sequences, "
        f"avg {avg_size:.2f}, {singletons} singletons ({singleton_pct:.1f}%)"
    )
    if singleton_pct > 90 and layer > 1:
        console.print(
            f"    [yellow]⚠ High singleton rate – parsing may be incorrect![/]"
        )


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
            console.log(f"Layer {layer}: found {len(files)} chunk files")
        if not files:
            console.print(f"[yellow]Warning: No files found for layer {layer}[/]")
            continue

        combined = {}
        for f in files:
            try:
                cmap = parse_clstr(f, args.format, verbose=args.verbose)
                if args.verbose:
                    console.log(
                        f"  Parsed {len(cmap)} clusters from {os.path.basename(f)}"
                    )
                combined.update(cmap)
            except Exception as e:
                console.print(f"[red]Error parsing {f}: {e}[/]")
                if args.verbose:
                    console.print_exception()
        if combined:
            layer_maps[layer] = combined
            print_layer_summary(layer, combined, args.verbose)

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
        console.log(f"Final layer: {len(final_reps)} representatives")
        print_layer_summary(args.cluster_iter, final_rep_map, args.verbose)

    # Flatten
    flat = flatten_clusters(layer_maps, final_reps, verbose=args.verbose)

    # Write output
    with open(args.output, "w") as out:
        for rep, members in flat.items():
            out.write(f"{rep}\t{','.join(members)}\n")

    console.print(f"[green]✓ Flattened {len(flat)} clusters to {args.output}[/]")


if __name__ == "__main__":
    main()
