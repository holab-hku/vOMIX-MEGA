#!/usr/bin/env python
"""Utility functions for vOMIX-MEGA Snakemake workflow."""

import os
import sys
import datetime
import json
from typing import List, Union

from rich.console import Console
from rich.panel import Panel

# ----------------------------------------------------------------------
# Global objects (set from the Snakefile)
# ----------------------------------------------------------------------
console = Console()
config = None  # set to the Snakemake config dict
outdir = None  # set to the output directory
targets = []  # global list of target logfiles
current_module = None


def _log(message: str) -> None:
    """Log a message if config['verbose'] is True."""
    if config is not None and config.get("verbose", False):
        console.log(message)


# ----------------------------------------------------------------------
# Path helpers
# ----------------------------------------------------------------------
def cleanpath(input_path: str) -> str:
    """Remove trailing slashes from a directory path."""
    _log(f"Cleaning path: {input_path}")
    if input_path.endswith(os.sep):
        return input_path.rstrip(os.sep)
    return input_path


def relpath(input_path: str) -> str:
    """Join output directory with a relative path."""
    if outdir is None:
        raise RuntimeError("outdir is not set. Set vomix_utils.outdir first.")
    result = os.path.join(outdir, input_path)
    _log(f"Relative path '{input_path}' -> '{result}'")
    return result


def fullpath(input_path: str) -> str:
    """Expand `~` and convert to absolute path."""
    result = os.path.abspath(os.path.expanduser(input_path))
    _log(f"Full path '{input_path}' -> '{result}'")
    return result


# ----------------------------------------------------------------------
# FASTA input helpers
# ----------------------------------------------------------------------
def readfasta(f: str) -> str:
    """Validate and return a FASTA file path, setting sample-name if needed."""
    _log(f"Reading fasta file: {f}")
    _, extension = os.path.splitext(f)
    console.print(
        f"\n[dim]The config['fasta'] parameter is not empty, using '{f}' as input."
    )
    if extension.lower() not in [".fa", ".fasta", ".fna"]:
        console.print(
            Panel.fit(
                "File path does not end with .fa, .fasta, or .fna",
                title="Error",
                subtitle="Input not fasta file",
            )
        )
        sys.exit(1)

    if not os.path.exists(f):
        console.print(
            Panel.fit(
                "[dim]The fasta file path provided does not exist.",
                title="Error",
                subtitle="Fasta File Path",
            )
        )
        sys.exit(1)

    if config.get("sample-name", "") == "":
        sample_id = os.path.splitext(os.path.basename(f))[0]
        console.print(
            Panel.fit(
                f"[dim]config['sample-name'] not provided, using base name {sample_id} for output naming.",
                title="Warning",
                subtitle="Sample Name",
            )
        )
        config["sample-name"] = sample_id
        _log(f"Set sample-name to '{sample_id}'")
    else:
        sample_id = config["sample-name"]
        config["sample-name"] = sample_id
    return f


def readfastadir(f: str) -> List[str]:
    """Read all FASTA files from a directory and return their full paths."""
    _log(f"Reading fasta directory: {f}")
    fastadir = cleanpath(f)
    if not os.path.exists(fastadir):
        console.print(
            Panel.fit(
                f"The input file path '{fastadir}' does not exist.",
                title="Error",
                subtitle="Fasta Directory Not Found.",
            )
        )
        sys.exit(1)

    fastafs = [
        fname
        for fname in os.listdir(fastadir)
        if fname.endswith((".fa", ".fasta", ".fna"))
    ]
    if len(fastafs) == 0:
        console.print(
            Panel.fit(
                f"There are no files ending with .fa, .fasta, or .fna in '{fastadir}'.",
                title="Error",
                subtitle="No Fasta Files",
            )
        )
        sys.exit(1)

    assembly_ids = [os.path.basename(fname).rsplit(".", 1)[0] for fname in fastafs]
    fastafs_str = "\n".join(fastafs)
    console.print(
        Panel.fit(
            f"Notice:[dim] Reading input as fasta files from config['fastadir'], "
            f"the following files have been parsed correctly.\n{fastafs_str}",
            title="Notice",
            subtitle="Reading Fasta Directory",
        )
    )
    config["assembly-ids"] = assembly_ids
    _log(f"Detected {len(fastafs)} fasta files; assembly IDs: {assembly_ids}")
    return [os.path.join(fastadir, fname) for fname in fastafs]


# ----------------------------------------------------------------------
# Configuration preparation
# ----------------------------------------------------------------------
def update_db_paths(input_list: list) -> None:
    """Make all database paths absolute using the workflow basedir."""
    _log(f"Updating database paths for keys: {input_list}")
    for key in input_list:
        config[key] = os.path.join(config["basedir"], config[key])
        _log(f"  {key} -> {config[key]}")


def clean_config_paths(configkeys: list) -> None:
    """Convert the specified config keys to clean, absolute paths."""
    _log(f"Cleaning config paths for keys: {configkeys}")
    for configkey in configkeys:
        if configkey in config and config[configkey] != "":
            full_path = fullpath(cleanpath(config[configkey]))
            config[configkey] = full_path
            _log(f"  {configkey} -> {full_path}")


def setup_output_dir(outdir: str) -> None:
    """Create the output directory and its .vomix subdirectory."""
    _log(f"Setting up output directory: {outdir}")
    if not (os.path.exists(outdir) and os.path.exists(os.path.join(outdir, ".vomix"))):
        os.makedirs(os.path.join(outdir, ".vomix"), exist_ok=True)
        _log("  Created .vomix subdirectory")


def validate_input_count(config: dict) -> None:
    """Exit if the number of provided input types is not exactly one."""
    _log("Validating input count")
    inputkeys = ["fasta", "samplelist", "fastadir"]
    emptykeys = sum(1 for key in inputkeys if config.get(key, "") != "")
    if emptykeys != 1 and (config["module"] not in ["setup-database", "containerize"]):
        console.print(
            Panel.fit(
                "[dim] Multiple or no inputs detected. You can EITHER provide config['samplelist'] "
                "OR config['fastadir'] OR config['fasta']. You must provide at least one input format. "
                "Please read more about the permitted inputs on our wiki page.",
                title="Error",
                subtitle="Multiple Input Formats",
            )
        )
        sys.exit(1)
    _log(f"  Number of input types: {emptykeys}")


def save_config_file(config: dict, outdir: str) -> None:
    """Save the current config as a JSON file under .vomix/log/vomix<timestamp>."""
    _log("Saving configuration file")
    now = datetime.datetime.now()
    nowstr = now.strftime("%Y%m%d_%H%M%S")
    config["latest-run"] = nowstr

    logdir = os.path.join(outdir, ".vomix", "log", "vomix" + nowstr)
    os.makedirs(logdir, exist_ok=True)
    with open(os.path.join(logdir, "config.json"), "w") as configf:
        json.dump(config, configf)
    _log(f"  Saved to {os.path.join(logdir, 'config.json')}")


# ----------------------------------------------------------------------
# Module target management
# ----------------------------------------------------------------------
def set_module_target(logfiles: Union[str, List[str]]) -> None:
    """Add one or more logfile paths to the global targets list.

    If config['reset'] is True, the logfiles are deleted before being added.
    """
    logfile_list = [logfiles] if isinstance(logfiles, str) else logfiles
    _log(f"Adding targets: {logfile_list}")

    if config.get("reset", False):
        module_name = config.get("module", "")
        _log(f"Reset mode enabled for module '{module_name}'")
        for logfile in logfile_list:
            if not os.path.exists(logfile):
                console.print()
                console.print(
                    Panel.fit(
                        f"[bold yellow]Warning[/bold yellow]: [cyan]config['reset'] = True[/cyan] "
                        f"for module [bold]{module_name}[/bold], but no done.log file "
                        f"was found at [dim]{logfile}[/dim]. No action taken.",
                        title=f"Warning: {module_name}",
                        subtitle="No Modules to Reset",
                        border_style="yellow",
                    )
                )
            else:
                try:
                    os.remove(logfile)
                    console.print()
                    console.print(
                        Panel.fit(
                            f"[bold red]Warning[/bold red]: [cyan]config['reset'] = True[/cyan] "
                            f"for module [bold]{module_name}[/bold]. The done.log file "
                            f"at [dim]{logfile}[/dim] has been successfully deleted.",
                            title=f"Warning: {module_name}",
                            subtitle="Module Reset",
                            border_style="red",
                        )
                    )
                    _log(f"  Deleted {logfile}")
                except OSError as e:
                    console.print(
                        f"[bold red]Error[/bold red]: Could not delete {logfile}. Reason: {e}"
                    )

    targets.extend(logfile_list)
