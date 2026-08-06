import rich_click as click
import sys
import logging
import os
import platform
import time
from rich.logging import RichHandler
from rich.console import Console

from vomix.snakemakeFlags import SnakemakeFlags
from vomix.vomix_actions import vomix_actions
from vomix.modules import (
    PreProcessingModule,
    AssemblyCoAssemblyModule,
    ViralIdentifyModule,
    ViralBenchmarkModule,
    ViralTaxonomyModule,
    ViralHostModule,
    ViralCommunityModule,
    ViralAnnotateModule,
    ProkaryoticCommunityModule,
    ProkaryoticBinningModule,
    ProkaryoticAnnotateModule,
    ViralEndToEndModule,
    ClusterFastModule,
    CheckVPyHMMERModule,
    SetupDatabaseModule,
)

# ---------------------------------------------------------
# Rich Logging Configuration
# ---------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, markup=True, show_path=False)],
)

log = logging.getLogger("rich")
console = Console()

modules_list = [
    "assembly",
    "checkv-pyhmmer",
    "checkv",
    "clustering-fast",
    "clustering-sensitive",
    "host-cherry",
    "host",
    "preprocess",
    "prok-annotate",
    "prok-binning",
    "prok-community",
    "refilter-genomad",
    "setup-database",
    "symlink",
    "viral-annotate",
    "viral-benchmark",
    "viral-binning",
    "viral-community",
    "viral-host",
    "viral-identify",
    "viral-refilter",
    "viral-taxonomy",
]

END_MODULE_RUN_LOG = "Module (or dry-run) execution completed successfully."

# ---------------------------------------------------------
# CLI Help Formatting Configurations (Rich Click)
# ---------------------------------------------------------

click.rich_click.USE_RICH_MARKUP = True
click.rich_click.USE_MARKDOWN = True
click.rich_click.SHOW_ARGUMENTS = True
click.rich_click.GROUP_ARGUMENTS_OPTIONS = True

click.rich_click.COMMAND_GROUPS = {
    "vomix": [
        {
            "name": "1. Data Preparation & Assembly",
            "commands": ["preprocess", "assembly"],
        },
        {
            "name": "2. Viral Analytics",
            "commands": [
                "viral-identify",
                "viral-benchmark",
                "viral-taxonomy",
                "viral-host",
                "viral-community",
                "viral-annotate",
                "viral-end-to-end",
            ],
        },
        {
            "name": "3. Prokaryotic Analytics",
            "commands": ["prok-binning", "prok-annotate", "prok-community"],
        },
        {
            "name": "4. Utilities & Setup",
            "commands": ["cluster-fast", "checkv-pyhmmer", "setup-database"],
        },
    ]
}

_common_opts = [
    "--workdir",
    "--outdir",
    "--datadir",
    "--samplelist",
    "--fasta",
    "--fastadir",
    "--sample-name",
    "--assembly-ids",
    "--latest-run",
    "--splits",
    "--keep-intermediates",
    "--setup-database",
    "--max-cores",
    "--NCBI-email",
    "--NCBI-API-key",
    "--custom-config",
]
_smk_opts = [
    "--dry-run",
    "--forceall",
    "--configfile",
    "--unlock",
    "--cores",
    "--jobs",
    "--latency-wait",
    "--rerun-incomplete",
    "--rerun-triggers",
    "--software-deployment-method",
    "--executor",
    "--cluster-generic-submit-cmd",
    "--printshellcmds",
    "--quiet",
    "--snakemake-args",
    "--dryrun",
]

click.rich_click.OPTION_GROUPS = {}
for cmd in modules_list + ["cluster-fast", "viral-end-to-end"]:
    click.rich_click.OPTION_GROUPS[f"vomix {cmd}"] = [
        {"name": "Core Pipeline Options", "options": _common_opts},
        {"name": "Snakemake Backend Options", "options": _smk_opts},
    ]

# ---------------------------------------------------------
# Dynamic Helper Functions for Efficient Execution
# ---------------------------------------------------------


def log_system_info(module_name):
    """Logs system and execution environment info for robust provenance."""
    log.info("vOMIX-MEGA initialized.")
    log.info(f"Starting Module: {module_name}")
    log.info(f"System: {platform.system()} {platform.release()} ({platform.machine()})")
    log.info(f"Python: {platform.python_version()}")


def setOptions(module_obj, kwargs):
    """Maps common Click options to the module object using kwargs."""
    module_obj.workdir = kwargs.get("workdir")
    module_obj.outdir = kwargs.get("outdir")
    module_obj.datadir = kwargs.get("datadir")
    module_obj.samplelist = kwargs.get("samplelist")
    module_obj.fasta = kwargs.get("fasta")
    module_obj.fastadir = kwargs.get("fastadir")
    module_obj.sample_name = kwargs.get("sample_name")
    module_obj.assembly_ids = kwargs.get("assembly_ids")
    module_obj.latest_run = kwargs.get("latest_run")
    module_obj.keep_intermediates = kwargs.get("keep_intermediates")
    module_obj.setup_database = kwargs.get("setup_database")
    module_obj.max_cores = kwargs.get("max_cores")
    module_obj.ncbi_email = kwargs.get("ncbi_email")
    module_obj.ncbi_api_key = kwargs.get("ncbi_api_key")
    module_obj.custom_config = kwargs.get("custom_config")
    return module_obj


def create_snakemake_flags(kwargs):
    """Generates the SnakemakeFlags object automatically from kwargs."""
    # Careful that click automatically returns the longest argument as the default name
    # hence you need to do software_deployment_method rather than sdm
    return SnakemakeFlags(
        kwargs.get("dry_run"),
        kwargs.get("forceall"),
        kwargs.get("configfile"),
        kwargs.get("unlock"),
        kwargs.get("cores"),
        kwargs.get("jobs"),
        kwargs.get("latency_wait"),
        kwargs.get("rerun_incomplete"),
        kwargs.get("rerun_triggers"),
        kwargs.get("software_deployment_method"),
        kwargs.get("executor"),
        kwargs.get("cluster_generic_submit_cmd"),
        kwargs.get("printshellcmds"),
        kwargs.get("quiet"),
        kwargs.get("snakemake_args"),
    )


def apply_module_options(module_obj, kwargs, specific_opts_mapping):
    """Dynamically applies custom options to a module object and triggers hasOptions."""
    for obj_attr, kwarg_key in specific_opts_mapping.items():
        if kwargs.get(kwarg_key) is not None:
            val = kwargs.get(kwarg_key)
            setattr(module_obj, obj_attr, val)
            module_obj.hasOptions = True
            log.info(f"  ↳ Set parameter {kwarg_key} = {val}")


def common_options(function):
    function = click.option(
        "--workdir",
        default=None,
        required=False,
        help="The working directory for the underlying Snakemake workflow back-end. [default: None]",
    )(function)
    function = click.option(
        "--outdir",
        default=None,
        required=False,
        help="Directory path where structured output results will be deposited. [default: None]",
    )(function)
    function = click.option(
        "--datadir",
        default=None,
        required=False,
        help="The path to raw FASTQ files. [default: None]",
    )(function)
    function = click.option(
        "--samplelist",
        default=None,
        required=False,
        help="The path to the sample_list.csv configuration file. [default: None]",
    )(function)
    function = click.option(
        "--fasta",
        default=None,
        required=False,
        help="The path to a single input FASTA file. [default: None]",
    )(function)
    function = click.option(
        "--fastadir",
        default=None,
        required=False,
        help="The path to a directory containing input FASTA files. [default: None]",
    )(function)
    function = click.option(
        "--sample-name",
        default=None,
        required=False,
        help="Explicit sample name utilized for output file naming. [default: None]",
    )(function)
    function = click.option(
        "--assembly-ids",
        default=None,
        required=False,
        help="A list format specified array mapping sample names to input files. [default: None]",
    )(function)
    function = click.option(
        "--latest-run",
        default=None,
        required=False,
        help="Internal logging parameter designating the timestamp. [default: None]",
    )(function)
    function = click.option(
        "--keep-intermediates",
        is_flag=True,
        default=None,
        required=False,
        help="Retain substantial intermediate processing files. [default: False]",
    )(function)
    function = click.option(
        "--setup-database",
        is_flag=True,
        default=None,
        required=False,
        help="Initialize or update databases. [default: False]",
    )(function)
    function = click.option(
        "--max-cores",
        default=None,
        required=False,
        help="Max CPU cores allocated dynamically across parallel tasks. [default: None]",
    )(function)
    function = click.option(
        "--NCBI-email",
        default=None,
        required=False,
        help="User email address provided to NCBI E-utilities. [default: None]",
    )(function)
    function = click.option(
        "--NCBI-API-key",
        default=None,
        required=False,
        help="NCBI API key required for higher throughput data retrieval. [default: None]",
    )(function)
    function = click.option(
        "--custom-config",
        default=None,
        required=False,
        help="Path to your custom config.yml [default: None]",
    )(function)
    return function


def snakemake_options(function):
    function = click.option(
        "--dry-run",
        "--dryrun",
        "-n",
        is_flag=True,
        flag_value=True,
        required=False,
        default=None,
        help="Do not execute anything, display what would be done. [default: False]",
    )(function)
    function = click.option(
        "--forceall",
        "-F",
        is_flag=True,
        flag_value=True,
        required=False,
        default=None,
        help="Force the execution of the selected rule and all dependencies. [default: False]",
    )(function)
    function = click.option(
        "--configfile",
        default=None,
        required=False,
        help="Specify or overwrite the config file of the workflow. [default: None]",
    )(function)
    function = click.option(
        "--unlock",
        is_flag=True,
        flag_value=True,
        required=False,
        default=None,
        help="Remove a lock on the working directory. [default: False]",
    )(function)
    function = click.option(
        "--cores",
        "-c",
        default=None,
        required=False,
        help="Use at most N CPU cores/jobs in parallel. [default: None]",
    )(function)
    function = click.option(
        "--jobs",
        "-j",
        default=None,
        required=False,
        help="Use at most N CPU cluster/cloud jobs in parallel. [default: None]",
    )(function)
    function = click.option(
        "--latency-wait",
        default=None,
        required=False,
        help="Wait given seconds if an output file of a job is not present. [default: None]",
    )(function)
    function = click.option(
        "--rerun-incomplete",
        "-ri",
        is_flag=True,
        flag_value=True,
        required=False,
        default=None,
        help="Re-run all jobs the output of which is recognized as incomplete. [default: False]",
    )(function)
    function = click.option(
        "--rerun-triggers",
        multiple=True,
        required=False,
        default=None,
        help="Define what triggers the rerunning of a job. [default: mtime]",
    )(function)
    function = click.option(
        "--software-deployment-method",
        "--deployment-method",
        "--deployment",
        "--sdm",
        multiple=True,
        required=False,
        default="conda",
        help="Specify software environment deployment method. [default: conda]",
    )(function)
    function = click.option(
        "--executor",
        "-e",
        required=False,
        default=None,
        help="Specify a custom executor, available via an executor plugin. [default: None]",
    )(function)
    function = click.option(
        "--cluster-generic-submit-cmd",
        required=False,
        default=None,
        help="Command for submitting jobs.",
    )(function)
    function = click.option(
        "--printshellcmds",
        "-p",
        is_flag=True,
        flag_value=True,
        required=False,
        default=None,
        help="Print out the shell commands that will be executed. (default: False)",
    )(function)
    function = click.option(
        "--quiet",
        "-q",
        required=False,
        default=None,
        flag_value=True,
        help="Do not output certain information. [default: False]",
    )(function)
    function = click.option(
        "--snakemake-args",
        required=False,
        default=None,
        help='Additional arguments to pass to the native snakemake command. Must be surround by double quotes `"`. [default: None]',
    )(function)
    return function


# ---------------------------------------------------------
# vOMIX-MEGA Command Line Interface
# ---------------------------------------------------------


@click.group(name="vomix", context_settings={"help_option_names": ["-h", "--help"]})
def cli():
    """vOMIX-MEGA - an ultra-fast & modular end-to-end pipeline for terabyte-scale viral metagenomics analysis."""


@cli.command(
    "preprocess",
    context_settings={"ignore_unknown_options": True},
    short_help="Run the Pre-processing module",
)
@common_options
@click.option(
    "--dwnld-only",
    is_flag=True,
    flag_value=True,
    default=None,
    required=False,
    help="Restrict execution exclusively to remote SRA file downloads. [default: False]",
)
@click.option(
    "--decontam-host/--no-decontam-host",
    is_flag=True,
    flag_value=True,
    default=None,
    required=False,
    help="Perform host decontamination post fastp quality trimming. [default: False]",
)
@click.option(
    "--dwnld-params",
    required=False,
    default=None,
    help="Optional configuration parameters used during raw FASTQ retrieval. [default: None]",
)
@click.option(
    "--pigz-params",
    required=False,
    default=None,
    help="Execution parameters passed directly to pigz. [default: None]",
)
@click.option(
    "--fastp-params",
    required=False,
    default=None,
    help="Additional runtime arguments supplied to the fastp quality control engine. [default: None]",
)
@click.option(
    "--hostile-params",
    required=False,
    default=None,
    help="Additional runtime arguments supplied to the Hostile host decontamination module. [default: None]",
)
@click.option(
    "--hostile-aligner",
    required=False,
    default=None,
    help="The short-read alignment backend algorithm employed for host decontamination. [default: None]",
)
@click.option(
    "--hostile-aligner-params",
    required=False,
    default=None,
    help="Additional runtime arguments supplied directly to the selected Hostile alignment tool. [default: None]",
)
@click.option(
    "--hostile-index-name",
    required=False,
    default=None,
    help="The name identifier of pre-built Hostile indices. [default: None]",
)
@click.option(
    "--hostile-index-db",
    required=False,
    default=None,
    help="The directory path where the Hostile database is installed. [default: None]",
)
@snakemake_options
def run_preprocess(**kwargs):
    log_system_info("preprocess")

    module_obj = setOptions(PreProcessingModule(), kwargs)
    module_obj.name = "preprocess"

    apply_module_options(
        module_obj,
        kwargs,
        {
            "dwnld_only": "dwnld_only",
            "decontam_host": "decontam_host",
            "dwnld_params": "dwnld_params",
            "pigz_params": "pigz_params",
            "fastp_params": "fastp_params",
            "hostile_params": "hostile_params",
            "hostile_aligner": "hostile_aligner",
            "aligner_params": "hostile_aligner_params",
            "hostile_index_name": "hostile_index_name",
            "hostile_index_db": "hostile_index_db",
        },
    )

    vomix_actions().run_module("preprocess", module_obj, create_snakemake_flags(kwargs))
    log.info(END_MODULE_RUN_LOG)


@cli.command(
    "assembly",
    context_settings={"ignore_unknown_options": True},
    short_help="Run the Assembly & Co-assembly module",
)
@common_options
@click.option(
    "--assembler",
    default=None,
    required=False,
    help="The primary assembly tool engine selected for the assembly module. [default: None]",
)
@click.option(
    "--megahit-min-len",
    required=False,
    default=None,
    help="The minimum contig length threshold used during the MEGAHIT assembly filter steps. [default: 300]",
)
@click.option(
    "--megahit-params",
    required=False,
    default=None,
    help="Additional runtime parameters supplied directly to the MEGAHIT execution pipeline. [default: None]",
)
@click.option(
    "--spades-params",
    required=False,
    default=None,
    help="Additional runtime parameters supplied directly to the SPAdes metagenomic assembler execution line. [default: None]",
)
@click.option(
    "--spades-memory",
    required=False,
    default=None,
    help="The upper threshold of RAM memory (in gigabytes) allocated for SPAdes assembly execution. [default: 250]",
)
@snakemake_options
def run_assembly(**kwargs):
    log_system_info("assembly")

    module_obj = setOptions(AssemblyCoAssemblyModule(), kwargs)
    module_obj.name = "assembly"

    apply_module_options(
        module_obj,
        kwargs,
        {
            "assembler": "assembler",
            "megahit_min_len": "megahit_min_len",
            "megahit_params": "megahit_params",
            "spades_params": "spades_params",
            "spades_memory": "spades_memory",
        },
    )

    vomix_actions().run_module("assembly", module_obj, create_snakemake_flags(kwargs))
    log.info(END_MODULE_RUN_LOG)


@cli.command(
    "viral-identify",
    context_settings={"ignore_unknown_options": True},
    short_help="Run the Viral Identify module",
)
@common_options
@click.option(
    "--contig-min-len",
    required=False,
    default=None,
    help="The absolute minimum length constraint for contig inclusion. [default: 0]",
)
@click.option(
    "--contig-splits",
    required=False,
    default=None,
    help="Number of clustering and pooling iterations (L) for vOTU clustering. Initially splits data into 2^(L-1) partitions for memory efficiency. Set to 1 to disable partitioning. [default: 0]",
)
@click.option(
    "--contig-splits",
    required=False,
    default=None,
    help="Number of clustering and pooling iterations (L) for vOTU clustering. Initially splits data into 2^(L-1) partitions for memory efficiency. Set to 1 to disable partitioning. [default: 0]",
)
@click.option(
    "--genomad-db",
    required=False,
    default=None,
    help="The directory path where the geNomad database is installed. [default: None]",
)
@click.option(
    "--genomad-min-len",
    required=False,
    default=None,
    help="The minimum contig length evaluated by geNomad. [default: 1500]",
)
@click.option(
    "--genomad-params",
    required=False,
    default=None,
    help="Additional runtime command line arguments supplied to geNomad execution. [default: None]",
)
@click.option(
    "--genomad-cutoff",
    required=False,
    default=None,
    help="The minimal numeric confidence threshold required by geNomad. [default: 0.7]",
)
@click.option(
    "--genomad-cutoff-s",
    required=False,
    default=None,
    help="The minimum confidence threshold applied during geNomad secondary filtering. [default: 0]",
)
@click.option(
    "--checkv-original/--checkv-pyhmmer",
    is_flag=True,
    flag_value=True,
    required=False,
    default=None,
    help="Flag allowing execution of standard CheckV. [default: False]",
)
@click.option(
    "--checkv-splits",
    required=False,
    default=None,
    help="Number of input data splits for CheckV or CheckV-PyHMMER analysis, splitting data into n+1 chunks. Higher values reduce memory usage but increase runtime. Set to 0 for no splitting. [default: 0]",
)
@click.option(
    "--checkv-params",
    required=False,
    default=None,
    help="Additional operational arguments supplied directly to the CheckV pipeline. [default: None]",
)
@click.option(
    "--checkv-database",
    required=False,
    default=None,
    help="The directory path where the CheckV database is installed. [default: None]",
)
@click.option(
    "--clustering-fast/--clustering-cdhit",
    is_flag=True,
    flag_value=True,
    required=False,
    default=None,
    help="Flag triggering an accelerated MEGABlast-based clustering protocol. [default: True]",
)
@click.option(
    "--cluster-iter",
    required=False,
    default=None,
    help="Number of clustering and pooling iterations (L) for vOTU clustering. Initially splits data into 2^(L-1) partitions for memory efficiency. Set to 1 to disable partitioning. [default: 1]",
)
@click.option(
    "--cdhit-params",
    required=False,
    default=None,
    help="Additional operational runtime values supplied directly to the CD-HIT clustering utility. [default: None]",
)
@click.option(
    "--vOTU-ani",
    required=False,
    default=None,
    help="The average nucleotide identity (ANI) clustering percentage threshold. [default: 95]",
)
@click.option(
    "--vOTU-targetcov",
    required=False,
    default=None,
    help="The minimum target coverage alignment coverage percentage. [default: 85]",
)
@click.option(
    "--vOTU-querycov",
    required=False,
    default=None,
    help="The target query alignment coverage percentage criteria. [default: 0]",
)
@snakemake_options
def run_viral_identify(**kwargs):
    log_system_info("viral-identify")

    module_obj = setOptions(ViralIdentifyModule(), kwargs)
    module_obj.name = "viral-identify"

    apply_module_options(
        module_obj,
        kwargs,
        {
            "contig_min_len": "contig_min_len",
            "contig_splits": "contig_splits",
            "genomad_db": "genomad_db",
            "genomad_min_len": "genomad_min_len",
            "genomad_params": "genomad_params",
            "genomad_cutoff": "genomad_cutoff",
            "genomad_cutoff_s": "genomad_cutoff_s",
            "checkv_original": "checkv_original",
            "checkv_splits": "checkv_splits",
            "checkv_params": "checkv_params",
            "checkv_database": "checkv_database",
            "clustering_fast": "clustering_fast",
            "cluster_iter": "cluster_iter",
            "cdhit_params": "cdhit_params",
            "vOTU_ani": "votu_ani",
            "vOTU_targetcov": "votu_targetcov",
            "vOTU_querycov": "votu_querycov",
        },
    )

    vomix_actions().run_module(
        "viral-identify", module_obj, create_snakemake_flags(kwargs)
    )
    log.info(END_MODULE_RUN_LOG)


@cli.command(
    "viral-benchmark",
    context_settings={"ignore_unknown_options": True},
    short_help="Run the Viral Benchmark module",
)
@common_options
@click.option(
    "--PhaBox2-db",
    required=False,
    default=None,
    help="The directory path where the PhaBox2 database is installed. [default: None]",
)
@click.option(
    "--phabox2-db-name",
    required=False,
    default=None,
    help="The designated database name or identifier file package required for PhaBox2. [default: None]",
)
@click.option(
    "--phabox2-db-baselink",
    required=False,
    default=None,
    help="The primary remote server base link URL used to fetch PhaBox2. [default: None]",
)
@click.option(
    "--genomad-db",
    required=False,
    default=None,
    help="The directory path where the geNomad database is installed. [default: None]",
)
@click.option(
    "--virsorter2-db",
    required=False,
    default=None,
    help="The directory path where the VirSorter2 database is installed. [default: None]",
)
@click.option(
    "--vibrant-db",
    required=False,
    default=None,
    help="The directory path where the VIBRANT database is installed. [default: None]",
)
@click.option(
    "--contig-min-len",
    required=False,
    default=None,
    help="The absolute minimum length constraint for contig inclusion. [default: 0]",
)
@click.option(
    "--contig-splits",
    required=False,
    default=None,
    help="Number of clustering and pooling iterations (L) for vOTU clustering. Initially splits data into 2^(L-1) partitions for memory efficiency. Set to 1 to disable partitioning. [default: 0]",
)
@click.option(
    "--genomad-min-len",
    required=False,
    default=None,
    help="The minimum contig length evaluated by geNomad. [default: 1000]",
)
@click.option(
    "--genomad-params",
    required=False,
    default=None,
    help="Additional runtime command line arguments supplied to geNomad. [default: None]",
)
@click.option(
    "--genomad-cutoff",
    required=False,
    default=None,
    help="The minimal numeric confidence threshold required by geNomad. [default: 0.7]",
)
@click.option(
    "--genomad-cutoff-s",
    required=False,
    default=None,
    help="The minimum confidence threshold applied during geNomad secondary filtering. [default: 0]",
)
@click.option(
    "--dvf-min-len",
    required=False,
    default=None,
    help="The lower bound contig length cut-off implemented during DeepVirFinder evaluation. [default: 1500]",
)
@click.option(
    "--phamer-min-len",
    required=False,
    default=None,
    help="The lower bound contig length cut-off implemented during PhaMer evaluation. [default: 2000]",
)
@click.option(
    "--dvf-params",
    required=False,
    default=None,
    help="Additional system parameters passed directly to the DeepVirFinder tool. [default: None]",
)
@click.option(
    "--phamer-params",
    required=False,
    default=None,
    help="Additional system parameters passed directly to the PhaMer tool. [default: None]",
)
@click.option(
    "--virsorter2-params",
    required=False,
    default=None,
    help="Additional system parameters passed directly to the VirSorter2 tool. [default: None]",
)
@click.option(
    "--vf-params",
    required=False,
    default=None,
    help="Additional system parameters passed directly to the VirFinder tool. [default: None]",
)
@click.option(
    "--seeker-params",
    required=False,
    default=None,
    help="Additional system parameters passed directly to the Seeker tool. [default: None]",
)
@click.option(
    "--PPR-params",
    required=False,
    default=None,
    help="Additional system parameters passed directly to the PPR-META tool. [default: None]",
)
@click.option(
    "--dvf-cutoff",
    required=False,
    default=None,
    help="The minimal confidence score metric required by DeepVirFinder. [default: 0.7]",
)
@click.option(
    "--dvf-pval",
    required=False,
    default=None,
    help="The maximum critical p-value threshold permitted by DeepVirFinder. [default: 0.05]",
)
@click.option(
    "--phamer-pred",
    required=False,
    default=None,
    help="The taxonomic classification category targeted by PhaMer prediction routines. [default: None]",
)
@click.option(
    "--phamer-cutoff",
    required=False,
    default=None,
    help="The minimal confidence threshold value required for a positive viral determination within PhaMer. [default: 0]",
)
@click.option(
    "--vf-cutoff",
    required=False,
    default=None,
    help="The minimal confidence threshold value required for a positive viral determination within VirFinder. [default: 0]",
)
@click.option(
    "--virsorter2-cutoff",
    required=False,
    default=None,
    help="The minimal confidence threshold value required for a positive viral determination within VirSorter2. [default: 0]",
)
@click.option(
    "--seeker-cutoff",
    required=False,
    default=None,
    help="The minimal confidence threshold value required for a positive viral determination within Seeker. [default: 0]",
)
@click.option(
    "--ppr-cutoff",
    required=False,
    default=None,
    help="The minimal confidence threshold value required for a positive viral determination within PPR-META. [default: 0]",
)
@click.option(
    "--vibrant-cutoff",
    required=False,
    default=None,
    help="The minimal confidence threshold value required for a positive viral determination within VIBRANT. [default: 0]",
)
@click.option(
    "--checkv-original/--checkv-pyhmmer",
    is_flag=True,
    flag_value=True,
    required=False,
    default=None,
    help="Flag allowing execution of standard CheckV. [default: False]",
)
@click.option(
    "--checkv-splits",
    required=False,
    default=None,
    help="Number of input data splits for CheckV or CheckV-PyHMMER analysis, splitting data into n+1 chunks. Higher values reduce memory usage but increase runtime. Set to 0 for no splitting. [default: 0]",
)
@click.option(
    "--checkv-params",
    required=False,
    default=None,
    help="Additional operational arguments supplied directly to the CheckV pipeline execution. [default: None]",
)
@click.option(
    "--checkv-database",
    required=False,
    default=None,
    help="The directory path where the CheckV database is installed. [default: None]",
)
@click.option(
    "--clustering-fast/--clustering-cdhit",
    is_flag=True,
    flag_value=True,
    required=False,
    default=None,
    help="Flag triggering an accelerated MEGABlast-based clustering protocol. [default: True]",
)
@click.option(
    "--cluster-iter",
    required=False,
    default=None,
    help="Number of clustering and pooling iterations (L) for vOTU clustering. Initially splits data into 2^(L-1) partitions for memory efficiency. Set to 1 to disable partitioning. [default: 1]",
)
@click.option(
    "--cdhit-params",
    required=False,
    default=None,
    help="Additional operational runtime values supplied directly to the CD-HIT clustering utility. [default: None]",
)
@click.option(
    "--vOTU-ani",
    required=False,
    default=None,
    help="The average nucleotide identity (ANI) clustering percentage threshold. [default: 95]",
)
@click.option(
    "--vOTU-targetcov",
    required=False,
    default=None,
    help="The minimum target coverage alignment coverage percentage. [default: 85]",
)
@click.option(
    "--vOTU-querycov",
    required=False,
    default=None,
    help="The target query alignment coverage percentage criteria. [default: 0]",
)
@snakemake_options
def run_viral_benchmark(**kwargs):
    log_system_info("viral-benchmark")

    module_obj = setOptions(ViralBenchmarkModule(), kwargs)
    module_obj.name = "viral-benchmark"

    apply_module_options(
        module_obj,
        kwargs,
        {
            "contig_min_len": "contig_min_len",
            "contig_splits": "contig_splits",
            "PhaBox2_db": "phabox2_db",
            "phabox2_db_name": "phabox2_db_name",
            "phabox2_db_baselink": "phabox2_db_baselink",
            "genomad_db": "genomad_db",
            "virsorter2_db": "virsorter2_db",
            "vibrant_db": "vibrant_db",
            "genomad_min_len": "genomad_min_len",
            "genomad_params": "genomad_params",
            "genomad_cutoff": "genomad_cutoff",
            "genomad_cutoff_s": "genomad_cutoff_s",
            "dvf_min_len": "dvf_min_len",
            "phamer_min_len": "phamer_min_len",
            "dvf_params": "dvf_params",
            "phamer_params": "phamer_params",
            "virsorter2_params": "virsorter2_params",
            "vf_params": "vf_params",
            "seeker_params": "seeker_params",
            "ppr_params": "ppr_params",
            "dvf_cutoff": "dvf_cutoff",
            "dvf_pval": "dvf_pval",
            "phamer_pred": "phamer_pred",
            "phamer_cutoff": "phamer_cutoff",
            "vf_cutoff": "vf_cutoff",
        },
    )

    vomix_actions().run_module(
        "viral-benchmark", module_obj, create_snakemake_flags(kwargs)
    )
    log.info(END_MODULE_RUN_LOG)


@cli.command(
    "viral-taxonomy",
    context_settings={"ignore_unknown_options": True},
    short_help="Run the Viral Taxonomy module",
)
@common_options
@click.option(
    "--genomad-db",
    required=False,
    default=None,
    help="Path to geNomad database directory. [default: None]",
)
@click.option(
    "--PhaBox2-db",
    required=False,
    default=None,
    help="The directory path where the PhaBox2 database is installed. [default: None]",
)
@click.option(
    "--phabox2-db-name",
    required=False,
    default=None,
    help="The designated database name or identifier file package required for PhaBox2. [default: None]",
)
@click.option(
    "--phabox2-db-baselink",
    required=False,
    default=None,
    help="The primary remote server base link URL used to fetch PhaBox2. [default: None]",
)
@click.option(
    "--phagcn-min-len",
    required=False,
    default=None,
    help="The minimum allowed contig length for evaluation using PhaGCN. [default: 1500]",
)
@click.option(
    "--phagcn-params",
    required=False,
    default=None,
    help="Additional operational arguments passed to the PhaGCN classification instance. [default: None]",
)
@click.option(
    "--genomad-params-tax",
    required=False,
    default=None,
    help="Additional operational configurations passed to geNomad during viral taxonomic assignment. [default: None]",
)
@snakemake_options
def run_viral_taxonomy(**kwargs):
    log_system_info("viral-taxonomy")

    module_obj = setOptions(ViralTaxonomyModule(), kwargs)
    module_obj.name = "viral-taxonomy"

    apply_module_options(
        module_obj,
        kwargs,
        {
            "PhaBox2_db": "phabox2_db",
            "phabox2_db_name": "phabox2_db_name",
            "phabox2_db_baselink": "phabox2_db_baselink",
            "genomad_db": "genomad_db",
            "phagcn_min_len": "phagcn_min_len",
            "phagcn_params": "phagcn_params",
            "genomad_params_tax": "genomad_params_tax",
        },
    )

    vomix_actions().run_module(
        "viral-taxonomy", module_obj, create_snakemake_flags(kwargs)
    )
    log.info(END_MODULE_RUN_LOG)


@cli.command(
    "viral-host",
    context_settings={"ignore_unknown_options": True},
    short_help="Run the Viral Host module",
)
@common_options
@click.option(
    "--PhaBox2-db",
    required=False,
    default=None,
    help="The directory path where the PhaBox2 database is installed. [default: None]",
)
@click.option(
    "--phabox2-db-name",
    required=False,
    default=None,
    help="The designated database name or identifier file package required for PhaBox2. [default: None]",
)
@click.option(
    "--phabox2-db-baselink",
    required=False,
    default=None,
    help="The primary remote server base link URL used to fetch PhaBox2. [default: None]",
)
@click.option(
    "--CHERRY-params",
    required=False,
    default=None,
    help="Additional execution parameters configured for the CHERRY host prediction algorithm. [default: None]",
)
@click.option(
    "--PhaTYP-params",
    required=False,
    default=None,
    help="Additional custom parameters passed directly to the PhaTYP lifestyle prediction module. [default: None]",
)
@click.option(
    "--iphop-host/--cherry-host",
    is_flag=True,
    flag_value=True,
    required=False,
    default=None,
    help="Flag indicating whether to perform iPHoP-based viral host prediction instead of CHERRY. [default: False]",
)
@click.option(
    "--iphop-cutoff",
    required=False,
    default=None,
    help="The minimum confidence threshold required by iPHoP to assign a host classification profile. [default: 90]",
)
@click.option(
    "--iphop-params",
    required=False,
    default=None,
    help="Additional configuration arguments supplied directly to the iPHoP platform interface. [default: None]",
)
@click.option(
    "--iphop-db",
    required=False,
    default=None,
    help="The directory path where the iPHoP database is installed. [default: None]",
)
@click.option(
    "--iphop-db-version",
    required=False,
    default=None,
    help="The version identifier for the iPHoP database. [default: None]",
)
@click.option(
    "--iphop-db-basename",
    required=False,
    default=None,
    help="The primary base name of the iPHoP database. [default: None]",
)
@snakemake_options
def run_viral_host(**kwargs):
    log_system_info("viral-host")

    module_obj = setOptions(ViralHostModule(), kwargs)
    module_obj.name = "viral-host"

    apply_module_options(
        module_obj,
        kwargs,
        {
            "PhaBox2_db": "phabox2_db",
            "phabox2_db_name": "phabox2_db_name",
            "phabox2_db_baselink": "phabox2_db_baselink",
            "CHERRY_params": "cherry_params",
            "PhaTYP_params": "phatyp_params",
            "iphop_host": "iphop_host",
            "iphop_cutoff": "iphop_cutoff",
            "iphop_db": "iphop_db",
            "iphop_db_version": "iphop_db_version",
            "iphop_db_basename": "iphop_db_basename",
            "iphop_params": "iphop_params",
        },
    )

    vomix_actions().run_module("viral-host", module_obj, create_snakemake_flags(kwargs))
    log.info(END_MODULE_RUN_LOG)


@cli.command(
    "viral-community",
    context_settings={"ignore_unknown_options": True},
    short_help="Run the Viral Community module",
)
@common_options
@click.option(
    "--coverm-params",
    required=False,
    default=None,
    help="Additional mapping or calculation flags passed to the CoverM coverage engine. [default: None]",
)
@click.option(
    "--coverm-methods",
    required=False,
    default=None,
    help="The calculation metric outputs selected for CoverM. [default: None]",
)
@snakemake_options
def run_viral_community(**kwargs):
    log_system_info("viral-community")

    module_obj = setOptions(ViralCommunityModule(), kwargs)
    module_obj.name = "viral-community"

    apply_module_options(
        module_obj,
        kwargs,
        {"coverm_params": "coverm_params", "coverm_methods": "coverm_methods"},
    )

    vomix_actions().run_module(
        "viral-community", module_obj, create_snakemake_flags(kwargs)
    )
    log.info(END_MODULE_RUN_LOG)


@cli.command(
    "viral-annotate",
    context_settings={"ignore_unknown_options": True},
    short_help="Run the Viral Annotate module",
)
@common_options
@click.option(
    "--PhaBox2-db",
    required=False,
    default=None,
    help="The directory path where the PhaBox2 database is installed. [default: None]",
)
@click.option(
    "--phabox2-db-name",
    required=False,
    default=None,
    help="The designated database name or identifier file package required for PhaBox2. [default: None]",
)
@click.option(
    "--phabox2-db-baselink",
    required=False,
    default=None,
    help="The primary remote server base link URL used to fetch PhaBox2. [default: None]",
)
@click.option(
    "--eggNOG-params",
    required=False,
    default=None,
    help="Parameters for running eggNOG-mapper v2. [default: None]",
)
@click.option(
    "--PhaVIP-params",
    required=False,
    default=None,
    help="Minimum contig length to filter BEFORE viral identification. [default: None]",
)
@click.option(
    "--metacerberus-db",
    required=False,
    default=None,
    help="The directory path where the MetaCerberus database is installed. [default: None]",
)
@click.option(
    "--metacerberus-setup-params",
    required=False,
    default=None,
    help="Operational configurations supplied to initialize build or index the MetaCerberus database. [default: None]",
)
@click.option(
    "--metacerberus-params",
    required=False,
    default=None,
    help="Parameters for running the MetaCerberus database. [default: None]",
)
@click.option(
    "--pharokka-db",
    required=False,
    default=None,
    help="The directory path where the pharokka database is installed. [default: None]",
)
@click.option(
    "--pharokka-params",
    required=False,
    default=None,
    help="Additional execution parameters passed directly to the pharokka bacteriophage annotation framework. [default: None]",
)
@snakemake_options
def run_viral_annotate(**kwargs):
    log_system_info("viral-annotate")

    module_obj = setOptions(ViralAnnotateModule(), kwargs)
    module_obj.name = "viral-annotate"

    apply_module_options(
        module_obj,
        kwargs,
        {
            "PhaBox2_db": "phabox2_db",
            "phabox2_db_name": "phabox2_db_name",
            "phabox2_db_baselink": "phabox2_db_baselink",
            "eggNOG_params": "eggnog_params",
            "PhaVIP_params": "phavip_params",
            "metacerberus_db": "metacerberus_db",
            "metacerberus_setup_params": "metacerberus_setup_params",
            "metacerberus_params": "metacerberus_params",
            "pharokka_db": "pharokka_db",
            "pharokka_params": "pharokka_params",
        },
    )

    vomix_actions().run_module(
        "viral-annotate", module_obj, create_snakemake_flags(kwargs)
    )
    log.info(END_MODULE_RUN_LOG)


@cli.command(
    "prok-community",
    context_settings={"ignore_unknown_options": True},
    short_help="Run the Prokaryotic Community module",
)
@common_options
@click.option(
    "--mpa-params",
    required=False,
    default=None,
    help="Parameters for metaphlan function. [default: None]",
)
@click.option(
    "--mpa-indexv",
    required=False,
    default=None,
    help="Database version for metaphlan to use. [default: None]",
)
@snakemake_options
def run_prok_community(**kwargs):
    log_system_info("prok-community")

    module_obj = setOptions(ProkaryoticCommunityModule(), kwargs)
    module_obj.name = "prok-community"

    apply_module_options(
        module_obj, kwargs, {"mpa_params": "mpa_params", "mpa_indexv": "mpa_indexv"}
    )

    vomix_actions().run_module(
        "prok-community", module_obj, create_snakemake_flags(kwargs)
    )
    log.info(END_MODULE_RUN_LOG)


@cli.command(
    "prok-binning",
    context_settings={"ignore_unknown_options": True},
    short_help="Run the Prokaryotic Binning module",
)
@common_options
@click.option(
    "--checkm2-db",
    required=False,
    default=None,
    help="The directory path where the CheckM2 database is installed. [default: None]",
)
@click.option(
    "--GTDBTk-db",
    required=False,
    default=None,
    help="The directory path where the GTDB-Tk database is installed. [default: None]",
)
@click.option(
    "--GTDBTk-db-version",
    required=False,
    default=None,
    type=str,
    help="The reference version of the GTDB-Tk database. [default: None]",
)
@click.option(
    "--GTDBTk-identify-params",
    required=False,
    default=None,
    help="Additional parameter variables supplied directly to the GTDB-Tk identify command execution. [default: None]",
)
@click.option(
    "--GTDBTk-align-params",
    required=False,
    default=None,
    help="Additional parameter variables supplied directly to the GTDB-Tk align command execution. [default: None]",
)
@click.option(
    "--GTDBTk-classify-params",
    required=False,
    default=None,
    help="Additional parameter variables supplied directly to the GTDB-Tk classify command execution. [default: None]",
)
@click.option(
    "--VAMB-params",
    required=False,
    default=None,
    help="Additional parameter variables supplied directly to the VAMB command execution. [default: None]",
)
@click.option(
    "--binning-consensus/--binning-gpu",
    is_flag=True,
    flag_value=True,
    required=False,
    default=None,
    help="Flag enabling a consensus-based metagenomic binning protocol. [default: True]",
)
@click.option(
    "--strobealign-params",
    required=False,
    default=None,
    help="Additional alignment flags or scoring rules passed to the strobealign tool backend. [default: None]",
)
@click.option(
    "--MetaBAT2-params",
    required=False,
    default=None,
    help="Additional parameters for the MetaBAT2 tool. [default: None]",
)
@click.option(
    "--MaxBin2-params",
    required=False,
    default=None,
    help="Additional parameters for the MaxBin2 tool. [default: None]",
)
@click.option(
    "--CONCOCT-params",
    required=False,
    default=None,
    help="Additional parameters for the CONCOCT tool. [default: None]",
)
@click.option(
    "--jgi-summarize-params",
    required=False,
    default=None,
    help="Additional runtime parameters supplied to the jgi_summarize_bam_contig_depth command. [default: None]",
)
@click.option(
    "--DASTool-params",
    required=False,
    default=None,
    help="Additional parameters for the DASTool tool. [default: None]",
)
@click.option(
    "--checkm2-params",
    required=False,
    default=None,
    help="Additional parameter fields provided directly to the CheckM2 bin validation pipeline. [default: None]",
)
@click.option(
    "--galah-params",
    required=False,
    default=None,
    help="Additional parameters for the Galah tool. [default: None]",
)
@snakemake_options
def run_prok_binning(**kwargs):
    log_system_info("prok-binning")

    module_obj = setOptions(ProkaryoticBinningModule(), kwargs)
    module_obj.name = "prok-binning"

    apply_module_options(
        module_obj,
        kwargs,
        {
            "checkm2_db": "checkm2_db",
            "GTDBTk_db": "gtdbtk_db",
            "GTDBTk_db_version": "gtdbtk_db_version",
            "GTDBTk_identify_params": "gtdbtk_identify_params",
            "GTDBTk_align_params": "gtdbtk_align_params",
            "GTDBTk_classify_params": "gtdbtk_classify_params",
            "VAMB_params": "vamb_params",
            "binning_consensus": "binning_consensus",
            "strobealign_params": "strobealign_params",
            "MetaBAT2_params": "metabat2_params",
            "MaxBin2_params": "maxbin2_params",
            "CONCOCT_params": "concoct_params",
            "jgi_summarize_params": "jgi_summarize_params",
            "DASTool_params": "dastool_params",
            "checkm2_params": "checkm2_params",
            "galah_params": "galah_params",
        },
    )

    vomix_actions().run_module(
        "prok-binning", module_obj, create_snakemake_flags(kwargs)
    )
    log.info(END_MODULE_RUN_LOG)


@cli.command(
    "prok-annotate",
    context_settings={"ignore_unknown_options": True},
    short_help="Run the Prokaryotic Annotate module",
)
@common_options
@click.option(
    "--humann-params",
    required=False,
    default=None,
    help="Additional software parameters directed to the HUMAnN functional annotation pipeline. [default: None]",
)
@snakemake_options
def run_prok_annotate(**kwargs):
    log_system_info("prok-annotate")

    module_obj = setOptions(ProkaryoticAnnotateModule(), kwargs)
    module_obj.name = "prok-annotate"

    if kwargs.get("humann_params") is not None:
        module_obj.humann_params = kwargs.get("humann_params")
        log.info(f"  ↳ Set parameter humann_params = {module_obj.humann_params}")

    vomix_actions().run_module(
        "prok-annotate", module_obj, create_snakemake_flags(kwargs)
    )
    log.info(END_MODULE_RUN_LOG)


@cli.command(
    "viral-end-to-end",
    context_settings={"ignore_unknown_options": True},
    short_help="Run the Viral End-To-End module",
)
@common_options
@click.option(
    "--dwnld-only",
    default=None,
    required=False,
    help="Restrict execution exclusively to remote SRA file downloads. [default: False]",
)
@click.option(
    "--decontam-host/--no-decontam-host",
    is_flag=True,
    flag_value=True,
    default=None,
    required=False,
    help="Perform host decontamination post fastp quality trimming. [default: False]",
)
@click.option(
    "--dwnld-params",
    required=False,
    default=None,
    help="Optional configuration parameters used during raw FASTQ retrieval. [default: None]",
)
@click.option(
    "--pigz-params",
    required=False,
    default=None,
    help="Execution parameters passed directly to pigz. [default: None]",
)
@click.option(
    "--fastp-params",
    required=False,
    default=None,
    help="Additional runtime arguments supplied to the fastp quality control engine. [default: None]",
)
@click.option(
    "--hostile-params",
    required=False,
    default=None,
    help="Additional runtime arguments supplied to the Hostile host decontamination module. [default: None]",
)
@click.option(
    "--hostile-aligner",
    required=False,
    default=None,
    help="The short-read alignment backend algorithm employed for host decontamination. [default: None]",
)
@click.option(
    "--hostile-aligner-params",
    required=False,
    default=None,
    help="Additional runtime arguments supplied directly to the selected Hostile alignment tool. [default: None]",
)
@click.option(
    "--hostile-index-name",
    required=False,
    default=None,
    help="The name identifier of pre-built Hostile indices. [default: None]",
)
@click.option(
    "--hostile-index-db",
    required=False,
    default=None,
    help="The directory path where the Hostile database is installed. [default: None]",
)
@click.option(
    "--assembler",
    default=None,
    required=False,
    help="The primary assembly tool engine selected for the assembly module. [default: None]",
)
@click.option(
    "--megahit-min-len",
    required=False,
    default=None,
    help="The minimum contig length threshold used during the MEGAHIT assembly filter steps. [default: 300]",
)
@click.option(
    "--megahit-params",
    required=False,
    default=None,
    help="Additional runtime parameters supplied directly to the MEGAHIT execution pipeline. [default: None]",
)
@click.option(
    "--spades-params",
    required=False,
    default=None,
    help="Additional runtime parameters supplied directly to the SPAdes metagenomic assembler execution line. [default: None]",
)
@click.option(
    "--spades-memory",
    required=False,
    default=None,
    help="The upper threshold of RAM memory (in gigabytes) allocated for SPAdes assembly execution. [default: 250]",
)
@click.option(
    "--contig-min-len",
    required=False,
    default=None,
    help="The absolute minimum length constraint for contig inclusion. [default: 0]",
)
@click.option(
    "--contig-splits",
    required=False,
    default=None,
    help="Number of clustering and pooling iterations (L) for vOTU clustering. Initially splits data into 2^(L-1) partitions for memory efficiency. Set to 1 to disable partitioning. [default: 0]",
)
@click.option(
    "--contig-splits",
    required=False,
    default=None,
    help="Number of clustering and pooling iterations (L) for vOTU clustering. Initially splits data into 2^(L-1) partitions for memory efficiency. Set to 1 to disable partitioning. [default: 0]",
)
@click.option(
    "--genomad-db",
    required=False,
    default=None,
    help="The directory path where the geNomad database is installed. [default: None]",
)
@click.option(
    "--genomad-min-len",
    required=False,
    default=None,
    help="The minimum contig length evaluated by geNomad. [default: 1500]",
)
@click.option(
    "--genomad-params",
    required=False,
    default=None,
    help="Additional runtime command line arguments supplied to geNomad execution. [default: None]",
)
@click.option(
    "--genomad-cutoff",
    required=False,
    default=None,
    help="The minimal numeric confidence threshold required by geNomad. [default: 0.7]",
)
@click.option(
    "--genomad-cutoff-s",
    required=False,
    default=None,
    help="The minimum confidence threshold applied during geNomad secondary filtering. [default: 0]",
)
@click.option(
    "--checkv-original/--checkv-pyhmmer",
    is_flag=True,
    flag_value=True,
    required=False,
    default=None,
    help="Flag allowing execution of standard CheckV. [default: False]",
)
@click.option(
    "--checkv-splits",
    required=False,
    default=None,
    help="Number of input data splits for CheckV or CheckV-PyHMMER analysis, splitting data into n+1 chunks. Higher values reduce memory usage but increase runtime. Set to 0 for no splitting. [default: 0]",
)
@click.option(
    "--checkv-params",
    required=False,
    default=None,
    help="Additional operational arguments supplied directly to the CheckV pipeline. [default: None]",
)
@click.option(
    "--checkv-database",
    required=False,
    default=None,
    help="The directory path where the CheckV database is installed. [default: None]",
)
@click.option(
    "--clustering-fast/--clustering-cdhit",
    is_flag=True,
    flag_value=True,
    required=False,
    default=None,
    help="Flag triggering an accelerated MEGABlast-based clustering protocol. [default: True]",
)
@click.option(
    "--cluster-iter",
    required=False,
    default=None,
    help="Number of clustering and pooling iterations (L) for vOTU clustering. Initially splits data into 2^(L-1) partitions for memory efficiency. Set to 1 to disable partitioning. [default: 1]",
)
@click.option(
    "--cdhit-params",
    required=False,
    default=None,
    help="Additional operational runtime values supplied directly to the CD-HIT clustering utility. [default: None]",
)
@click.option(
    "--vOTU-ani",
    required=False,
    default=None,
    help="The average nucleotide identity (ANI) clustering percentage threshold. [default: 95]",
)
@click.option(
    "--vOTU-targetcov",
    required=False,
    default=None,
    help="The minimum target coverage alignment coverage percentage. [default: 85]",
)
@click.option(
    "--vOTU-querycov",
    required=False,
    default=None,
    help="The target query alignment coverage percentage criteria. [default: 0]",
)
@click.option(
    "--PhaBox2-db",
    required=False,
    default=None,
    help="The directory path where the PhaBox2 database is installed. [default: None]",
)
@click.option(
    "--phabox2-db-name",
    required=False,
    default=None,
    help="The designated database name or identifier file package required for PhaBox2. [default: None]",
)
@click.option(
    "--phabox2-db-baselink",
    required=False,
    default=None,
    help="The primary remote server base link URL used to fetch PhaBox2. [default: None]",
)
@click.option(
    "--phagcn-min-len",
    required=False,
    default=None,
    help="The minimum allowed contig length for evaluation using PhaGCN. [default: 1500]",
)
@click.option(
    "--phagcn-params",
    required=False,
    default=None,
    help="Additional operational arguments passed to the PhaGCN classification instance. [default: None]",
)
@click.option(
    "--genomad-params-tax",
    required=False,
    default=None,
    help="Additional operational configurations passed to geNomad during viral taxonomic assignment. [default: None]",
)
@click.option(
    "--CHERRY-params",
    required=False,
    default=None,
    help="Additional execution parameters configured for the CHERRY host prediction algorithm. [default: None]",
)
@click.option(
    "--PhaTYP-params",
    required=False,
    default=None,
    help="Additional custom parameters passed directly to the PhaTYP lifestyle prediction module. [default: None]",
)
@click.option(
    "--iphop-host/--cherry-host",
    is_flag=True,
    flag_value=True,
    required=False,
    default=None,
    help="Flag indicating whether to perform iPHoP-based viral host prediction instead of CHERRY. [default: False]",
)
@click.option(
    "--iphop-cutoff",
    required=False,
    default=None,
    help="The minimum confidence threshold required by iPHoP to assign a host classification profile. [default: 90]",
)
@click.option(
    "--iphop-params",
    required=False,
    default=None,
    help="Additional configuration arguments supplied directly to the iPHoP platform interface. [default: None]",
)
@click.option(
    "--iphop-db",
    required=False,
    default=None,
    help="The directory path where the iPHoP database is installed. [default: None]",
)
@click.option(
    "--eggNOG-params",
    required=False,
    default=None,
    help="Parameters for running eggNOG-mapper v2. [default: None]",
)
@click.option(
    "--PhaVIP-params",
    required=False,
    default=None,
    help="Minimum contig length to filter BEFORE viral identification. [default: None]",
)
@click.option(
    "--metacerberus-db",
    required=False,
    default=None,
    help="The directory path where the MetaCerberus database is installed. [default: None]",
)
@click.option(
    "--metacerberus-setup-params",
    required=False,
    default=None,
    help="Operational configurations supplied to initialize build or index the MetaCerberus database. [default: None]",
)
@click.option(
    "--metacerberus-params",
    required=False,
    default=None,
    help="Parameters for running the MetaCerberus database. [default: None]",
)
@click.option(
    "--pharokka-db",
    required=False,
    default=None,
    help="The directory path where the pharokka database is installed. [default: None]",
)
@click.option(
    "--pharokka-params",
    required=False,
    default=None,
    help="Additional execution parameters passed directly to the pharokka bacteriophage annotation framework. [default: None]",
)
@click.option(
    "--eggNOG-params",
    required=False,
    default=None,
    help="Parameters for running eggNOG-mapper v2. [default: None]",
)
@click.option(
    "--PhaVIP-params",
    required=False,
    default=None,
    help="Minimum contig length to filter BEFORE viral identification. [default: None]",
)
@click.option(
    "--metacerberus-db",
    required=False,
    default=None,
    help="The directory path where the MetaCerberus database is installed. [default: None]",
)
@click.option(
    "--metacerberus-setup-params",
    required=False,
    default=None,
    help="Operational configurations supplied to initialize build or index the MetaCerberus database. [default: None]",
)
@click.option(
    "--metacerberus-params",
    required=False,
    default=None,
    help="Parameters for running the MetaCerberus database. [default: None]",
)
@click.option(
    "--pharokka-db",
    required=False,
    default=None,
    help="The directory path where the pharokka database is installed. [default: None]",
)
@click.option(
    "--pharokka-params",
    required=False,
    default=None,
    help="Additional execution parameters passed directly to the pharokka bacteriophage annotation framework. [default: None]",
)
@click.option(
    "--coverm-params",
    required=False,
    default=None,
    help="Additional mapping or calculation flags passed to the CoverM coverage engine. [default: None]",
)
@click.option(
    "--coverm-methods",
    required=False,
    default=None,
    help="The calculation metric outputs selected for CoverM. [default: None]",
)
@snakemake_options
def run_viral_end_to_end(**kwargs):
    log_system_info("viral-end-to-end")

    module_obj = setOptions(ViralEndToEndModule(), kwargs)
    module_obj.name = "viral-end-to-end"

    apply_module_options(
        module_obj,
        kwargs,
        {
            "dwnld_only": "dwnld_only",
            "decontam_host": "decontam_host",
            "dwnld_params": "dwnld_params",
            "pigz_params": "pigz_params",
            "fastp_params": "fastp_params",
            "hostile_params": "hostile_params",
            "hostile_aligner": "hostile_aligner",
            "aligner_params": "hostile_aligner_params",
            "hostile_index_name": "hostile_index_name",
            "hostile_index_db": "hostile_index_db",
            "assembler": "assembler",
            "megahit_min_len": "megahit_min_len",
            "megahit_params": "megahit_params",
            "spades_params": "spades_params",
            "spades_memory": "spades_memory",
            "contig_min_len": "contig_min_len",
            "contig_splits": "contig_splits",
            "genomad_db": "genomad_db",
            "genomad_min_len": "genomad_min_len",
            "genomad_params": "genomad_params",
            "genomad_cutoff": "genomad_cutoff",
            "genomad_cutoff_s": "genomad_cutoff_s",
            "checkv_original": "checkv_original",
            "checkv_splits": "checkv_splits",
            "checkv_params": "checkv_params",
            "checkv_database": "checkv_database",
            "clustering_fast": "clustering_fast",
            "cluster_iter": "cluster_iter",
            "cdhit_params": "cdhit_params",
            "vOTU_ani": "votu_ani",
            "vOTU_targetcov": "votu_targetcov",
            "vOTU_querycov": "votu_querycov",
            "phagcn_min_len": "phagcn_min_len",
            "phagcn_params": "phagcn_params",
            "genomad_params_tax": "genomad_params_tax",
            "eggNOG_params": "eggnog_params",
            "PhaVIP_params": "phavip_params",
            "metacerberus_db": "metacerberus_db",
            "metacerberus_setup_params": "metacerberus_setup_params",
            "metacerberus_params": "metacerberus_params",
            "pharokka_db": "pharokka_db",
            "pharokka_params": "pharokka_params",
            "CHERRY_params": "cherry_params",
            "PhaTYP_params": "phatyp_params",
            "iphop_host": "iphop_host",
            "iphop_cutoff": "iphop_cutoff",
            "iphop_db": "iphop_db",
            "iphop_db_version": "iphop_db_version",
            "iphop_db_basename": "iphop_db_basename",
            "iphop_params": "iphop_params",
            "coverm_params": "coverm_params",
            "coverm_methods": "coverm_methods",
        },
    )

    vomix_actions().run_module(
        "viral-end-to-end", module_obj, create_snakemake_flags(kwargs)
    )
    log.info(END_MODULE_RUN_LOG)


@cli.command(
    "cluster-fast",
    context_settings={"ignore_unknown_options": True},
    short_help="Run the Cluster Fast module",
)
@common_options
@click.option(
    "--clustering-fast/--clustering-cdhit",
    is_flag=True,
    flag_value=True,
    required=False,
    default=None,
    help="Flag to run fast clustering using CheckV's MEGABLAST approach. [default: True]",
)
@click.option(
    "--cluster-iter",
    required=False,
    default=None,
    help="Number of clustering and pooling iterations (L) for vOTU clustering. Initially splits data into 2^(L-1) partitions for memory efficiency. Set to 1 to disable partitioning. [default: 1]",
)
@click.option(
    "--cdhit-params",
    required=False,
    default=None,
    help="Additional parameters to pass on to CD-HIT if clustering-fast is set to False. [default: None]",
)
@click.option(
    "--vOTU-ani",
    required=False,
    default=None,
    help="Minimum average nucleotide identity for fast clustering algorithm of viral contigs. [default: 95]",
)
@click.option(
    "--vOTU-targetcov",
    required=False,
    default=None,
    help="Minimum target coverage for fast clustering algorithm of viral contigs. [default: 85]",
)
@click.option(
    "--vOTU-querycov",
    required=False,
    default=None,
    help="Minimum query coverage for fast clustering algorithm of viral contigs. [default: 0]",
)
@snakemake_options
def run_cluster_fast(**kwargs):
    log_system_info("cluster-fast")

    module_obj = setOptions(ClusterFastModule(), kwargs)
    module_obj.name = "cluster-fast"

    apply_module_options(
        module_obj,
        kwargs,
        {
            "clustering_fast": "clustering_fast",
            "cluster_iter": "cluster_iter",
            "cdhit_params": "cdhit_params",
            "vOTU_ani": "votu_ani",
            "vOTU_targetcov": "votu_targetcov",
            "vOTU_querycov": "votu_querycov",
        },
    )

    vomix_actions().run_module(
        "cluster-fast", module_obj, create_snakemake_flags(kwargs)
    )
    log.info(END_MODULE_RUN_LOG)


@cli.command(
    "checkv-pyhmmer",
    context_settings={"ignore_unknown_options": True},
    short_help="Run the CheckV PyHMMER module",
)
@common_options
@click.option(
    "--checkv-original/--checkv-pyhmmer",
    is_flag=True,
    flag_value=True,
    required=False,
    default=None,
    help="Flag to use CheckV original instead of CheckV-PyHMMER. [default: False]",
)
@click.option(
    "--checkv-splits",
    required=False,
    default=None,
    help="Number of input data splits for CheckV or CheckV-PyHMMER analysis, splitting data into n+1 chunks. Higher values reduce memory usage but increase runtime. Set to 0 for no splitting. [default: 0]",
)
@click.option(
    "--checkv-params",
    required=False,
    default=None,
    help="Additional parameters to pass on to CheckV. [default: None]",
)
@click.option(
    "--checkv-database",
    required=False,
    default=None,
    help="Path to CheckV database. [default: None]",
)
@snakemake_options
def run_checkv_pyhmmer(**kwargs):
    log_system_info("checkv-pyhmmer")

    module_obj = setOptions(CheckVPyHMMERModule(), kwargs)
    module_obj.name = "checkv-pyhmmer"

    apply_module_options(
        module_obj,
        kwargs,
        {
            "checkv_original": "checkv_original",
            "checkv_splits": "checkv_splits",
            "checkv_params": "checkv_params",
            "checkv_database": "checkv_database",
        },
    )

    vomix_actions().run_module(
        "checkv-pyhmmer", module_obj, create_snakemake_flags(kwargs)
    )
    log.info(END_MODULE_RUN_LOG)


@cli.command(
    "setup-database",
    context_settings={"ignore_unknown_options": True},
    short_help="Run the Setup Database module",
)
@common_options
@click.option(
    "--hostile-index-db",
    required=False,
    default=None,
    help="The directory path where the Hostile database is installed. [default: None]",
)
@click.option(
    "--PhaBox2-db",
    required=False,
    default=None,
    help="Path to PhaBox2 database for download. [default: None]",
)
@click.option(
    "--phabox2-db-name",
    required=False,
    default=None,
    help="The designated database name or identifier file package required for PhaBox2. [default: None]",
)
@click.option(
    "--phabox2-db-baselink",
    required=False,
    default=None,
    help="The primary remote server base link URL used to fetch PhaBox2. [default: None]",
)
@click.option(
    "--genomad-db",
    required=False,
    default=None,
    help="Path to geNomad database for download. [default: None]",
)
@click.option(
    "--virsorter2-db",
    required=False,
    default=None,
    help="The directory path where the VirSorter2 database is installed. [default: None]",
)
@click.option(
    "--vibrant-db",
    required=False,
    default=None,
    help="The directory path where the VIBRANT database is installed. [default: None]",
)
@click.option(
    "--checkv-db",
    required=False,
    default=None,
    help="Path to CheckV database for download. [default: None]",
)
@click.option(
    "--eggNOG-db",
    required=False,
    default=None,
    help="Path to eggNOG v2 database for download. [default: None]",
)
@click.option(
    "--eggNOG-db-params",
    required=False,
    default=None,
    help="Parameters for downloading eggNOG v2 database. [default: None]",
)
@click.option(
    "--checkm2-db",
    required=False,
    default=None,
    help="The directory path where the CheckM2 database is installed. [default: None]",
)
@click.option(
    "--GTDBTk-db",
    required=False,
    default=None,
    help="The directory path where the GTDB-Tk database is installed. [default: None]",
)
@click.option(
    "--GTDBTk-db-version",
    required=False,
    default=None,
    type=str,
    help="The reference version of the GTDB-Tk database. [default: None]",
)
@click.option(
    "--iphop-db",
    required=False,
    default=None,
    help="Path to iPHoP database for download. [default: None]",
)
@click.option(
    "--iphop-db-version",
    required=False,
    default=None,
    help="The version identifier for the iPHoP database. [default: None]",
)
@click.option(
    "--iphop-db-basename",
    required=False,
    default=None,
    help="The primary base name of the iPHoP database. [default: None]",
)
@click.option(
    "--humann-db",
    required=False,
    default=None,
    help="Path to HUMAnN3 databases for download. [default: None]",
)
@click.option(
    "--metacerberus-db",
    required=False,
    default=None,
    help="The directory path where the MetaCerberus database is installed. [default: None]",
)
@click.option(
    "--metacerberus-setup-params",
    required=False,
    default=None,
    help="Operational configurations supplied to initialize build or index the MetaCerberus database. [default: None]",
)
@click.option(
    "--pharokka-db",
    required=False,
    default=None,
    help="The directory path where the pharokka database is installed. [default: None]",
)
@snakemake_options
def run_setup_database(**kwargs):
    log_system_info("setup-database")

    module_obj = setOptions(SetupDatabaseModule(), kwargs)
    module_obj.name = "setup-database"

    apply_module_options(
        module_obj,
        kwargs,
        {
            "PhaBox2_db": "phabox2_db",
            "phabox2_db_name": "phabox2_db_name",
            "phabox2_db_baselink": "phabox2_db_baselink",
            "genomad_db": "genomad_db",
            "checkv_db": "checkv_db",
            "eggNOG_db": "eggnog_db",
            "eggNOG_db_params": "eggnog_db_params",
            "virsorter2_db": "virsorter2_db",
            "iphop_db": "iphop_db",
            "humann_db": "humann_db",
            "hostile_index_db": "hostile_index_db",
            "metacerberus_db": "metacerberus_db",
            "metacerberus_setup_params": "metacerberus_setup_params",
            "pharokka_db": "pharokka_db",
            "checkm2_db": "checkm2_db",
            "GTDBTk_db": "gtdbtk_db",
            "GTDBTk_db_version": "gtdbtk_db_version",
            "iphop_db_version": "iphop_db_version",
            "iphop_db_basename": "iphop_db_basename",
            "vibrant_db": "vibrant_db",
        },
    )

    vomix_actions().run_module(
        "setup-database", module_obj, create_snakemake_flags(kwargs)
    )
    log.info(END_MODULE_RUN_LOG)
