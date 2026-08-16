import rich_click as click
import sys
import logging
import os
import platform
import time
from importlib.metadata import version
from rich.logging import RichHandler
from rich.console import Console

try:
    __version__ = version("vomix")
except Exception:
    __version__ = "unknown (development)"

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
    "--reset",
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
    "--list-conda-envs",
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
    log.info(f"vOMIX-MEGA v{__version__} initialized.")
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
    module_obj.reset = kwargs.get("reset")
    return module_obj


def create_snakemake_flags(kwargs):
    """Generates the SnakemakeFlags object automatically from kwargs."""
    from vomix.snakemakeFlags import SnakemakeFlags

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
        kwargs.get("list_conda_envs"),
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
        help="Working directory for the Snakemake workflow backend. Only modify for advanced customization or debugging purposes. (default: .)",
    )(function)
    function = click.option(
        "--outdir",
        default=None,
        required=False,
        help="Output directory for all pipeline results. New directories are automatically created; existing directories will be overwritten or appended. (default: results/)",
    )(function)
    function = click.option(
        "--datadir",
        default=None,
        required=False,
        help="Directory path for raw FASTQ files. Used to verify existing files or as the target destination for new downloads from NCBI SRA. (default: fastq/)",
    )(function)
    function = click.option(
        "--samplelist",
        default=None,
        required=False,
        help="Path to the sample_list.csv configuration file. This file defines input files and sample metadata. See documentation at https://github.com/holab-hku/vOMIX-MEGA/wiki for formatting specifications. (default: )",
    )(function)
    function = click.option(
        "--fasta",
        default=None,
        required=False,
        help="Path to a single FASTA input file for modules that accept single-file input. File extension must be .fasta, .fa, or .fna. (default: )",
    )(function)
    function = click.option(
        "--fastadir",
        default=None,
        required=False,
        help="Path to a directory containing multiple FASTA files. All files with .fasta, .fa, or .fna extensions will be automatically selected and processed. (default: )",
    )(function)
    function = click.option(
        "--sample-name",
        default=None,
        required=False,
        help="Sample name for output file naming when providing input via --fasta or config['fasta']. (default: )",
    )(function)
    function = click.option(
        "--assembly-ids",
        default=None,
        required=False,
        help='JSON-formatted array mapping sample names to input files (e.g., \'["sampleA", "SampleB"]\') when using --fasta-dir. This feature is currently under evaluation. (default: )',
    )(function)
    function = click.option(
        "--latest-run",
        default=None,
        required=False,
        help="Internal timestamp tracking the current execution run for history logging within the .vomix subdirectory. (default: )",
    )(function)
    function = click.option(
        "--keep-intermediates",
        is_flag=True,
        default=None,
        required=False,
        help="Retain intermediate processing files (e.g., fastp-cleaned raw FASTQ files before host decontamination). Useful for debugging or manual inspection. (default: False)",
    )(function)
    function = click.option(
        "--setup-database",
        is_flag=True,
        default=None,
        required=False,
        help="Initialize or update databases when executing modules other than 'setup-database'. Existing databases won't be reinstalled unless forced with Snakemake parameters like --forcerun or -F. (default: True)",
    )(function)
    function = click.option(
        "--max-cores",
        default=None,
        required=False,
        help="Maximum number of CPU cores allocated across parallel Snakemake tasks. Distinct from execution parameters like -j or -n. Currently under development. (default: 4)",
    )(function)
    function = click.option(
        "--NCBI-email",
        default=None,
        required=False,
        help="Email address for NCBI E-utilities API access. Required for data retrieval and download verification. (default: vomixtest@gmail.com)",
    )(function)
    function = click.option(
        "--NCBI-API-key",
        default=None,
        required=False,
        help="NCBI API key for higher throughput data retrieval. Obtain from NCBI following instructions at https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities. (default: )",
    )(function)
    function = click.option(
        "--custom-config",
        default=None,
        required=False,
        help="Path to your custom config.yml (default: None)",
    )(function)
    function = click.option(
        "--reset",
        is_flag=True,
        flag_value=True,
        default=False,
        required=False,
        help="Delete the log files to reset previously completed module run and its metadata, allowing module to be run again. (default: False)",
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
        help="Do not execute anything, display what would be done. (default: )",
    )(function)
    function = click.option(
        "--forceall",
        "-F",
        is_flag=True,
        flag_value=True,
        required=False,
        default=None,
        help="Force the execution of the selected rule and all dependencies. (default: )",
    )(function)
    function = click.option(
        "--configfile",
        default=None,
        required=False,
        help="Specify or overwrite the config file of the workflow. (default: )",
    )(function)
    function = click.option(
        "--unlock",
        is_flag=True,
        flag_value=True,
        required=False,
        default=None,
        help="Remove a lock on the working directory. (default: )",
    )(function)
    function = click.option(
        "--cores",
        "-c",
        default=4,
        required=False,
        help="Use at most N CPU cores/jobs in parallel. (default: 4)",
    )(function)
    function = click.option(
        "--jobs",
        "-j",
        default=4,
        required=False,
        help="Use at most N CPU cluster/cloud jobs in parallel. (default: 4)",
    )(function)
    function = click.option(
        "--latency-wait",
        default=20,
        required=False,
        help="Wait given seconds if an output file of a job is not present. (default: 20)",
    )(function)
    function = click.option(
        "--rerun-incomplete",
        "-ri",
        is_flag=True,
        flag_value=True,
        required=False,
        default=False,
        help="Re-run all jobs the output of which is recognized as incomplete. (default: )",
    )(function)
    function = click.option(
        "--rerun-triggers",
        multiple=True,
        required=False,
        default=None,
        help="Define what triggers the rerunning of a job. (default: )",
    )(function)
    function = click.option(
        "--software-deployment-method",
        "--deployment-method",
        "--deployment",
        "--sdm",
        multiple=True,
        required=False,
        default=["conda"],
        help="Specify software environment deployment method. (default: )",
    )(function)
    function = click.option(
        "--executor",
        "-e",
        required=False,
        default=None,
        help="Specify a custom executor, available via an executor plugin. (default: )",
    )(function)
    function = click.option(
        "--cluster-generic-submit-cmd",
        required=False,
        default=None,
        help="Command for submitting jobs. (default: )",
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
        "--list-conda-envs",
        is_flag=True,
        flag_value=True,
        required=False,
        default=None,
        help="List all conda environments and their location on disk. (default: False)",
    )(function)
    function = click.option(
        "--quiet",
        "-q",
        required=False,
        default=None,
        flag_value=True,
        help="Do not output certain information. (default: )",
    )(function)
    function = click.option(
        "--snakemake-args",
        required=False,
        default=None,
        help='Additional arguments to pass to the native snakemake command. Must be surround by double quotes `"`. (default: "")',
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
    help="Only download SRA data from NCBI without performing preprocessing steps. Downloads FASTQ files and stops. (default: False)",
)
@click.option(
    "--decontam-host/--no-decontam-host",
    is_flag=True,
    flag_value=True,
    default=None,
    required=False,
    help="Perform host decontamination after fastp quality trimming. Requires host indexes specified via hostile-index-name. (default: False)",
)
@click.option(
    "--dwnld-params",
    required=False,
    default=None,
    help="Additional parameters for NCBI SRA data retrieval via entrez-fetch. Customize download behavior and options. (default: )",
)
@click.option(
    "--pigz-params",
    required=False,
    default=None,
    help="Multi-threaded compression parameters for pigz when compressing downloaded FASTQ data. (default: )",
)
@click.option(
    "--fastp-params",
    required=False,
    default=None,
    help="Additional quality control parameters for fastp read trimming and filtering. (default: )",
)
@click.option(
    "--hostile-params",
    required=False,
    default=None,
    help="Additional parameters for Hostile host decontamination module. (default: )",
)
@click.option(
    "--hostile-aligner",
    required=False,
    default=None,
    help="Alignment algorithm for host decontamination. Options: 'bowtie2' (recommended for short reads) or 'minimap2'. (default: bowtie2)",
)
@click.option(
    "--hostile-aligner-params",
    required=False,
    default=None,
    help="Additional parameters for the selected Hostile alignment tool (bowtie2 or minimap2). (default: )",
)
@click.option(
    "--hostile-index-name",
    required=False,
    default=None,
    help="Name of pre-built Hostile indices. Available indices depend on installed Hostile version. See https://github.com/bede/hostile for available options. (default: human-t2t-hla)",
)
@click.option(
    "--hostile-index-db",
    required=False,
    default=None,
    help="Path to Hostile database directory for host decontamination. Database will be downloaded automatically if not present. (default: database/hostile)",
)
@snakemake_options
def run_preprocess(**kwargs):
    from vomix.snakemakeFlags import SnakemakeFlags
    from vomix.vomix_actions import vomix_actions
    from vomix.modules import PreProcessingModule

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
    help="Assembly engine for metagenome assembly. Currently supports MEGAHIT (recommended) and SPAdes (in development). (default: megahit)",
)
@click.option(
    "--megahit-min-len",
    required=False,
    default=None,
    help="Minimum contig length for MEGAHIT assembly filtering. Shorter contigs are excluded from output. (default: 300)",
)
@click.option(
    "--megahit-params",
    required=False,
    default=None,
    help="Additional parameters for MEGAHIT assembler. See MEGAHIT documentation for available options. (default: --prune-level 3)",
)
@click.option(
    "--spades-params",
    required=False,
    default=None,
    help="Additional parameters for SPAdes assembler. Currently supports metagenomic mode with --meta. (default: --meta)",
)
@click.option(
    "--spades-memory",
    required=False,
    default=None,
    help="Maximum RAM allocation in GB for SPAdes assembly. Adjust based on available system resources. (default: 250)",
)
@snakemake_options
def run_assembly(**kwargs):
    from vomix.snakemakeFlags import SnakemakeFlags
    from vomix.vomix_actions import vomix_actions
    from vomix.modules import AssemblyCoAssemblyModule

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
    help="Minimum contig length for viral identification analysis. Sequences shorter than this are excluded. (default: 0)",
)
@click.option(
    "--contig-splits",
    required=False,
    default=None,
    help="Number of parallel data partitions in viral identification tools to produce n+1 chunks of data. Used for reducing memory overhead (e.g. geNomad, VIBRANT, VirFinder, DeepVirFinder, PhaMer, VirSorter2). Set to 0 to disable partitioning. (default: 0)",
)
@click.option(
    "--genomad-db",
    required=False,
    default=None,
    help="Path to geNomad database directory for viral identification. Downloaded automatically if not present. (default: database/genomad)",
)
@click.option(
    "--genomad-min-len",
    required=False,
    default=None,
    help="Minimum contig length for geNomad analysis. Shorter sequences are excluded from viral classification. (default: 1000)",
)
@click.option(
    "--genomad-params",
    required=False,
    default=None,
    help="Additional parameters for geNomad viral identification. Optimize sensitivity/specificity using score calibration and relaxed settings. (default: --enable-score-calibration --relaxed)",
)
@click.option(
    "--genomad-cutoff",
    required=False,
    default=None,
    help="Minimum confidence threshold for geNomad viral classification. Higher values increase specificity at the cost of sensitivity. (default: 0.7)",
)
@click.option(
    "--genomad-cutoff-s",
    required=False,
    default=None,
    help="Secondary filtering confidence threshold for geNomad. Set to 0 to disable secondary filtering. (default: 0)",
)
@click.option(
    "--checkv-original/--checkv-pyhmmer",
    is_flag=True,
    flag_value=True,
    required=False,
    default=None,
    help="Use original CheckV implementation instead of the faster CheckV-PyHMMER version. CheckV-PyHMMER is recommended for large datasets. (default: False)",
)
@click.option(
    "--checkv-splits",
    required=False,
    default=None,
    help="Number of input data splits for CheckV or CheckV-PyHMMER analysis, splitting data into n+1 chunks. Higher values reduce memory usage but increase runtime. Set to 0 for no splitting. (default: 0)",
)
@click.option(
    "--checkv-params",
    required=False,
    default=None,
    help="Additional parameters for CheckV viral genome quality assessment. See CheckV documentation for available options. (default: )",
)
@click.option(
    "--checkv-database",
    required=False,
    default=None,
    help="Path to CheckV database directory for viral quality control and genome completeness assessment. (default: database/checkv)",
)
@click.option(
    "--clustering-fast/--clustering-cdhit",
    is_flag=True,
    flag_value=True,
    required=False,
    default=None,
    help="Use MEGABLAST-based fast clustering for viral operational taxonomic units (vOTUs). Disable to use CD-HIT (slower but sometimes more precise). (default: True)",
)
@click.option(
    "--cluster-iter",
    required=False,
    default=None,
    help="Number of clustering and pooling iterations (L) for vOTU clustering. Initially splits data into 2^(L-1) partitions for memory efficiency. Set to 1 to disable partitioning. (default: 1)",
)
@click.option(
    "--cdhit-params",
    required=False,
    default=None,
    help="Additional parameters for CD-HIT clustering. Used when clustering-fast is disabled. (default: -c 0.95 -aS 0.85 -d 400 -M 0 -n 5)",
)
@click.option(
    "--vOTU-ani",
    required=False,
    default=None,
    help="Average Nucleotide Identity (ANI) threshold for vOTU clustering (default 95% per MIUViG guidelines). (default: 95)",
)
@click.option(
    "--vOTU-targetcov",
    required=False,
    default=None,
    help="Minimum target coverage percentage for vOTU clustering (default 85% per MIUViG guidelines). (default: 85)",
)
@click.option(
    "--vOTU-querycov",
    required=False,
    default=None,
    help="Minimum query coverage percentage for vOTU clustering. Adjust for more permissive or stringent clustering. (default: 0)",
)
@snakemake_options
def run_viral_identify(**kwargs):
    from vomix.snakemakeFlags import SnakemakeFlags
    from vomix.vomix_actions import vomix_actions
    from vomix.modules import ViralIdentifyModule

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
    help="Path to PhaBox2 database directory for phage classification and analysis. Automatically downloaded if missing. (default: database/phabox_db_v2)",
)
@click.option(
    "--phabox2-db-name",
    required=False,
    default=None,
    help="Database name identifier for PhaBox2 tool execution. Automatically determined from the database version. (default: phabox_db_v2)",
)
@click.option(
    "--phabox2-db-baselink",
    required=False,
    default=None,
    help="Base URL for PhaBox2 database downloads. Used for updates and reinstallation. (default: https://github.com/KennthShang/PhaBOX/releases/download/v2)",
)
@click.option(
    "--genomad-db",
    required=False,
    default=None,
    help="Path to geNomad database directory for viral identification. Downloaded automatically if not present. (default: database/genomad)",
)
@click.option(
    "--virsorter2-db",
    required=False,
    default=None,
    help="Path to VirSorter2 database directory for viral sequence detection. Downloaded automatically if missing. (default: database/virsorter2)",
)
@click.option(
    "--vibrant-db",
    required=False,
    default=None,
    help="Path to VIBRANT database directory for viral functional annotation. Downloaded automatically if not present. (default: database/vibrant)",
)
@click.option(
    "--contig-min-len",
    required=False,
    default=None,
    help="Minimum contig length for viral identification analysis. Sequences shorter than this are excluded. (default: 0)",
)
@click.option(
    "--contig-splits",
    required=False,
    default=None,
    help="Number of parallel data partitions in viral identification tools to produce n+1 chunks of data. Used for reducing memory overhead (e.g. geNomad, VIBRANT, VirFinder, DeepVirFinder, PhaMer, VirSorter2). Set to 0 to disable partitioning. (default: 0)",
)
@click.option(
    "--genomad-min-len",
    required=False,
    default=None,
    help="Minimum contig length for geNomad analysis. Shorter sequences are excluded from viral classification. (default: 1000)",
)
@click.option(
    "--genomad-params",
    required=False,
    default=None,
    help="Additional parameters for geNomad viral identification. Optimize sensitivity/specificity using score calibration and relaxed settings. (default: --enable-score-calibration --relaxed)",
)
@click.option(
    "--genomad-cutoff",
    required=False,
    default=None,
    help="Minimum confidence threshold for geNomad viral classification. Higher values increase specificity at the cost of sensitivity. (default: 0.7)",
)
@click.option(
    "--genomad-cutoff-s",
    required=False,
    default=None,
    help="Secondary filtering confidence threshold for geNomad. Set to 0 to disable secondary filtering. (default: 0)",
)
@click.option(
    "--dvf-min-len",
    required=False,
    default=None,
    help="Minimum contig length for DeepVirFinder analysis. Shorter contigs are excluded from evaluation. (default: 1500)",
)
@click.option(
    "--phamer-min-len",
    required=False,
    default=None,
    help="Minimum contig length for PhaMer analysis. Shorter sequences are excluded from evaluation. (default: 2000)",
)
@click.option(
    "--dvf-params",
    required=False,
    default=None,
    help="Additional parameters for DeepVirFinder viral detection. Customize sensitivity and specificity. (default: )",
)
@click.option(
    "--phamer-params",
    required=False,
    default=None,
    help="Additional parameters for PhaMer execution. (default: )",
)
@click.option(
    "--virsorter2-params",
    required=False,
    default=None,
    help="Additional parameters for VirSorter2 viral detection. (default: )",
)
@click.option(
    "--vf-params",
    required=False,
    default=None,
    help="Additional parameters for VirFinder viral detection tool. (default: )",
)
@click.option(
    "--seeker-params",
    required=False,
    default=None,
    help="Additional parameters for Seeker deep-learning viral detection. (default: )",
)
@click.option(
    "--PPR-params",
    required=False,
    default=None,
    help="Additional parameters for PPR-META viral detection. (default: )",
)
@click.option(
    "--dvf-cutoff",
    required=False,
    default=None,
    help="Minimum confidence score for DeepVirFinder viral classification. Higher values increase specificity. (default: 0.7)",
)
@click.option(
    "--dvf-pval",
    required=False,
    default=None,
    help="Maximum p-value threshold for DeepVirFinder viral classification. More significant p-values indicate higher confidence. (default: 0.05)",
)
@click.option(
    "--phamer-pred",
    required=False,
    default=None,
    help="PhaMer prediction category. Options: 'phage' (detect phages) or 'non-phage'. (default: phage)",
)
@click.option(
    "--phamer-cutoff",
    required=False,
    default=None,
    help="Minimum confidence threshold for PhaMer viral classification. Adjust based on desired sensitivity/specificity. (default: 0)",
)
@click.option(
    "--vf-cutoff",
    required=False,
    default=None,
    help="Minimum confidence threshold for VirFinder viral classification. (default: 0)",
)
@click.option(
    "--virsorter2-cutoff",
    required=False,
    default=None,
    help="Minimum confidence score for VirSorter2 viral classification. (default: 0)",
)
@click.option(
    "--seeker-cutoff",
    required=False,
    default=None,
    help="Minimum confidence threshold for Seeker viral classification. (default: 0)",
)
@click.option(
    "--ppr-cutoff",
    required=False,
    default=None,
    help="Minimum confidence threshold for PPR-META viral classification. (default: 0)",
)
@click.option(
    "--vibrant-cutoff",
    required=False,
    default=None,
    help="Minimum confidence threshold for VIBRANT viral classification. (default: 0)",
)
@click.option(
    "--checkv-original/--checkv-pyhmmer",
    is_flag=True,
    flag_value=True,
    required=False,
    default=None,
    help="Use original CheckV implementation instead of the faster CheckV-PyHMMER version. CheckV-PyHMMER is recommended for large datasets. (default: False)",
)
@click.option(
    "--checkv-splits",
    required=False,
    default=None,
    help="Number of input data splits for CheckV or CheckV-PyHMMER analysis, splitting data into n+1 chunks. Higher values reduce memory usage but increase runtime. Set to 0 for no splitting. (default: 0)",
)
@click.option(
    "--checkv-params",
    required=False,
    default=None,
    help="Additional parameters for CheckV viral genome quality assessment. See CheckV documentation for available options. (default: )",
)
@click.option(
    "--checkv-database",
    required=False,
    default=None,
    help="Path to CheckV database directory for viral quality control and genome completeness assessment. (default: database/checkv)",
)
@click.option(
    "--clustering-fast/--clustering-cdhit",
    is_flag=True,
    flag_value=True,
    required=False,
    default=None,
    help="Use MEGABLAST-based fast clustering for viral operational taxonomic units (vOTUs). Disable to use CD-HIT (slower but sometimes more precise). (default: True)",
)
@click.option(
    "--cluster-iter",
    required=False,
    default=None,
    help="Number of clustering and pooling iterations (L) for vOTU clustering. Initially splits data into 2^(L-1) partitions for memory efficiency. Set to 1 to disable partitioning. (default: 1)",
)
@click.option(
    "--cdhit-params",
    required=False,
    default=None,
    help="Additional parameters for CD-HIT clustering. Used when clustering-fast is disabled. (default: -c 0.95 -aS 0.85 -d 400 -M 0 -n 5)",
)
@click.option(
    "--vOTU-ani",
    required=False,
    default=None,
    help="Average Nucleotide Identity (ANI) threshold for vOTU clustering (default 95% per MIUViG guidelines). (default: 95)",
)
@click.option(
    "--vOTU-targetcov",
    required=False,
    default=None,
    help="Minimum target coverage percentage for vOTU clustering (default 85% per MIUViG guidelines). (default: 85)",
)
@click.option(
    "--vOTU-querycov",
    required=False,
    default=None,
    help="Minimum query coverage percentage for vOTU clustering. Adjust for more permissive or stringent clustering. (default: 0)",
)
@snakemake_options
def run_viral_benchmark(**kwargs):
    from vomix.snakemakeFlags import SnakemakeFlags
    from vomix.vomix_actions import vomix_actions
    from vomix.modules import ViralBenchmarkModule

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
    help="Path to geNomad database directory for viral identification. Downloaded automatically if not present. (default: database/genomad)",
)
@click.option(
    "--PhaBox2-db",
    required=False,
    default=None,
    help="Path to PhaBox2 database directory for phage classification and analysis. Automatically downloaded if missing. (default: database/phabox_db_v2)",
)
@click.option(
    "--phabox2-db-name",
    required=False,
    default=None,
    help="Database name identifier for PhaBox2 tool execution. Automatically determined from the database version. (default: phabox_db_v2)",
)
@click.option(
    "--phabox2-db-baselink",
    required=False,
    default=None,
    help="Base URL for PhaBox2 database downloads. Used for updates and reinstallation. (default: https://github.com/KennthShang/PhaBOX/releases/download/v2)",
)
@click.option(
    "--phagcn-min-len",
    required=False,
    default=None,
    help="Minimum contig length for PhaGCN taxonomic classification. Shorter sequences are excluded. (default: 1000)",
)
@click.option(
    "--phagcn-params",
    required=False,
    default=None,
    help="Additional parameters for PhaGCN taxonomy classification. (default: )",
)
@click.option(
    "--genomad-params-tax",
    required=False,
    default=None,
    help="Additional parameters for geNomad taxonomy classification. Customize taxonomic assignment behavior. (default: --enable-score-calibration --relaxed)",
)
@snakemake_options
def run_viral_taxonomy(**kwargs):
    from vomix.snakemakeFlags import SnakemakeFlags
    from vomix.vomix_actions import vomix_actions
    from vomix.modules import ViralTaxonomyModule

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
    help="Path to PhaBox2 database directory for phage classification and analysis. Automatically downloaded if missing. (default: database/phabox_db_v2)",
)
@click.option(
    "--phabox2-db-name",
    required=False,
    default=None,
    help="Database name identifier for PhaBox2 tool execution. Automatically determined from the database version. (default: phabox_db_v2)",
)
@click.option(
    "--phabox2-db-baselink",
    required=False,
    default=None,
    help="Base URL for PhaBox2 database downloads. Used for updates and reinstallation. (default: https://github.com/KennthShang/PhaBOX/releases/download/v2)",
)
@click.option(
    "--CHERRY-params",
    required=False,
    default=None,
    help="Additional parameters for CHERRY host prediction algorithm. (default: )",
)
@click.option(
    "--PhaTYP-params",
    required=False,
    default=None,
    help="Additional parameters for PhaTYP phage lifestyle prediction. (default: )",
)
@click.option(
    "--iphop-host/--cherry-host",
    is_flag=True,
    flag_value=True,
    required=False,
    default=None,
    help="Use iPHoP for viral host prediction instead of CHERRY. iPHoP requires more memory (>100GB) and larger database. (default: False)",
)
@click.option(
    "--iphop-cutoff",
    required=False,
    default=None,
    help="Minimum confidence threshold for iPHoP host classification. Higher values increase accuracy at the cost of assignment coverage. (default: 90)",
)
@click.option(
    "--iphop-params",
    required=False,
    default=None,
    help="Additional parameters for iPHoP host prediction. (default: )",
)
@click.option(
    "--iphop-db",
    required=False,
    default=None,
    help="Path to iPHoP database directory for viral host prediction. Note this database is large (>100GB) and requires substantial disk space. (default: database/iphop)",
)
@click.option(
    "--iphop-db-version",
    required=False,
    default=None,
    help="Version identifier for iPHoP database. Verify compatibility at https://bitbucket.org/srouxjgi/iphop/src/main/#markdown-header-host-databases-and-versions. (default: iPHoP_db_Aug23_rw)",
)
@click.option(
    "--iphop-db-basename",
    required=False,
    default=None,
    help="Base name of the iPHoP database. Must match the version specified by iphop-db-version. See documentation for compatibility. (default: Aug_2023_pub_rw)",
)
@snakemake_options
def run_viral_host(**kwargs):
    from vomix.snakemakeFlags import SnakemakeFlags
    from vomix.vomix_actions import vomix_actions
    from vomix.modules import ViralHostModule

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
    help="Additional parameters for CoverM read mapping and coverage calculation. Customize mapping sensitivity and coverage metrics. (default: --mapper minimap2-sr --min-read-percent-identity 95 --min-read-aligned-percent 75 --trim-min 10 --trim-max 90)",
)
@click.option(
    "--coverm-methods",
    required=False,
    default=None,
    help="CoverM coverage metrics to calculate. Available options include tpm, rpkm, and relative abundance. (default: tpm rpkm)",
)
@snakemake_options
def run_viral_community(**kwargs):
    from vomix.snakemakeFlags import SnakemakeFlags
    from vomix.vomix_actions import vomix_actions
    from vomix.modules import ViralCommunityModule

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
    help="Path to PhaBox2 database directory for phage classification and analysis. Automatically downloaded if missing. (default: database/phabox_db_v2)",
)
@click.option(
    "--phabox2-db-name",
    required=False,
    default=None,
    help="Database name identifier for PhaBox2 tool execution. Automatically determined from the database version. (default: phabox_db_v2)",
)
@click.option(
    "--phabox2-db-baselink",
    required=False,
    default=None,
    help="Base URL for PhaBox2 database downloads. Used for updates and reinstallation. (default: https://github.com/KennthShang/PhaBOX/releases/download/v2)",
)
@click.option(
    "--eggNOG-params",
    required=False,
    default=None,
    help="Additional parameters for eggNOG-mapper viral functional annotation. Customize database search and annotation. (default: )",
)
@click.option(
    "--PhaVIP-params",
    required=False,
    default=None,
    help="Additional parameters for PhaVIP viral pathogenicity prediction. (default: )",
)
@click.option(
    "--metacerberus-db",
    required=False,
    default=None,
    help="Path to MetaCerberus database directory for comprehensive viral annotation. (default: database/metacerberus)",
)
@click.option(
    "--metacerberus-setup-params",
    required=False,
    default=None,
    help="Additional parameters for MetaCerberus database setup and indexing. (default: )",
)
@click.option(
    "--metacerberus-params",
    required=False,
    default=None,
    help="Additional parameters for MetaCerberus viral annotation. Customize HMM database and annotation options. (default: --hmm ALL)",
)
@click.option(
    "--pharokka-db",
    required=False,
    default=None,
    help="Path to pharokka database directory for bacteriophage gene annotation. (default: database/pharokka)",
)
@click.option(
    "--pharokka-params",
    required=False,
    default=None,
    help="Additional parameters for pharokka phage annotation. Customize gene prediction and functional annotation. (default: -g prodigal-gv --meta)",
)
@snakemake_options
def run_viral_annotate(**kwargs):
    from vomix.snakemakeFlags import SnakemakeFlags
    from vomix.vomix_actions import vomix_actions
    from vomix.modules import ViralAnnotateModule

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
    help="Additional parameters for MetaPhlAn prokaryotic community profiling. (default: --ignore_eukaryotes)",
)
@click.option(
    "--mpa-indexv",
    required=False,
    default=None,
    help="MetaPhlAn database version. The database will be automatically downloaded if not present. (default: mpa_vOct22_CHOCOPhlAnSGB_202212)",
)
@snakemake_options
def run_prok_community(**kwargs):
    from vomix.snakemakeFlags import SnakemakeFlags
    from vomix.vomix_actions import vomix_actions
    from vomix.modules import ProkaryoticCommunityModule

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
    help="Path to CheckM2 database directory for genome quality assessment and completeness evaluation. (default: database/checkm2)",
)
@click.option(
    "--GTDBTk-db",
    required=False,
    default=None,
    help="Path to GTDB-Tk database directory for taxonomic classification of prokaryotic genomes. (default: database/GTDB-Tk)",
)
@click.option(
    "--GTDBTk-db-version",
    required=False,
    default=None,
    type=str,
    help="GTDB-Tk database version. Must match the installed GTDB-Tk version. See https://ecogenomics.github.io/GTDBTk/installing/index.html#gtdb-tk-reference-data. (default: 232)",
)
@click.option(
    "--GTDBTk-identify-params",
    required=False,
    default=None,
    help="Additional parameters for GTDB-Tk identify step. Customize marker gene identification and selection. (default: )",
)
@click.option(
    "--GTDBTk-align-params",
    required=False,
    default=None,
    help="Additional parameters for GTDB-Tk align step. Customize multiple sequence alignment. (default: )",
)
@click.option(
    "--GTDBTk-classify-params",
    required=False,
    default=None,
    help="Additional parameters for GTDB-Tk classify step. Customize taxonomic classification. (default: )",
)
@click.option(
    "--VAMB-params",
    required=False,
    default=None,
    help="Additional parameters for VAMB binning. Customize neural network clustering parameters. (default: -m 100)",
)
@click.option(
    "--binning-consensus/--binning-gpu",
    is_flag=True,
    flag_value=True,
    required=False,
    default=None,
    help="Use consensus binning with MetaBAT2, MaxBin2, and CONCOCT via DASTool. Disable to use GPU-accelerated VAMB binning. (default: True)",
)
@click.option(
    "--strobealign-params",
    required=False,
    default=None,
    help="Additional parameters for strobealign read alignment during binning. (default: )",
)
@click.option(
    "--MetaBAT2-params",
    required=False,
    default=None,
    help="Additional parameters for MetaBAT2 metagenomic binning. (default: -m 1500)",
)
@click.option(
    "--MaxBin2-params",
    required=False,
    default=None,
    help="Additional parameters for MaxBin2 metagenomic binning. Customize iterative binning and probability thresholds. (default: -min_contig_length 1500 -max_iteration 50 -prob_threshold 0.9)",
)
@click.option(
    "--CONCOCT-params",
    required=False,
    default=None,
    help="Additional parameters for CONCOCT metagenomic binning. (default: )",
)
@click.option(
    "--jgi-summarize-params",
    required=False,
    default=None,
    help="Additional parameters for jgi_summarize_bam_contig_depth. Customize read depth calculation and identity thresholds. (default: --percentIdentity 97)",
)
@click.option(
    "--DASTool-params",
    required=False,
    default=None,
    help="Additional parameters for DASTool bin consolidation. Customize scoring and bin selection. (default: )",
)
@click.option(
    "--checkm2-params",
    required=False,
    default=None,
    help="Additional parameters for CheckM2 genome quality assessment. (default: )",
)
@click.option(
    "--galah-params",
    required=False,
    default=None,
    help="Additional parameters for galah bin dereplication. Customize ANI threshold and alignment parameters. (default: --ani 95 --min-aligned-fraction 15 --fragment-length 3000)",
)
@snakemake_options
def run_prok_binning(**kwargs):
    from vomix.snakemakeFlags import SnakemakeFlags
    from vomix.vomix_actions import vomix_actions
    from vomix.modules import ProkaryoticBinningModule

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
    help="Additional parameters for HUMAnN3 functional annotation. Customize mapping and statistical analysis. (default: --remove-temp-output)",
)
@snakemake_options
def run_prok_annotate(**kwargs):
    from vomix.snakemakeFlags import SnakemakeFlags
    from vomix.vomix_actions import vomix_actions
    from vomix.modules import ProkaryoticAnnotateModule

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
    help="Only download SRA data from NCBI without performing preprocessing steps. Downloads FASTQ files and stops. (default: False)",
)
@click.option(
    "--decontam-host/--no-decontam-host",
    is_flag=True,
    flag_value=True,
    default=None,
    required=False,
    help="Perform host decontamination after fastp quality trimming. Requires host indexes specified via hostile-index-name. (default: False)",
)
@click.option(
    "--dwnld-params",
    required=False,
    default=None,
    help="Additional parameters for NCBI SRA data retrieval via entrez-fetch. Customize download behavior and options. (default: )",
)
@click.option(
    "--pigz-params",
    required=False,
    default=None,
    help="Multi-threaded compression parameters for pigz when compressing downloaded FASTQ data. (default: )",
)
@click.option(
    "--fastp-params",
    required=False,
    default=None,
    help="Additional quality control parameters for fastp read trimming and filtering. (default: )",
)
@click.option(
    "--hostile-params",
    required=False,
    default=None,
    help="Additional parameters for Hostile host decontamination module. (default: )",
)
@click.option(
    "--hostile-aligner",
    required=False,
    default=None,
    help="Alignment algorithm for host decontamination. Options: 'bowtie2' (recommended for short reads) or 'minimap2'. (default: bowtie2)",
)
@click.option(
    "--hostile-aligner-params",
    required=False,
    default=None,
    help="Additional parameters for the selected Hostile alignment tool (bowtie2 or minimap2). (default: )",
)
@click.option(
    "--hostile-index-name",
    required=False,
    default=None,
    help="Name of pre-built Hostile indices. Available indices depend on installed Hostile version. See https://github.com/bede/hostile for available options. (default: human-t2t-hla)",
)
@click.option(
    "--hostile-index-db",
    required=False,
    default=None,
    help="Path to Hostile database directory for host decontamination. Database will be downloaded automatically if not present. (default: database/hostile)",
)
@click.option(
    "--assembler",
    default=None,
    required=False,
    help="Assembly engine for metagenome assembly. Currently supports MEGAHIT (recommended) and SPAdes (in development). (default: megahit)",
)
@click.option(
    "--megahit-min-len",
    required=False,
    default=None,
    help="Minimum contig length for MEGAHIT assembly filtering. Shorter contigs are excluded from output. (default: 300)",
)
@click.option(
    "--megahit-params",
    required=False,
    default=None,
    help="Additional parameters for MEGAHIT assembler. See MEGAHIT documentation for available options. (default: --prune-level 3)",
)
@click.option(
    "--spades-params",
    required=False,
    default=None,
    help="Additional parameters for SPAdes assembler. Currently supports metagenomic mode with --meta. (default: --meta)",
)
@click.option(
    "--spades-memory",
    required=False,
    default=None,
    help="Maximum RAM allocation in GB for SPAdes assembly. Adjust based on available system resources. (default: 250)",
)
@click.option(
    "--contig-min-len",
    required=False,
    default=None,
    help="Minimum contig length for viral identification analysis. Sequences shorter than this are excluded. (default: 0)",
)
@click.option(
    "--contig-splits",
    required=False,
    default=None,
    help="Number of parallel data partitions in viral identification tools to produce n+1 chunks of data. Used for reducing memory overhead (e.g. geNomad, VIBRANT, VirFinder, DeepVirFinder, PhaMer, VirSorter2). Set to 0 to disable partitioning. (default: 0)",
)
@click.option(
    "--genomad-db",
    required=False,
    default=None,
    help="Path to geNomad database directory for viral identification. Downloaded automatically if not present. (default: database/genomad)",
)
@click.option(
    "--genomad-min-len",
    required=False,
    default=None,
    help="Minimum contig length for geNomad analysis. Shorter sequences are excluded from viral classification. (default: 1000)",
)
@click.option(
    "--genomad-params",
    required=False,
    default=None,
    help="Additional parameters for geNomad viral identification. Optimize sensitivity/specificity using score calibration and relaxed settings. (default: --enable-score-calibration --relaxed)",
)
@click.option(
    "--genomad-cutoff",
    required=False,
    default=None,
    help="Minimum confidence threshold for geNomad viral classification. Higher values increase specificity at the cost of sensitivity. (default: 0.7)",
)
@click.option(
    "--genomad-cutoff-s",
    required=False,
    default=None,
    help="Secondary filtering confidence threshold for geNomad. Set to 0 to disable secondary filtering. (default: 0)",
)
@click.option(
    "--checkv-original/--checkv-pyhmmer",
    is_flag=True,
    flag_value=True,
    required=False,
    default=None,
    help="Use original CheckV implementation instead of the faster CheckV-PyHMMER version. CheckV-PyHMMER is recommended for large datasets. (default: False)",
)
@click.option(
    "--checkv-splits",
    required=False,
    default=None,
    help="Number of input data splits for CheckV or CheckV-PyHMMER analysis, splitting data into n+1 chunks. Higher values reduce memory usage but increase runtime. Set to 0 for no splitting. (default: 0)",
)
@click.option(
    "--checkv-params",
    required=False,
    default=None,
    help="Additional parameters for CheckV viral genome quality assessment. See CheckV documentation for available options. (default: )",
)
@click.option(
    "--checkv-database",
    required=False,
    default=None,
    help="Path to CheckV database directory for viral quality control and genome completeness assessment. (default: database/checkv)",
)
@click.option(
    "--clustering-fast/--clustering-cdhit",
    is_flag=True,
    flag_value=True,
    required=False,
    default=None,
    help="Use MEGABLAST-based fast clustering for viral operational taxonomic units (vOTUs). Disable to use CD-HIT (slower but sometimes more precise). (default: True)",
)
@click.option(
    "--cluster-iter",
    required=False,
    default=None,
    help="Number of clustering and pooling iterations (L) for vOTU clustering. Initially splits data into 2^(L-1) partitions for memory efficiency. Set to 1 to disable partitioning. (default: 1)",
)
@click.option(
    "--cdhit-params",
    required=False,
    default=None,
    help="Additional parameters for CD-HIT clustering. Used when clustering-fast is disabled. (default: -c 0.95 -aS 0.85 -d 400 -M 0 -n 5)",
)
@click.option(
    "--vOTU-ani",
    required=False,
    default=None,
    help="Average Nucleotide Identity (ANI) threshold for vOTU clustering (default 95% per MIUViG guidelines). (default: 95)",
)
@click.option(
    "--vOTU-targetcov",
    required=False,
    default=None,
    help="Minimum target coverage percentage for vOTU clustering (default 85% per MIUViG guidelines). (default: 85)",
)
@click.option(
    "--vOTU-querycov",
    required=False,
    default=None,
    help="Minimum query coverage percentage for vOTU clustering. Adjust for more permissive or stringent clustering. (default: 0)",
)
@click.option(
    "--PhaBox2-db",
    required=False,
    default=None,
    help="Path to PhaBox2 database directory for phage classification and analysis. Automatically downloaded if missing. (default: database/phabox_db_v2)",
)
@click.option(
    "--phabox2-db-name",
    required=False,
    default=None,
    help="Database name identifier for PhaBox2 tool execution. Automatically determined from the database version. (default: phabox_db_v2)",
)
@click.option(
    "--phabox2-db-baselink",
    required=False,
    default=None,
    help="Base URL for PhaBox2 database downloads. Used for updates and reinstallation. (default: https://github.com/KennthShang/PhaBOX/releases/download/v2)",
)
@click.option(
    "--phagcn-min-len",
    required=False,
    default=None,
    help="Minimum contig length for PhaGCN taxonomic classification. Shorter sequences are excluded. (default: 1000)",
)
@click.option(
    "--phagcn-params",
    required=False,
    default=None,
    help="Additional parameters for PhaGCN taxonomy classification. (default: )",
)
@click.option(
    "--genomad-params-tax",
    required=False,
    default=None,
    help="Additional parameters for geNomad taxonomy classification. Customize taxonomic assignment behavior. (default: --enable-score-calibration --relaxed)",
)
@click.option(
    "--CHERRY-params",
    required=False,
    default=None,
    help="Additional parameters for CHERRY host prediction algorithm. (default: )",
)
@click.option(
    "--PhaTYP-params",
    required=False,
    default=None,
    help="Additional parameters for PhaTYP phage lifestyle prediction. (default: )",
)
@click.option(
    "--iphop-host/--cherry-host",
    is_flag=True,
    flag_value=True,
    required=False,
    default=None,
    help="Use iPHoP for viral host prediction instead of CHERRY. iPHoP requires more memory (>100GB) and larger database. (default: False)",
)
@click.option(
    "--iphop-cutoff",
    required=False,
    default=None,
    help="Minimum confidence threshold for iPHoP host classification. Higher values increase accuracy at the cost of assignment coverage. (default: 90)",
)
@click.option(
    "--iphop-params",
    required=False,
    default=None,
    help="Additional parameters for iPHoP host prediction. (default: )",
)
@click.option(
    "--iphop-db",
    required=False,
    default=None,
    help="Path to iPHoP database directory for viral host prediction. Note this database is large (>100GB) and requires substantial disk space. (default: database/iphop)",
)
@click.option(
    "--eggNOG-params",
    required=False,
    default=None,
    help="Additional parameters for eggNOG-mapper viral functional annotation. Customize database search and annotation. (default: )",
)
@click.option(
    "--PhaVIP-params",
    required=False,
    default=None,
    help="Additional parameters for PhaVIP viral pathogenicity prediction. (default: )",
)
@click.option(
    "--metacerberus-db",
    required=False,
    default=None,
    help="Path to MetaCerberus database directory for comprehensive viral annotation. (default: database/metacerberus)",
)
@click.option(
    "--metacerberus-setup-params",
    required=False,
    default=None,
    help="Additional parameters for MetaCerberus database setup and indexing. (default: )",
)
@click.option(
    "--metacerberus-params",
    required=False,
    default=None,
    help="Additional parameters for MetaCerberus viral annotation. Customize HMM database and annotation options. (default: --hmm ALL)",
)
@click.option(
    "--pharokka-db",
    required=False,
    default=None,
    help="Path to pharokka database directory for bacteriophage gene annotation. (default: database/pharokka)",
)
@click.option(
    "--pharokka-params",
    required=False,
    default=None,
    help="Additional parameters for pharokka phage annotation. Customize gene prediction and functional annotation. (default: -g prodigal-gv --meta)",
)
@click.option(
    "--coverm-params",
    required=False,
    default=None,
    help="Additional parameters for CoverM read mapping and coverage calculation. Customize mapping sensitivity and coverage metrics. (default: --mapper minimap2-sr --min-read-percent-identity 95 --min-read-aligned-percent 75 --trim-min 10 --trim-max 90)",
)
@click.option(
    "--coverm-methods",
    required=False,
    default=None,
    help="CoverM coverage metrics to calculate. Available options include tpm, rpkm, and relative abundance. (default: tpm rpkm)",
)
@snakemake_options
def run_viral_end_to_end(**kwargs):
    from vomix.snakemakeFlags import SnakemakeFlags
    from vomix.vomix_actions import vomix_actions
    from vomix.modules import ViralEndToEndModule

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
    help="Use MEGABLAST-based fast clustering for viral operational taxonomic units (vOTUs). Disable to use CD-HIT (slower but sometimes more precise). (default: True)",
)
@click.option(
    "--cluster-iter",
    required=False,
    default=None,
    help="Number of clustering and pooling iterations (L) for vOTU clustering. Initially splits data into 2^(L-1) partitions for memory efficiency. Set to 1 to disable partitioning. (default: 1)",
)
@click.option(
    "--cdhit-params",
    required=False,
    default=None,
    help="Additional parameters for CD-HIT clustering. Used when clustering-fast is disabled. (default: -c 0.95 -aS 0.85 -d 400 -M 0 -n 5)",
)
@click.option(
    "--vOTU-ani",
    required=False,
    default=None,
    help="Average Nucleotide Identity (ANI) threshold for vOTU clustering (default 95% per MIUViG guidelines). (default: 95)",
)
@click.option(
    "--vOTU-targetcov",
    required=False,
    default=None,
    help="Minimum target coverage percentage for vOTU clustering (default 85% per MIUViG guidelines). (default: 85)",
)
@click.option(
    "--vOTU-querycov",
    required=False,
    default=None,
    help="Minimum query coverage percentage for vOTU clustering. Adjust for more permissive or stringent clustering. (default: 0)",
)
@snakemake_options
def run_cluster_fast(**kwargs):
    from vomix.snakemakeFlags import SnakemakeFlags
    from vomix.vomix_actions import vomix_actions
    from vomix.modules import ClusterFastModule

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
    help="Use original CheckV implementation instead of the faster CheckV-PyHMMER version. CheckV-PyHMMER is recommended for large datasets. (default: False)",
)
@click.option(
    "--checkv-splits",
    required=False,
    default=None,
    help="Number of input data splits for CheckV or CheckV-PyHMMER analysis, splitting data into n+1 chunks. Higher values reduce memory usage but increase runtime. Set to 0 for no splitting. (default: 0)",
)
@click.option(
    "--checkv-params",
    required=False,
    default=None,
    help="Additional parameters for CheckV viral genome quality assessment. See CheckV documentation for available options. (default: )",
)
@click.option(
    "--checkv-database",
    required=False,
    default=None,
    help="Path to CheckV database directory for viral quality control and genome completeness assessment. (default: database/checkv)",
)
@snakemake_options
def run_checkv_pyhmmer(**kwargs):
    from vomix.snakemakeFlags import SnakemakeFlags
    from vomix.vomix_actions import vomix_actions
    from vomix.modules import CheckVPyHMMERModule

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
    help="Path to Hostile database directory for host decontamination. Database will be downloaded automatically if not present. (default: database/hostile)",
)
@click.option(
    "--PhaBox2-db",
    required=False,
    default=None,
    help="Path to PhaBox2 database directory for phage classification and analysis. Automatically downloaded if missing. (default: database/phabox_db_v2)",
)
@click.option(
    "--phabox2-db-name",
    required=False,
    default=None,
    help="Database name identifier for PhaBox2 tool execution. Automatically determined from the database version. (default: phabox_db_v2)",
)
@click.option(
    "--phabox2-db-baselink",
    required=False,
    default=None,
    help="Base URL for PhaBox2 database downloads. Used for updates and reinstallation. (default: https://github.com/KennthShang/PhaBOX/releases/download/v2)",
)
@click.option(
    "--genomad-db",
    required=False,
    default=None,
    help="Path to geNomad database directory for viral identification. Downloaded automatically if not present. (default: database/genomad)",
)
@click.option(
    "--virsorter2-db",
    required=False,
    default=None,
    help="Path to VirSorter2 database directory for viral sequence detection. Downloaded automatically if missing. (default: database/virsorter2)",
)
@click.option(
    "--vibrant-db",
    required=False,
    default=None,
    help="Path to VIBRANT database directory for viral functional annotation. Downloaded automatically if not present. (default: database/vibrant)",
)
@click.option(
    "--checkv-db",
    required=False,
    default=None,
    help="Path to CheckV database directory for viral quality control and genome completeness assessment. (default: database/checkv)",
)
@click.option(
    "--eggNOG-db",
    required=False,
    default=None,
    help="Path to eggNOG database directory for functional annotation. Supports eggNOG-mapper v2. (default: database/eggNOGv2)",
)
@click.option(
    "--eggNOG-db-params",
    required=False,
    default=None,
    help="Additional parameters for eggNOG database download and setup. See https://github.com/eggnogdb/eggnog-mapper for options. (default: )",
)
@click.option(
    "--checkm2-db",
    required=False,
    default=None,
    help="Path to CheckM2 database directory for genome quality assessment and completeness evaluation. (default: database/checkm2)",
)
@click.option(
    "--GTDBTk-db",
    required=False,
    default=None,
    help="Path to GTDB-Tk database directory for taxonomic classification of prokaryotic genomes. (default: database/GTDB-Tk)",
)
@click.option(
    "--GTDBTk-db-version",
    required=False,
    default=None,
    type=str,
    help="GTDB-Tk database version. Must match the installed GTDB-Tk version. See https://ecogenomics.github.io/GTDBTk/installing/index.html#gtdb-tk-reference-data. (default: 232)",
)
@click.option(
    "--iphop-db",
    required=False,
    default=None,
    help="Path to iPHoP database directory for viral host prediction. Note this database is large (>100GB) and requires substantial disk space. (default: database/iphop)",
)
@click.option(
    "--iphop-db-version",
    required=False,
    default=None,
    help="Version identifier for iPHoP database. Verify compatibility at https://bitbucket.org/srouxjgi/iphop/src/main/#markdown-header-host-databases-and-versions. (default: iPHoP_db_Aug23_rw)",
)
@click.option(
    "--iphop-db-basename",
    required=False,
    default=None,
    help="Base name of the iPHoP database. Must match the version specified by iphop-db-version. See documentation for compatibility. (default: Aug_2023_pub_rw)",
)
@click.option(
    "--humann-db",
    required=False,
    default=None,
    help="Path to HUMAnN3 database directory for functional profiling of prokaryotic communities. (default: database/humann)",
)
@click.option(
    "--metacerberus-db",
    required=False,
    default=None,
    help="Path to MetaCerberus database directory for comprehensive viral annotation. (default: database/metacerberus)",
)
@click.option(
    "--metacerberus-setup-params",
    required=False,
    default=None,
    help="Additional parameters for MetaCerberus database setup and indexing. (default: )",
)
@click.option(
    "--pharokka-db",
    required=False,
    default=None,
    help="Path to pharokka database directory for bacteriophage gene annotation. (default: database/pharokka)",
)
@snakemake_options
def run_setup_database(**kwargs):
    from vomix.snakemakeFlags import SnakemakeFlags
    from vomix.vomix_actions import vomix_actions
    from vomix.modules import SetupDatabaseModule

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
