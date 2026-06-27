import rich_click as click
import sys
import logging
import os
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

logging.basicConfig(level=logging.INFO)

modules_list = [
    "assembly",
    "checkv-pyhmmer",
    "checkv",
    "clustering-fast",
    "clustering-sensitive",
    "host-cherry",
    "host",
    "preprocess-decontam",
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
    "viral-benchmark",
]

END_MODULE_RUN_LOG = "End module run"


def use_last_options_check(ctx, value):
    i = 2
    total_params = len(ctx.command.params)
    if value == True:
        while i < total_params:
            ctx.command.params[i].prompt_required = False
            ctx.command.params[i].required = False
            i += 1


# common options decorator
def common_options(function):
    function = click.option(
        "--workdir",
        default=None,
        required=False,
        help="The working directory for the underlying Snakemake workflow back-end. Modify only for advanced customization or debugging purposes.",
    )(function)
    function = click.option(
        "--outdir",
        default="result/",
        required=False,
        help="The directory path where structured output results from vOMIX-MEGA will be deposited. New directories will be automatically created (including nested structures). Existing directories will be overwritten or appended",
    )(function)
    function = click.option(
        "--datadir",
        default=None,
        required=False,
        help="The path to raw FASTQ files. Used to verify if files are pre-downloaded or as the target destination for new downloads.",
    )(function)
    function = click.option(
        "--samplelist",
        default=None,
        required=False,
        help="The path to the sample_list.csv configuration file for inputs and file paths. For detailed specifications visit https://vomix-mega.readthedocs.io/en/latest/",
    )(function)
    function = click.option(
        "--fasta",
        default=None,
        required=False,
        help="The path to a single input FASTA file for modules that accept a single file as valid input. File path must end with .fasta, .fa, or .fna",
    )(function)
    function = click.option(
        "--fastadir",
        default=None,
        required=False,
        help="The path to a directory containing input FASTA files. The workflow automatically selects all files within this directory ending with .fasta, .fa, or .fna",
    )(function)
    function = click.option(
        "--sample-name",
        default=None,
        required=False,
        help='The explicit sample name utilized for output file naming when providing inputs via --fasta or config["fasta"].',
    )(function)
    function = click.option(
        "--assembly-ids",
        default=None,
        required=False,
        help='A list format specified array (e.g. ["sampleA", "SampleB"]) mapping sample names to input files when using --fasta-dir or config["fastadir"]. Note that this feature is undergoing evaluation.',
    )(function)
    function = click.option(
        "--latest-run",
        default=None,
        required=False,
        help="An internal logging parameter designating the timestamp of the current execution for tracking vOMIX-MEGA history inside the .vomix subdirectory. See documentation at https://vomix-mega.readthedocs.io/en/latest/.",
    )(function)
    function = click.option(
        "--splits",
        default=0,
        required=False,
        help="The number of parallel groups to partition data into during processing to reduce memory overhead (e.g. when running geNomad). A value of 0 disables partitioning.",
    )(function)
    function = click.option(
        "--keep-intermediates",
        is_flag=True,
        default=False,
        required=False,
        help="Specifies whether to retain substantial intermediate processing files such as fastp-cleaned raw FASTQ files prior to host decontamination.",
    )(function)
    function = click.option(
        "--setup-database",
        is_flag=True,
        default=False,
        required=False,
        help="Determines whether to initialize or update databases when executing modules other than vomix setup-database. Existing databases will not be re-installed unless forced using Snakemake parameters like --forcerun or -F.",
    )(function)
    function = click.option(
        "--max-cores",
        default=4,
        required=False,
        help="The maximum number of CPU cores allocated dynamically across parallel Snakemake tasks. This feature is distinct from execution parameters like -j or -n and is currently in development.",
    )(function)
    function = click.option(
        "--NCBI-email",
        default=None,
        required=False,
        help="The user email address provided to NCBI E-utilities for data retrieval and download verification.",
    )(function)
    function = click.option(
        "--NCBI-API-key",
        default=None,
        required=False,
        help="The NCBI API key required for higher throughput data retrieval. For setup details visit https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities.",
    )(function)
    function = click.option(
        "--custom-config",
        default=None,
        required=False,
        help="Path to your custom config.yml",
    )(function)
    return function


# common snakemake options decorator
def snakemake_options(function):
    function = click.option(
        "--dry-run",
        "--dryrun",
        "-n",
        required=False,
        default=False,
        flag_value=True,
        help="Do not execute anything, and display what would be done. If you have a very large workflow, use --dry-run --quiet to just print a summary of the DAG of jobs. Snakemake local function, please run snakemake -h for more information.",
    )(function)
    function = click.option(
        "--forceall",
        "-F",
        required=False,
        default=False,
        flag_value=True,
        help="Force the execution of the selected (or the first) rule and all rules it is dependent on regardless of already created output. Snakemake local function, please run snakemake -h for more information.",
    )(function)
    function = click.option(
        "--configfile",
        default=None,
        required=False,
        help="Specify or overwrite the config file of the workflow (see the docs). Values specified in JSON or YAML format are available in the global config dictionary inside the workflow. Multiple files overwrite each other in the given order. Thereby missing keys in previous config files are extended by following configfiles. Note that this order also includes a config file defined in the workflow definition itself (which will come first). Snakemake local function, please run snakemake -h for more information.",
    )(function)
    function = click.option(
        "--unlock",
        required=False,
        default=False,
        flag_value=True,
        help="Remove a lock on the working directory. Snakemake local function, please run snakemake -h for more information.",
    )(function)
    function = click.option(
        "--cores",
        "-c",
        default=1,
        required=False,
        help='Use at most N CPU cores/jobs in parallel. If N is omitted or "all", the limit is set to the number of available CPU cores. In case of cluster/cloud execution, this argument sets the maximum number of cores requested from the cluster or cloud scheduler. (See https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#resources-remote-execution for more info.) This number is available to rules via workflow.cores. Snakemake local function, please run snakemake -h for more information.',
    )(function)
    function = click.option(
        "--jobs",
        "-j",
        default=4,
        required=False,
        help='Use at most N CPU cluster/cloud jobs in parallel. For local execution this is an alias for --cores (it is though recommended to use --cores in that case). Note: Set to "unlimited" to allow any number of parallel jobs. Snakemake local function, please run snakemake -h for more information.',
    )(function)
    function = click.option(
        "--latency-wait",
        default=20,
        required=False,
        help="Wait given seconds if an output file of a job is not present after the job finished. This helps if your filesystem suffers from latency. Snakemake local function, please run snakemake -h for more information.",
    )(function)
    function = click.option(
        "--rerun-incomplete",
        "-ri",
        required=False,
        default=False,
        flag_value=True,
        help="Re-run all jobs the output of which is recognized as incomplete. Snakemake local function, please run snakemake -h for more information.",
    )(function)
    function = click.option(
        "--rerun-triggers",
        required=False,
        default=None,
        help='Define what triggers the rerunning of a job. By default, all triggers are used, which guarantees that results are consistent with the workflow code and configuration. If you rather prefer the traditional way of just considering file modification dates, use "--rerun- trigger mtime". (default: code input mtime params software-env). Snakemake local function, please run snakemake -h for more information.',
    )(function)
    function = click.option(
        "--sdm",
        required=False,
        default=None,
        help="Specify software environment deployment method. Snakemake local function, please run snakemake -h for more information.",
    )(function)
    function = click.option(
        "--executor",
        "-e",
        required=False,
        default=None,
        help="Specify a custom executor, available via an executor plugin: snakemake_executor_<name>. Snakemake local function, please run snakemake -h for more information.",
    )(function)
    function = click.option(
        "--quiet",
        required=False,
        default=False,
        flag_value=True,
        help='Do not output certain information. If used without arguments, do not output any progress or rule information. Defining "all" results in no information being printed at all. Snakemake local function, please run snakemake -h for more information.',
    )(function)
    function = click.option(
        "--snakemake-args",
        required=False,
        default=None,
        help='Additional arguments to pass to the native snakemake command. (default: "")',
    )(function)

    return function


def setOptions(
    module_obj,
    workdir,
    outdir,
    datadir,
    samplelist,
    fasta,
    fastadir,
    sample_name,
    assembly_ids,
    latest_run,
    splits,
    keep_intermediates,
    setup_database,
    max_cores,
    ncbi_email,
    ncbi_api_key,
    custom_config,
):
    module_obj.workdir = workdir
    module_obj.outdir = outdir
    module_obj.datadir = datadir
    module_obj.samplelist = samplelist
    module_obj.fasta = fasta
    module_obj.fastadir = fastadir
    module_obj.sample_name = sample_name
    module_obj.assembly_ids = assembly_ids
    module_obj.latest_run = latest_run
    module_obj.splits = splits
    module_obj.keep_intermediates = keep_intermediates
    module_obj.setup_database = setup_database
    module_obj.max_cores = max_cores
    module_obj.ncbi_email = ncbi_email
    module_obj.ncbi_api_key = ncbi_api_key
    module_obj.custom_config = custom_config

    return module_obj


# vOMIX-MEGA command line interface


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def cli():
    """
    vOMIX-MEGA - an ultra-fast end-to-end pipeline for terabyte-scale viral metagenomics analysis.
    """


# Preprocess Module


@cli.command(
    "preprocess",
    context_settings={"ignore_unknown_options": True},
    short_help="Run the Pre-processing module",
)
@common_options
@click.option(
    "--dwnld-only",
    default=False,
    required=False,
    help="Flag specifying whether to restrict execution exclusively to remote SRA file downloads from NCBI without proceeding to pre-processing steps. (default: False)",
)
@click.option(
    "--decontam-host",
    is_flag=True,
    flag_value=True,
    default=False,
    required=False,
    help="Flag specifying whether to perform host decontamination post fastp quality trimming. Host indexes are provided via the --hostile-index-name CLI flag. (default: False)",
)
@click.option(
    "--dwnld-params",
    required=False,
    default=None,
    help='Optional configuration parameters used during raw FASTQ retrieval from NCBI via entrez-fetch. (default: "")',
)
@click.option(
    "--pigz-params",
    required=False,
    default=None,
    help='Execution parameters passed directly to pigz for multi-threaded compression of downloaded SRA FASTQ data. (default: "")',
)
@click.option(
    "--fastp-params",
    required=False,
    default=None,
    help='Additional runtime arguments supplied to the fastp quality control engine. (default: "")',
)
@click.option(
    "--hostile-params",
    required=False,
    default=None,
    help='Additional runtime arguments supplied to the Hostile host decontamination module. (default: "")',
)
@click.option(
    "--hostile-aligner",
    required=False,
    default=None,
    help='The short-read alignment backend algorithm employed for host decontamination (either "bowtie2" or "minimap2"). Note that minimap2 is not recommended for short-read sequencing datasets. (default: "bowtie2")',
)
@click.option(
    "--hostile-aligner-params",
    required=False,
    default=None,
    help='Additional runtime arguments supplied directly to the selected Hostile alignment tool. (default: "")',
)
@click.option(
    "--hostile-index-name",
    required=False,
    default=None,
    help='The name identifier of pre-built Hostile indices. Verify availability within the installed Hostile tool version; consult https://github.com/bede/hostile for details. (default: "human-t2t-hla")',
)
@click.option(
    "--hostile-index-db",
    required=False,
    default=None,
    help='The directory path where the Hostile database is installed or will be downloaded. Defaults to the Snakemake base directory under workflow/databases. (default: "database/hostile")',
)
@snakemake_options
def run_preprocess(
    workdir,
    outdir,
    datadir,
    samplelist,
    custom_config,
    fasta,
    fastadir,
    sample_name,
    assembly_ids,
    latest_run,
    splits,
    dwnld_only,
    keep_intermediates,
    setup_database,
    max_cores,
    ncbi_email,
    ncbi_api_key,
    decontam_host,
    dwnld_params,
    pigz_params,
    fastp_params,
    hostile_params,
    hostile_aligner,
    hostile_aligner_params,
    hostile_index_name,
    hostile_index_db,
    dry_run,
    forceall,
    configfile,
    unlock,
    cores,
    jobs,
    latency_wait,
    rerun_incomplete,
    rerun_triggers,
    sdm,
    executor,
    quiet,
    snakemake_args,
):
    cwd = os.path.abspath(os.getcwd())
    logging.info("Running module: preprocess")
    logging.info(f"User directory: {cwd}")
    logging.info(
        f"decontamHost: {decontam_host}, outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}"
    )
    # logging.info("Run with snakemake flags: " + "dry_run:" + str(dry_run) + ", forceall:" + str(forceall) + ", configfile:" + str(configfile))

    module_obj = PreProcessingModule()
    module_obj.name = "preprocess"
    # Set the attributes of the module object
    module_obj = setOptions(
        module_obj=module_obj,
        workdir=workdir,
        outdir=outdir,
        datadir=datadir,
        samplelist=samplelist,
        fasta=fasta,
        fastadir=fastadir,
        sample_name=sample_name,
        assembly_ids=assembly_ids,
        latest_run=latest_run,
        splits=splits,
        keep_intermediates=keep_intermediates,
        setup_database=setup_database,
        max_cores=max_cores,
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        custom_config=custom_config,
    )

    module_obj.decontam_host = decontam_host

    # optional params
    if dwnld_only:
        module_obj.dwnld_only = dwnld_only
        module_obj.hasOptions = True
    if dwnld_params:
        module_obj.dwnld_params = dwnld_params
        module_obj.hasOptions = True
    if pigz_params:
        module_obj.pigz_params = pigz_params
        module_obj.hasOptions = True
    if fastp_params:
        module_obj.fastp_params = fastp_params
        module_obj.hasOptions = True
    if hostile_params:
        module_obj.hostile_params = hostile_params
        module_obj.hasOptions = True
    if hostile_aligner:
        module_obj.hostile_aligner = hostile_aligner
        module_obj.hasOptions = True
    if hostile_aligner_params:
        module_obj.aligner_params = hostile_aligner_params
        module_obj.hasOptions = True
    if hostile_index_name:
        module_obj.hostile_index_name = hostile_index_name
        module_obj.hasOptions = True

    # snakemake options
    snakemake_obj = SnakemakeFlags(
        dry_run,
        forceall,
        configfile,
        unlock,
        cores,
        jobs,
        latency_wait,
        rerun_incomplete,
        rerun_triggers,
        sdm,
        executor,
        quiet,
        snakemake_args,
    )

    vomix_actions_instance = vomix_actions()
    vomix_actions_instance.run_module("preprocess", module_obj, snakemake_obj)
    logging.info(END_MODULE_RUN_LOG)


# Assembly Module


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
    help='The primary assembly tool engine selected for the assembly module (spades support is currently in development). (default: "megahit")',
)
@click.option(
    "--megahit-min-len",
    required=False,
    default=300,
    help="The minimum contig length threshold used during the MEGAHIT assembly filter steps. (default: 300)",
)
@click.option(
    "--megahit-params",
    required=False,
    default=None,
    help='Additional runtime parameters supplied directly to the MEGAHIT execution pipeline. (default: "--prune-level 3")',
)
@click.option(
    "--spades-params",
    required=False,
    default=None,
    help='Additional runtime parameters supplied directly to the SPAdes metagenomic assembler execution line. (default: "--meta")',
)
@click.option(
    "--spades-memory",
    required=False,
    default=250,
    help="The upper threshold of RAM memory (in gigabytes) allocated for SPAdes assembly execution. (default: 250)",
)
@snakemake_options
def run_assembly(
    workdir,
    outdir,
    datadir,
    samplelist,
    custom_config,
    fasta,
    fastadir,
    sample_name,
    assembly_ids,
    latest_run,
    splits,
    keep_intermediates,
    setup_database,
    max_cores,
    ncbi_email,
    ncbi_api_key,
    assembler,
    megahit_min_len,
    megahit_params,
    spades_params,
    spades_memory,
    dry_run,
    forceall,
    configfile,
    unlock,
    cores,
    jobs,
    latency_wait,
    rerun_incomplete,
    rerun_triggers,
    sdm,
    executor,
    quiet,
    snakemake_args,
):
    logging.info("Running module: assembly")
    logging.info(
        f"assembler: {assembler}, outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}"
    )

    module_obj = AssemblyCoAssemblyModule()
    module_obj.name = "assembly"
    # Set the attributes of the module object
    module_obj = setOptions(
        module_obj=module_obj,
        workdir=workdir,
        outdir=outdir,
        datadir=datadir,
        samplelist=samplelist,
        fasta=fasta,
        fastadir=fastadir,
        sample_name=sample_name,
        assembly_ids=assembly_ids,
        latest_run=latest_run,
        splits=splits,
        keep_intermediates=keep_intermediates,
        setup_database=setup_database,
        max_cores=max_cores,
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        custom_config=custom_config,
    )

    module_obj.assembler = assembler

    if megahit_min_len:
        module_obj.megahit_min_len = megahit_min_len
        module_obj.hasOptions = True
    if megahit_params:
        module_obj.megahit_params = megahit_params
        module_obj.hasOptions = True
    if spades_params:
        module_obj.spades_params = spades_params
        module_obj.hasOptions = True
    if spades_memory:
        module_obj.spades_memory = spades_memory
        module_obj.hasOptions = True

    snakemake_obj = SnakemakeFlags(
        dry_run,
        forceall,
        configfile,
        unlock,
        cores,
        jobs,
        latency_wait,
        rerun_incomplete,
        rerun_triggers,
        sdm,
        executor,
        quiet,
        snakemake_args,
    )

    vomix_actions_instance = vomix_actions()
    vomix_actions_instance.run_module("assembler", module_obj, snakemake_obj)
    logging.info(END_MODULE_RUN_LOG)


# Viral Identify Module


@cli.command(
    "viral-identify",
    context_settings={"ignore_unknown_options": True},
    short_help="Run the Viral Identify module",
)
@common_options
@click.option(
    "--contig-min-len",
    required=False,
    default=0,
    help="The absolute minimum length constraint for contig inclusion within the viral-identify module; shorter sequences are purged from analysis. (default: 0)",
)
@click.option(
    "--genomad-db",
    required=False,
    default=None,
    help='The directory path where the geNomad database is installed or will be downloaded. Defaults to the Snakemake base directory under workflow/databases. (default: "database/genomad")',
)
@click.option(
    "--genomad-min-len",
    required=False,
    default=None,
    help="The minimum contig length evaluated by geNomad; sequences falling below this parameter are excluded from classification steps. (default: 10000)",
)
@click.option(
    "--genomad-params",
    required=False,
    default=None,
    help='Additional runtime command line arguments supplied to geNomad execution within the viral-identify framework. (default: "--enable-score-calibration --relaxed")',
)
@click.option(
    "--genomad-cutoff",
    required=False,
    default=None,
    help="The minimal numeric confidence threshold required by geNomad to classify a contig sequence as viral. (default: 0.7)",
)
@click.option(
    "--genomad-cutoff-s",
    required=False,
    default=None,
    help="The minimum confidence threshold applied during geNomad secondary filtering. Setting this to 0 bypasses the secondary filter pipeline entirely.",
)
@click.option(
    "--checkv-original",
    is_flag=True,
    flag_value=True,
    required=False,
    default=False,
    help="Flag allowing execution of standard CheckV instead of the lower memory high-efficiency CheckV-PyHMMER implementation. (default: False)",
)
@click.option(
    "--checkv-params",
    required=False,
    default=None,
    help='Additional operational arguments supplied directly to the CheckV pipeline execution. (default: "")',
)
@click.option(
    "--checkv-database",
    required=False,
    default=None,
    help='The directory path where the CheckV database is installed or will be downloaded. (default: "database/checkv")',
)
@click.option(
    "--clustering-fast",
    is_flag=True,
    flag_value=True,
    required=False,
    default=True,
    help="Flag triggering an accelerated MEGABlast-based clustering protocol optimized for viral operational taxonomic unit (vOTU) compilation. (default: True)",
)
@click.option(
    "--cdhit-params",
    required=False,
    default=None,
    help='Additional operational runtime values supplied directly to the CD-HIT clustering utility. (default: "-c 0.95 -aS 0.85 -d 400 -M 0 -n 5")',
)
@click.option(
    "--vOTU-ani",
    required=False,
    default=None,
    help="The average nucleotide identity (ANI) clustering percentage threshold used during fast clustering workflows (default 95% per MIUViG guidelines: https://doi.org/10.1038/nbt.4306). (default: 95)",
)
@click.option(
    "--vOTU-targetcov",
    required=False,
    default=None,
    help="The minimum target coverage alignment coverage percentage used during fast clustering workflows (default 85% per MIUViG guidelines: https://doi.org/10.1038/nbt.4306). (default: 85)",
)
@click.option(
    "--vOTU-querycov",
    required=False,
    default=None,
    help="The target query alignment coverage percentage criteria implemented within the fast clustering protocol. (default: 0)",
)
@snakemake_options
def run_viral_identify(
    workdir,
    outdir,
    datadir,
    samplelist,
    custom_config,
    fasta,
    fastadir,
    sample_name,
    assembly_ids,
    latest_run,
    splits,
    keep_intermediates,
    setup_database,
    max_cores,
    ncbi_email,
    ncbi_api_key,
    contig_min_len,
    genomad_db,
    genomad_min_len,
    genomad_params,
    genomad_cutoff,
    genomad_cutoff_s,
    checkv_original,
    checkv_params,
    checkv_database,
    clustering_fast,
    cdhit_params,
    votu_ani,
    votu_targetcov,
    votu_querycov,
    dry_run,
    forceall,
    configfile,
    unlock,
    cores,
    jobs,
    latency_wait,
    rerun_incomplete,
    rerun_triggers,
    sdm,
    executor,
    quiet,
    snakemake_args,
):
    logging.info("Running module: viral-identify")
    logging.info(f"outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")

    module_obj = ViralIdentifyModule()
    module_obj.name = "viral-identify"
    # Set the attributes of the module object
    module_obj = setOptions(
        module_obj=module_obj,
        workdir=workdir,
        outdir=outdir,
        datadir=datadir,
        samplelist=samplelist,
        fasta=fasta,
        fastadir=fastadir,
        sample_name=sample_name,
        assembly_ids=assembly_ids,
        latest_run=latest_run,
        splits=splits,
        keep_intermediates=keep_intermediates,
        setup_database=setup_database,
        max_cores=max_cores,
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        custom_config=custom_config,
    )

    if contig_min_len:
        module_obj.contig_min_len = contig_min_len
        module_obj.hasOptions = True
    if genomad_db:
        module_obj.genomad_db = genomad_db
        module_obj.hasOptions = True
    if genomad_min_len:
        module_obj.genomad_min_len = genomad_min_len
        module_obj.hasOptions = True
    if genomad_params:
        module_obj.genomad_params = genomad_params
        module_obj.hasOptions = True
    if genomad_cutoff:
        module_obj.genomad_cutoff = genomad_cutoff
        module_obj.hasOptions = True
    if genomad_cutoff_s:
        module_obj.genomad_cutoff_s = genomad_cutoff_s
        module_obj.hasOptions = True
    if checkv_original:
        module_obj.checkv_original = checkv_original
        module_obj.hasOptions = True
    if checkv_params:
        module_obj.checkv_params = checkv_params
        module_obj.hasOptions = True
    if checkv_database:
        module_obj.checkv_database = checkv_database
        module_obj.hasOptions = True
    if clustering_fast:
        module_obj.clustering_fast = clustering_fast
        module_obj.hasOptions = True
    if cdhit_params:
        module_obj.cdhit_params = cdhit_params
        module_obj.hasOptions = True
    if votu_ani:
        module_obj.vOTU_ani = votu_ani
        module_obj.hasOptions = True
    if votu_targetcov:
        module_obj.vOTU_targetcov = votu_targetcov
        module_obj.hasOptions = True
    if votu_querycov:
        module_obj.vOTU_querycov = votu_querycov
        module_obj.hasOptions = True

    snakemake_obj = SnakemakeFlags(
        dry_run,
        forceall,
        configfile,
        unlock,
        cores,
        jobs,
        latency_wait,
        rerun_incomplete,
        rerun_triggers,
        sdm,
        executor,
        quiet,
        snakemake_args,
    )

    vomix_actions_instance = vomix_actions()
    vomix_actions_instance.run_module("viral-identify", module_obj, snakemake_obj)
    logging.info(END_MODULE_RUN_LOG)


# Viral Benchmark Module


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
    help='The directory path where the PhaBox2 database is installed or will be downloaded. Defaults to the Snakemake base directory under workflow/databases. (default: "database/phabox_db_v2")',
)
@click.option(
    "--phabox2-db-name",
    required=False,
    default=None,
    help='The designated database name or identifier file package required for the PhaBox2 classification tool execution. (default: "phabox_db_v2")',
)
@click.option(
    "--phabox2-db-baselink",
    required=False,
    default=None,
    help='The primary remote server base link URL used to fetch and download resource updates for the PhaBox2 database. (default: "https://github.com/KennthShang/PhaBOX/releases/download/v2")',
)
@click.option(
    "--genomad-db",
    required=False,
    default=None,
    help='The directory path where the geNomad database is installed or will be downloaded. Defaults to the Snakemake base directory under workflow/databases. (default: "database/genomad")',
)
@click.option(
    "--virsorter2-db",
    required=False,
    default=None,
    help='The directory path where the VirSorter2 database is installed or will be downloaded. Defaults to the Snakemake base directory under workflow/databases. (default: "database/virsorter2")',
)
@click.option(
    "--vibrant-db",
    required=False,
    default=None,
    help='The directory path where the VIBRANT database is installed or will be downloaded. Defaults to the Snakemake base directory under workflow/databases. (default: "database/virsorter2")',
)
@click.option(
    "--contig-min-len",
    required=False,
    default=0,
    help="The absolute minimum length constraint for contig inclusion within the viral-identify module; shorter sequences are purged from analysis. (default: 0)",
)
@click.option(
    "--genomad-min-len",
    required=False,
    default=1000,
    help="The minimum contig length evaluated by geNomad; sequences falling below this parameter are excluded from classification steps. (default: 10000)",
)
@click.option(
    "--genomad-params",
    required=False,
    default=None,
    help='Additional runtime command line arguments supplied to geNomad execution within the viral-identify framework. (default: "--enable-score-calibration --relaxed")',
)
@click.option(
    "--genomad-cutoff",
    required=False,
    default=0.7,
    help="The minimal numeric confidence threshold required by geNomad to classify a contig sequence as viral. (default: 0.7)",
)
@click.option(
    "--genomad-cutoff-s",
    required=False,
    default=0,
    help="The minimum confidence threshold applied during geNomad secondary filtering. Setting this to 0 bypasses the secondary filter pipeline entirely. (default: 0)",
)
@click.option(
    "--dvf-min-len",
    required=False,
    default=1500,
    help="The lower bound contig length cut-off implemented during DeepVirFinder evaluation; shorter contigs are ignored. (default: 1500)",
)
@click.option(
    "--phamer-min-len",
    required=False,
    default=2000,
    help="The lower bound contig length cut-off implemented during PhaMer evaluation; shorter sequences are omitted. (default: 2000)",
)
@click.option(
    "--dvf-params",
    required=False,
    default=None,
    help='Additional system parameters passed directly to the DeepVirFinder tool environment. (default: "")',
)
@click.option(
    "--phamer-params",
    required=False,
    default=None,
    help='Additional system parameters passed directly to the PhaMer tool environment. (default: "")',
)
@click.option(
    "--virsorter2-params",
    required=False,
    default=None,
    help='Additional system parameters passed directly to the VirSorter2 tool environment. (default: "")',
)
@click.option(
    "--vf-params",
    required=False,
    default=None,
    help='Additional system parameters passed directly to the VirFinder tool environment. (default: "")',
)
@click.option(
    "--seeker-params",
    required=False,
    default=None,
    help='Additional system parameters passed directly to the Seeker tool environment. (default: "")',
)
@click.option(
    "--PPR-params",
    required=False,
    default=None,
    help='Additional system parameters passed directly to the PPR-META tool environment. (default: "")',
)
@click.option(
    "--dvf-cutoff",
    required=False,
    default=0.7,
    help="The minimal confidence score metric required by DeepVirFinder to classify a sequence as viral. (default: 0.7)",
)
@click.option(
    "--dvf-pval",
    required=False,
    default=0.05,
    help="The maximum critical p-value threshold permitted by DeepVirFinder to confirm a sequence classification as viral. (default: 0.05)",
)
@click.option(
    "--phamer-pred",
    required=False,
    default=None,
    help='The taxonomic classification category targeted by PhaMer prediction routines. (default: "phage")',
)
@click.option(
    "--phamer-cutoff",
    required=False,
    default=0,
    help="The minimal confidence threshold value required for a positive viral determination within the PhaMer algorithm. (default: 0)",
)
@click.option(
    "--vf-cutoff",
    required=False,
    default=0,
    help="The minimal confidence threshold value required for a positive viral determination within the VirFinder algorithm. (default: 0)",
)
@click.option(
    "--virsorter2-cutoff",
    required=False,
    default=0,
    help="The minimal confidence threshold value required for a positive viral determination within the VirFinder algorithm. (default: 0)",
)
@click.option(
    "--seeker-cutoff",
    required=False,
    default=0,
    help="The minimal confidence threshold value required for a positive viral determination within the Seeker algorithm. (default: 0)",
)
@click.option(
    "--ppr-cutoff",
    required=False,
    default=0,
    help="The minimal confidence threshold value required for a positive viral determination within the PPR-META algorithm. (default: 0)",
)
@click.option(
    "--vibrant-cutoff",
    required=False,
    default=0,
    help="The minimal confidence threshold value required for a positive viral determination within the VIBRANT algorithm. (default: 0)",
)
@click.option(
    "--checkv-original",
    is_flag=True,
    flag_value=True,
    required=False,
    default=False,
    help="Flag allowing execution of standard CheckV instead of the lower memory high-efficiency CheckV-PyHMMER implementation. (default: False)",
)
@click.option(
    "--checkv-params",
    required=False,
    default=None,
    help='Additional operational arguments supplied directly to the CheckV pipeline execution. (default: "")',
)
@click.option(
    "--checkv-database",
    required=False,
    default=None,
    help='The directory path where the CheckV database is installed or will be downloaded. (default: "database/checkv")',
)
@click.option(
    "--clustering-fast",
    is_flag=True,
    flag_value=True,
    required=False,
    default=True,
    help="Flag triggering an accelerated MEGABlast-based clustering protocol optimized for viral operational taxonomic unit (vOTU) compilation. (default: True)",
)
@click.option(
    "--cdhit-params",
    required=False,
    default=None,
    help='Additional operational runtime values supplied directly to the CD-HIT clustering utility. (default: "-c 0.95 -aS 0.85 -d 400 -M 0 -n 5")',
)
@click.option(
    "--vOTU-ani",
    required=False,
    default=None,
    help="The average nucleotide identity (ANI) clustering percentage threshold used during fast clustering workflows (default 95% per MIUViG guidelines: https://doi.org/10.1038/nbt.4306). (default: 95)",
)
@click.option(
    "--vOTU-targetcov",
    required=False,
    default=None,
    help="The minimum target coverage alignment coverage percentage used during fast clustering workflows (default 85% per MIUViG guidelines: https://doi.org/10.1038/nbt.4306). (default: 85)",
)
@click.option(
    "--vOTU-querycov",
    required=False,
    default=None,
    help="The target query alignment coverage percentage criteria implemented within the fast clustering protocol. (default: 0)",
)
@snakemake_options
def run_viral_benchmark(
    workdir,
    outdir,
    datadir,
    samplelist,
    custom_config,
    fasta,
    fastadir,
    sample_name,
    assembly_ids,
    latest_run,
    splits,
    keep_intermediates,
    setup_database,
    max_cores,
    ncbi_email,
    ncbi_api_key,
    phabox2_db,
    phabox2_db_name,
    phabox2_db_baselink,
    genomad_db,
    virsorter2_db,
    vibrant_db,
    contig_min_len,
    genomad_min_len,
    genomad_params,
    genomad_cutoff,
    genomad_cutoff_s,
    dvf_min_len,
    phamer_min_len,
    dvf_params,
    phamer_params,
    virsorter2_params,
    vf_params,
    seeker_params,
    PPR_params,
    dvf_cutoff,
    dvf_pval,
    phamer_pred,
    phamer_cutoff,
    vf_cutoff,
    dry_run,
    forceall,
    configfile,
    unlock,
    cores,
    jobs,
    latency_wait,
    rerun_incomplete,
    rerun_triggers,
    sdm,
    executor,
    quiet,
    snakemake_args,
):
    logging.info("Running module: viral-benchmark")
    logging.info(f"outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")

    module_obj = ViralBenchmarkModule()
    module_obj.name = "viral-benchmark"
    # Set the attributes of the module object
    module_obj = setOptions(
        module_obj=module_obj,
        workdir=workdir,
        outdir=outdir,
        datadir=datadir,
        samplelist=samplelist,
        fasta=fasta,
        fastadir=fastadir,
        sample_name=sample_name,
        assembly_ids=assembly_ids,
        latest_run=latest_run,
        splits=splits,
        keep_intermediates=keep_intermediates,
        setup_database=setup_database,
        max_cores=max_cores,
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        custom_config=custom_config,
    )

    if contig_min_len:
        module_obj.contig_min_len = contig_min_len
        module_obj.hasOptions = True

    if phabox2_db:
        module_obj.PhaBox2_db = phabox2_db
        module_obj.hasOptions = True
    if phabox2_db_name:
        module_obj.phabox2_db_name = phabox2_db_name
        module_obj.hasOptions = True
    if phabox2_db_baselink:
        module_obj.phabox2_db_baselink = phabox2_db_baselink
        module_obj.hasOptions = True
    if genomad_db:
        module_obj.genomad_db = genomad_db
        module_obj.hasOptions = True
    if virsorter2_db:
        module_obj.virsorter2_db = virsorter2_db
        module_obj.hasOptions = True
    if vibrant_db:
        module_obj.vibrant_db = vibrant_db
        module_obj.hasOptions = True
    if genomad_min_len:
        module_obj.genomad_min_len = genomad_min_len
        module_obj.hasOptions = True
    if genomad_params:
        module_obj.genomad_params = genomad_params
        module_obj.hasOptions = True
    if genomad_cutoff:
        module_obj.genomad_cutoff = genomad_cutoff
        module_obj.hasOptions = True
    if genomad_cutoff_s:
        module_obj.genomad_cutoff_s = genomad_cutoff_s
        module_obj.hasOptions = True
    if dvf_min_len:
        module_obj.dvf_min_len = dvf_min_len
        module_obj.hasOptions = True
    if phamer_min_len:
        module_obj.phamer_min_len = phamer_min_len
        module_obj.hasOptions = True
    if dvf_params:
        module_obj.dvf_params = dvf_params
        module_obj.hasOptions = True
    if phamer_params:
        module_obj.phamer_params = phamer_params
        module_obj.hasOptions = True
    if virsorter2_params:
        module_obj.virsorter2_params = virsorter2_params
        module_obj.hasOptions = True
    if vf_params:
        module_obj.vf_params = vf_params
        module_obj.hasOptions = True
    if seeker_params:
        module_obj.seeker_params = seeker_params
        module_obj.hasOptions = True
    if PPR_params:
        module_obj.PPR_params = PPR_params
        module_obj.hasOptions = True
    if dvf_cutoff:
        module_obj.dvf_cutoff = dvf_cutoff
        module_obj.hasOptions = True
    if dvf_pval:
        module_obj.dvf_pval = dvf_pval
        module_obj.hasOptions = True
    if phamer_pred:
        module_obj.phamer_pred = phamer_pred
        module_obj.hasOptions = True
    if phamer_cutoff:
        module_obj.phamer_cutoff = phamer_cutoff
        module_obj.hasOptions = True
    if vf_cutoff:
        module_obj.vf_cutoff = vf_cutoff
        module_obj.hasOptions = True

    snakemake_obj = SnakemakeFlags(
        dry_run,
        forceall,
        configfile,
        unlock,
        cores,
        jobs,
        latency_wait,
        rerun_incomplete,
        rerun_triggers,
        sdm,
        executor,
        quiet,
        snakemake_args,
    )

    vomix_actions_instance = vomix_actions()
    vomix_actions_instance.run_module("viral-benchmark", module_obj, snakemake_obj)
    logging.info(END_MODULE_RUN_LOG)


# Viral Taxonomy Module


@cli.command(
    "viral-taxonomy",
    context_settings={"ignore_unknown_options": True},
    short_help="Run the Viral Taxonomy module",
)
@common_options
@click.option(
    "--PhaBox2-db",
    required=False,
    default=None,
    help='The directory path where the PhaBox2 database is installed or will be downloaded. Defaults to the Snakemake base directory under workflow/databases. (default: "database/phabox_db_v2")',
)
@click.option(
    "--phabox2-db-name",
    required=False,
    default=None,
    help='The designated database name or identifier file package required for the PhaBox2 classification tool execution. (default: "phabox_db_v2")',
)
@click.option(
    "--phabox2-db-baselink",
    required=False,
    default=None,
    help='The primary remote server base link URL used to fetch and download resource updates for the PhaBox2 database. (default: "https://github.com/KennthShang/PhaBOX/releases/download/v2")',
)
@click.option(
    "--phagcn-min-len",
    required=False,
    default=None,
    help="The minimum allowed contig length for evaluation using the PhaGCN taxonomic assignment engine. (default: 1000)",
)
@click.option(
    "--phagcn-params",
    required=False,
    default=None,
    help='Additional operational arguments passed to the PhaGCN classification instance. (default: "")',
)
@click.option(
    "--genomad-db",
    required=False,
    default=None,
    help='Path to geNomad database directory || default: "workflow/database/genomad" [STR]',
)
@click.option(
    "--genomad-params-tax",
    required=False,
    default=None,
    help='Additional operational configurations passed to geNomad during viral taxonomic assignment. (default: "--enable-score-calibration --relaxed")',
)
@snakemake_options
def run_viral_taxonomy(
    workdir,
    outdir,
    datadir,
    samplelist,
    custom_config,
    fasta,
    fastadir,
    sample_name,
    assembly_ids,
    latest_run,
    splits,
    keep_intermediates,
    setup_database,
    max_cores,
    ncbi_email,
    ncbi_api_key,
    phabox2_db,
    phabox2_db_name,
    phabox2_db_baselink,
    phagcn_min_len,
    phagcn_params,
    genomad_db,
    genomad_params_tax,
    dry_run,
    forceall,
    configfile,
    unlock,
    cores,
    jobs,
    latency_wait,
    rerun_incomplete,
    rerun_triggers,
    sdm,
    executor,
    quiet,
    snakemake_args,
):
    logging.info("Running module: viral-taxonomy")
    logging.info(f"fasta: {fasta}, outdir: {outdir}")

    module_obj = ViralTaxonomyModule()
    module_obj.name = "viral-taxonomy"
    # Set the attributes of the module object
    module_obj = setOptions(
        module_obj=module_obj,
        workdir=workdir,
        outdir=outdir,
        datadir=datadir,
        samplelist=samplelist,
        fasta=fasta,
        fastadir=fastadir,
        sample_name=sample_name,
        assembly_ids=assembly_ids,
        latest_run=latest_run,
        splits=splits,
        keep_intermediates=keep_intermediates,
        setup_database=setup_database,
        max_cores=max_cores,
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        custom_config=custom_config,
    )

    if phabox2_db:
        module_obj.PhaBox2_db = phabox2_db
        module_obj.hasOptions = True
    if phabox2_db_name:
        module_obj.phabox2_db_name = phabox2_db_name
        module_obj.hasOptions = True
    if phabox2_db_baselink:
        module_obj.phabox2_db_baselink = phabox2_db_baselink
        module_obj.hasOptions = True
    if phagcn_min_len:
        module_obj.phagcn_min_len = phagcn_min_len
        module_obj.hasOptions = True
    if phagcn_params:
        module_obj.phagcn_params = phagcn_params
        module_obj.hasOptions = True
    if genomad_db:
        module_obj.genomad_db = genomad_db
        module_obj.hasOptions = True
    if genomad_params_tax:
        module_obj.genomad_params_tax = genomad_params_tax
        module_obj.hasOptions = True

    snakemake_obj = SnakemakeFlags(
        dry_run,
        forceall,
        configfile,
        unlock,
        cores,
        jobs,
        latency_wait,
        rerun_incomplete,
        rerun_triggers,
        sdm,
        executor,
        quiet,
        snakemake_args,
    )

    vomix_actions_instance = vomix_actions()
    vomix_actions_instance.run_module("viral-taxonomy", module_obj, snakemake_obj)
    logging.info(END_MODULE_RUN_LOG)


# Viral Host Module


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
    help='The directory path where the PhaBox2 database is installed or will be downloaded. Defaults to the Snakemake base directory under workflow/databases. (default: "database/phabox_db_v2")',
)
@click.option(
    "--phabox2-db-name",
    required=False,
    default=None,
    help='The designated database name or identifier file package required for the PhaBox2 classification tool execution. (default: "phabox_db_v2")',
)
@click.option(
    "--phabox2-db-baselink",
    required=False,
    default=None,
    help='The primary remote server base link URL used to fetch and download resource updates for the PhaBox2 database. (default: "https://github.com/KennthShang/PhaBOX/releases/download/v2")',
)
@click.option(
    "--CHERRY-params",
    required=False,
    default=None,
    help='Additional execution parameters configured for the CHERRY host prediction algorithm. (default: "")',
)
@click.option(
    "--PhaTYP-params",
    required=False,
    default=None,
    help='Additional custom parameters passed directly to the PhaTYP lifestyle prediction module. (default: "")',
)
@click.option(
    "--iphop-host",
    is_flag=True,
    flag_value=True,
    required=False,
    default=False,
    help="Flag indicating whether to perform iPHoP-based viral host prediction instead of CHERRY. Note that iPHoP requires a significantly larger database and high memory allocation. (default: False)",
)
@click.option(
    "--iphop-cutoff",
    required=False,
    default=None,
    help="The minimum confidence threshold required by iPHoP to assign a host classification profile to a viral sequence. (default: 90)",
)
@click.option(
    "--iphop-params",
    required=False,
    default=None,
    help='Additional configuration arguments supplied directly to the iPHoP platform interface. (default: "")',
)
@click.option(
    "--iphop-db",
    required=False,
    default=None,
    help='The directory path where the iPHoP database is installed or will be downloaded. Defaults to the Snakemake base directory under workflow/databases. (default: "database/iphop")',
)
@click.option(
    "--iphop-db-version",
    required=False,
    default=None,
    help='The version identifier for the iPHoP database. Verify database compatibility mapping at https://bitbucket.org/srouxjgi/iphop/src/main/#markdown-header-host-databases-and-versions. (default: "iPHoP_db_Aug23_rw")',
)
@click.option(
    "--iphop-db-basename",
    required=False,
    default=None,
    help='The primary base name of the iPHoP database. Verify compatibility criteria at https://bitbucket.org/srouxjgi/iphop/src/main/#markdown-header-host-databases-and-versions. (default: "Aug_2023_pub_rw")',
)
@snakemake_options
def run_viral_host(
    workdir,
    outdir,
    datadir,
    samplelist,
    custom_config,
    fasta,
    fastadir,
    sample_name,
    assembly_ids,
    latest_run,
    splits,
    keep_intermediates,
    setup_database,
    max_cores,
    ncbi_email,
    ncbi_api_key,
    phabox2_db,
    phabox2_db_name,
    phabox2_db_baselink,
    cherry_params,
    phatyp_params,
    iphop_host,
    iphop_cutoff,
    iphop_params,
    iphop_db,
    iphop_db_version,
    iphop_db_basename,
    dry_run,
    forceall,
    configfile,
    unlock,
    cores,
    jobs,
    latency_wait,
    rerun_incomplete,
    rerun_triggers,
    sdm,
    executor,
    quiet,
    snakemake_args,
):
    logging.info("Running module: viral-host")
    logging.info(f"fasta: {fasta}, outdir: {outdir}")

    module_obj = ViralHostModule()
    module_obj.name = "viral-host"
    # Set the attributes of the module object
    module_obj = setOptions(
        module_obj=module_obj,
        workdir=workdir,
        outdir=outdir,
        datadir=datadir,
        samplelist=samplelist,
        fasta=fasta,
        fastadir=fastadir,
        sample_name=sample_name,
        assembly_ids=assembly_ids,
        latest_run=latest_run,
        splits=splits,
        keep_intermediates=keep_intermediates,
        setup_database=setup_database,
        max_cores=max_cores,
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        custom_config=custom_config,
    )

    module_obj.iphop_host = iphop_host

    if phabox2_db:
        module_obj.PhaBox2_db = phabox2_db
        module_obj.hasOptions = True
    if phabox2_db_name:
        module_obj.phabox2_db_name = phabox2_db_name
        module_obj.hasOptions = True
    if phabox2_db_baselink:
        module_obj.phabox2_db_baselink = phabox2_db_baselink
        module_obj.hasOptions = True
    if cherry_params:
        module_obj.CHERRY_params = cherry_params
        module_obj.hasOptions = True
    if phatyp_params:
        module_obj.PhaTYP_params = phatyp_params
        module_obj.hasOptions = True
    if iphop_cutoff:
        module_obj.iphop_cutoff = iphop_cutoff
        module_obj.hasOptions = True
    if iphop_db:
        module_obj.iphop_db = iphop_db
        module_obj.hasOptions = True
    if iphop_db_version:
        module_obj.iphop_db_version = iphop_db_version
        module_obj.hasOptions = True
    if iphop_db_basename:
        module_obj.iphop_db_basename = iphop_db_basename
        module_obj.hasOptions = True

    snakemake_obj = SnakemakeFlags(
        dry_run,
        forceall,
        configfile,
        unlock,
        cores,
        jobs,
        latency_wait,
        rerun_incomplete,
        rerun_triggers,
        sdm,
        executor,
        quiet,
        snakemake_args,
    )

    vomix_actions_instance = vomix_actions()
    vomix_actions_instance.run_module("viral-host", module_obj, snakemake_obj)
    logging.info(END_MODULE_RUN_LOG)


# Viral Community Module


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
    help='Additional mapping or calculation flags passed to the CoverM coverage engine. (default: "--mapper minimap2-sr --min-read-percent-identity 95 --min-read-aligned-percent 75 --trim-min 10 --trim-max 90")',
)
@click.option(
    "--coverm-methods",
    required=False,
    default=None,
    help='The calculation metric outputs selected for CoverM. Modify only for pipeline testing or debugging. (default: "tpm rpkm")',
)
@snakemake_options
def run_viral_community(
    workdir,
    outdir,
    datadir,
    samplelist,
    custom_config,
    fasta,
    fastadir,
    sample_name,
    assembly_ids,
    latest_run,
    splits,
    keep_intermediates,
    setup_database,
    max_cores,
    ncbi_email,
    ncbi_api_key,
    dry_run,
    forceall,
    configfile,
    unlock,
    cores,
    jobs,
    latency_wait,
    rerun_incomplete,
    rerun_triggers,
    sdm,
    executor,
    quiet,
    snakemake_args,
    coverm_params,
    coverm_methods,
):
    logging.info("Running module: viral-community")
    logging.info(f"outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")

    module_obj = ViralCommunityModule()
    module_obj.name = "viral-community"
    # Set the attributes of the module object
    module_obj = setOptions(
        module_obj=module_obj,
        workdir=workdir,
        outdir=outdir,
        datadir=datadir,
        samplelist=samplelist,
        fasta=fasta,
        fastadir=fastadir,
        sample_name=sample_name,
        assembly_ids=assembly_ids,
        latest_run=latest_run,
        splits=splits,
        keep_intermediates=keep_intermediates,
        setup_database=setup_database,
        max_cores=max_cores,
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        custom_config=custom_config,
    )

    if coverm_params:
        module_obj.coverm_params = coverm_params
        module_obj.hasOptions = True
    if coverm_methods:
        module_obj.coverm_methods = coverm_methods
        module_obj.hasOptions = True

    snakemake_obj = SnakemakeFlags(
        dry_run,
        forceall,
        configfile,
        unlock,
        cores,
        jobs,
        latency_wait,
        rerun_incomplete,
        rerun_triggers,
        sdm,
        executor,
        quiet,
        snakemake_args,
    )

    vomix_actions_instance = vomix_actions()
    vomix_actions_instance.run_module("viral-community", module_obj, snakemake_obj)
    logging.info(END_MODULE_RUN_LOG)


# Viral Annotate Module


@cli.command(
    "viral-annotate",
    context_settings={"ignore_unknown_options": True},
    short_help="Run the Viral Annotate module",
)
@common_options
@click.option(
    "--eggNOG-params",
    required=False,
    default=None,
    help='Parameters for running eggNOG-mapper v2. See more at https://github.com/eggnogdb/eggnog-mapper/wiki || default: "-m diamond --hmm_evalue 0.001 --hmm_score 60 --query-cover 20 --subject-cover 20 --tax_scope auto --target_orthologs all --go_evidence non-electronic --report_orthologs" [INT]',
)
@click.option(
    "--PhaVIP-params",
    required=False,
    default=None,
    help='Minimum contig length to filter BEFORE viral identification || default: "" [STR]',
)
@click.option(
    "--metacerberus-db",
    required=False,
    default=None,
    help='The directory path where the MetaCerberus database is installed or will be downloaded. Defaults to the Snakemake base directory under workflow/databases. (default: "database/metacerberus")',
)
@click.option(
    "--metacerberus-setup-params",
    required=False,
    default=None,
    help='Operational configurations supplied to initialize build or index the MetaCerberus database environment. (default: "")',
)
@click.option(
    "--metacerberus-params",
    required=False,
    default=None,
    help='Parameters for running the MetaCerberus database. (default: "--hmm ALL")',
)
@click.option(
    "pharokka-db",
    required=False,
    default=None,
    help='The directory path where the pharokka database is installed or will be downloaded. Defaults to the Snakemake base directory under workflow/databases. (default: "database/pharokka")',
)
@click.option(
    "pharokka-params",
    required=False,
    default=None,
    help='Additional execution parameters passed directly to the pharokka bacteriophage annotation framework. (default: "-g prodigal-gv --meta")',
)
@snakemake_options
def run_viral_annotate(
    workdir,
    outdir,
    datadir,
    samplelist,
    custom_config,
    fasta,
    fastadir,
    sample_name,
    assembly_ids,
    latest_run,
    splits,
    keep_intermediates,
    setup_database,
    max_cores,
    ncbi_email,
    ncbi_api_key,
    eggnog_params,
    phavip_params,
    metacerberus_db,
    metacerberus_setup_params,
    metacerberus_params,
    pharokka_db,
    pharokka_params,
    dry_run,
    forceall,
    configfile,
    unlock,
    cores,
    jobs,
    latency_wait,
    rerun_incomplete,
    rerun_triggers,
    sdm,
    executor,
    quiet,
    snakemake_args,
):
    logging.info("Running module: viral-annotate")
    logging.info(f"outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")

    module_obj = ViralAnnotateModule()
    module_obj.name = "viral-annotate"
    # Set the attributes of the module object
    module_obj = setOptions(
        module_obj=module_obj,
        workdir=workdir,
        outdir=outdir,
        datadir=datadir,
        samplelist=samplelist,
        fasta=fasta,
        fastadir=fastadir,
        sample_name=sample_name,
        assembly_ids=assembly_ids,
        latest_run=latest_run,
        splits=splits,
        keep_intermediates=keep_intermediates,
        setup_database=setup_database,
        max_cores=max_cores,
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        custom_config=custom_config,
    )

    if eggnog_params:
        module_obj.eggNOG_params = eggnog_params
        module_obj.hasOptions = True
    if phavip_params:
        module_obj.PhaVIP_params = phavip_params
        module_obj.hasOptions = True
    if metacerberus_db:
        module_obj.metacerberus_db = metacerberus_db
        module_obj.hasOptions = True
    if metacerberus_setup_params:
        module_obj.metacerberus_setup_params = metacerberus_setup_params
        module_obj.hasOptions = True
    if metacerberus_params:
        module_obj.metacerberus_params = metacerberus_params
        module_obj.hasOptions = True
    if pharokka_db:
        module_obj.pharokka_db = pharokka_db
        module_obj.hasOptions = True
    if pharokka_params:
        module_obj.pharokka_params = pharokka_params
        module_obj.hasOptions = True
    snakemake_obj = SnakemakeFlags(
        dry_run,
        forceall,
        configfile,
        unlock,
        cores,
        jobs,
        latency_wait,
        rerun_incomplete,
        rerun_triggers,
        sdm,
        executor,
        quiet,
        snakemake_args,
    )

    vomix_actions_instance = vomix_actions()
    vomix_actions_instance.run_module("viral-annotate", module_obj, snakemake_obj)
    logging.info(END_MODULE_RUN_LOG)


# Prok Community Module


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
    help='Parameters for metaphlan function. See more at https://huttenhower.sph.harvard.edu/metaphlan/ || default: "--ignore_eukaryotes" [STR]',
)
@click.option(
    "--mpa-indexv",
    required=False,
    default=None,
    help='Database version for metaphlan to use. See more at https://huttenhower.sph.harvard.edu/metaphlan/ || default: "mpa_vOct22_CHOCOPhlAnSGB_202212" [STR]',
)
@snakemake_options
def run_prok_community(
    workdir,
    outdir,
    datadir,
    samplelist,
    custom_config,
    fasta,
    fastadir,
    sample_name,
    assembly_ids,
    latest_run,
    splits,
    keep_intermediates,
    setup_database,
    max_cores,
    ncbi_email,
    ncbi_api_key,
    mpa_params,
    mpa_indexv,
    dry_run,
    forceall,
    configfile,
    unlock,
    cores,
    jobs,
    latency_wait,
    rerun_incomplete,
    rerun_triggers,
    sdm,
    executor,
    quiet,
    snakemake_args,
):
    logging.info("Running module: prok-community")
    logging.info(f"outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")

    module_obj = ProkaryoticCommunityModule()
    module_obj.name = "prok-community"
    # Set the attributes of the module object
    module_obj = setOptions(
        module_obj=module_obj,
        workdir=workdir,
        outdir=outdir,
        datadir=datadir,
        samplelist=samplelist,
        fasta=fasta,
        fastadir=fastadir,
        sample_name=sample_name,
        assembly_ids=assembly_ids,
        latest_run=latest_run,
        splits=splits,
        keep_intermediates=keep_intermediates,
        setup_database=setup_database,
        max_cores=max_cores,
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        custom_config=custom_config,
    )

    if mpa_params:
        module_obj.mpa_params = mpa_params
        module_obj.hasOptions = True
    if mpa_indexv:
        module_obj.mpa_indexv = mpa_indexv
        module_obj.hasOptions = True
    snakemake_obj = SnakemakeFlags(
        dry_run,
        forceall,
        configfile,
        unlock,
        cores,
        jobs,
        latency_wait,
        rerun_incomplete,
        rerun_triggers,
        sdm,
        executor,
        quiet,
        snakemake_args,
    )

    vomix_actions_instance = vomix_actions()
    vomix_actions_instance.run_module("prok-community", module_obj, snakemake_obj)
    logging.info(END_MODULE_RUN_LOG)


# Prok Binning Module


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
    help="The directory path where the CheckM2 database is installed or will be downloaded. Defaults to the Snakemake base directory under workflow/databases.",
)
@click.option(
    "--GTDBTk-db",
    required=False,
    default=None,
    help='The directory path where the GTDB-Tk database is installed or will be downloaded. Defaults to the Snakemake base directory under workflow/databases. (default: "database/GTDB-Tk")',
)
@click.option(
    "--GTDBTk-db-version",
    required=False,
    default=232,
    help="The reference version of the GTDB-Tk database. Ensure that the database version corresponds with the local GTDB-Tk installation environment as detailed at https://ecogenomics.github.io/GTDBTk/installing/index.html#gtdb-tk-reference-data. (default: 232)",
)
@click.option(
    "--GTDBTk-identify-params",
    required=False,
    default=None,
    help='Additional parameter variables supplied directly to the GTDB-Tk identify command execution. (default: "")',
)
@click.option(
    "--GTDBTk-align-params",
    required=False,
    default=None,
    help='Additional parameter variables supplied directly to the GTDB-Tk align command execution. (default: "")',
)
@click.option(
    "--GTDBTk-classify-params",
    required=False,
    default=None,
    help='Additional parameter variables supplied directly to the GTDB-Tk classify command execution. (default: "")',
)
@click.option(
    "--VAMB-params",
    required=False,
    default=None,
    help='Additional parameter variables supplied directly to the VAMB command execution. (default: "")',
)
@click.option(
    "--binning-consensus",
    is_flag=True,
    flag_value=True,
    required=False,
    default=True,
    help="Flag enabling a consensus-based metagenomic binning protocol combining MetaBAT2, MaxBin2, and CONCOCT via DASTool. Disabling this runs GPU-accelerated VAMB clustering instead. (default: True)",
)
@click.option(
    "--strobealign-params",
    required=False,
    default=None,
    help='Additional alignment flags or scoring rules passed to the strobealign tool backend. (default: "")',
)
@click.option(
    "--MetaBAT2-params",
    required=False,
    default=None,
    help='Additional parameters for the MetaBAT2 tool. (default: "-m 1500")',
)
@click.option(
    "--MaxBin2-params",
    required=False,
    default=None,
    help='Additional parameters for the MaxBin2 tool. (default: "-min_contig_length 1500 -max_iteration 50 -prob_threshold 0.9")',
)
@click.option(
    "--CONCOCT-params",
    required=False,
    default=None,
    help='Additional parameters for the CONCOCT tool. (default: "")',
)
@click.option(
    "--jgi-summarize-params",
    required=False,
    default=None,
    help='Additional runtime parameters supplied to the jgi_summarize_bam_contig_depth depth processing command. (default: "--percentIdentity 97")',
)
@click.option(
    "--DASTool-params",
    required=False,
    default=None,
    help='Additional parameters for the DASTool tool. (default: "")',
)
@click.option(
    "--checkm2-params",
    required=False,
    default=None,
    help='Additional parameter fields provided directly to the CheckM2 bin validation pipeline. (default: "")',
)
@click.option(
    "--galah-params",
    required=False,
    default=None,
    help='Additional parameters for the Galah tool. (default: "--ani 95 --min-aligned-fraction 15 --fragment-length 3000")',
)
@snakemake_options
def run_prok_binning(
    workdir,
    outdir,
    datadir,
    samplelist,
    custom_config,
    fasta,
    fastadir,
    sample_name,
    assembly_ids,
    latest_run,
    splits,
    keep_intermediates,
    setup_database,
    max_cores,
    ncbi_email,
    ncbi_api_key,
    checkm2_db,
    GTDBTk_db,
    GTDBTk_db_version,
    GTDBTk_identify_params,
    GTDBTk_align_params,
    GTDBTk_classify_params,
    VAMB_params,
    binning_consensus,
    strobealign_params,
    MetaBAT2_params,
    MaxBin2_params,
    CONCOCT_params,
    jgi_summarize_params,
    DASTool_params,
    checkm2_params,
    galah_params,
    dry_run,
    forceall,
    configfile,
    unlock,
    cores,
    jobs,
    latency_wait,
    rerun_incomplete,
    rerun_triggers,
    sdm,
    executor,
    quiet,
    snakemake_args,
):
    logging.info("Running module: prok-binning")
    logging.info(f"outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")

    module_obj = ProkaryoticBinningModule()
    module_obj.name = "prok-binning"
    # Set the attributes of the module object
    module_obj = setOptions(
        module_obj=module_obj,
        workdir=workdir,
        outdir=outdir,
        datadir=datadir,
        samplelist=samplelist,
        fasta=fasta,
        fastadir=fastadir,
        sample_name=sample_name,
        assembly_ids=assembly_ids,
        latest_run=latest_run,
        splits=splits,
        keep_intermediates=keep_intermediates,
        setup_database=setup_database,
        max_cores=max_cores,
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        custom_config=custom_config,
    )

    if checkm2_db:
        module_obj.checkm2_db = checkm2_db
        module_obj.hasOptions = True
    if GTDBTk_db:
        module_obj.GTDBTk_db = GTDBTk_db
        module_obj.hasOptions = True
    if GTDBTk_db_version:
        module_obj.GTDBTk_db_version = GTDBTk_db_version
        module_obj.hasOptions = True
    if GTDBTk_identify_params:
        module_obj.GTDBTk_identify_params = GTDBTk_identify_params
        module_obj.hasOptions = True
    if GTDBTk_align_params:
        module_obj.GTDBTk_align_params = GTDBTk_align_params
        module_obj.hasOptions = True
    if GTDBTk_classify_params:
        module_obj.GTDBTk_classify_params = GTDBTk_classify_params
        module_obj.hasOptions = True
    if VAMB_params:
        module_obj.VAMB_params = VAMB_params
        module_obj.hasOptions = True
    if binning_consensus is not None:
        module_obj.binning_consensus = binning_consensus
        module_obj.hasOptions = True
    if strobealign_params:
        module_obj.strobealign_params = strobealign_params
        module_obj.hasOptions = True
    if MetaBAT2_params:
        module_obj.MetaBAT2_params = MetaBAT2_params
        module_obj.hasOptions = True
    if MaxBin2_params:
        module_obj.MaxBin2_params = MaxBin2_params
        module_obj.hasOptions = True
    if CONCOCT_params:
        module_obj.CONCOCT_params = CONCOCT_params
        module_obj.hasOptions = True
    if jgi_summarize_params:
        module_obj.jgi_summarize_params = jgi_summarize_params
        module_obj.hasOptions = True
    if DASTool_params:
        module_obj.DASTool_params = DASTool_params
        module_obj.hasOptions = True
    if checkm2_params:
        module_obj.checkm2_params = checkm2_params
        module_obj.hasOptions = True
    if galah_params:
        module_obj.galah_params = galah_params
        module_obj.hasOptions = True

    snakemake_obj = SnakemakeFlags(
        dry_run,
        forceall,
        configfile,
        unlock,
        cores,
        jobs,
        latency_wait,
        rerun_incomplete,
        rerun_triggers,
        sdm,
        executor,
        quiet,
        snakemake_args,
    )

    vomix_actions_instance = vomix_actions()
    vomix_actions_instance.run_module("prok-binning", module_obj, snakemake_obj)
    logging.info(END_MODULE_RUN_LOG)


# Prok Annotate Module


@cli.command(
    "prok-annotate",
    context_settings={"ignore_unknown_options": True},
    short_help="Run the Prokaryotic Annotate module",
)
@common_options
@snakemake_options
@click.option("--humann-params", required=False, default=None, help="[STR]")
def run_prok_annotate(
    workdir,
    outdir,
    datadir,
    samplelist,
    fasta,
    fastadir,
    sample_name,
    assembly_ids,
    latest_run,
    splits,
    keep_intermediates,
    setup_database,
    max_cores,
    ncbi_email,
    ncbi_api_key,
    custom_config,
    dry_run,
    forceall,
    configfile,
    unlock,
    cores,
    jobs,
    latency_wait,
    rerun_incomplete,
    rerun_triggers,
    sdm,
    executor,
    quiet,
    snakemake_args,
    humann_params,
):
    logging.info("Running module: prok-annotate")
    logging.info(f"outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")

    module_obj = ProkaryoticAnnotateModule()
    module_obj.name = "prok-annotate"
    # Set the attributes of the module object
    module_obj = setOptions(
        module_obj=module_obj,
        workdir=workdir,
        outdir=outdir,
        datadir=datadir,
        samplelist=samplelist,
        fasta=fasta,
        fastadir=fastadir,
        sample_name=sample_name,
        assembly_ids=assembly_ids,
        latest_run=latest_run,
        splits=splits,
        keep_intermediates=keep_intermediates,
        setup_database=setup_database,
        max_cores=max_cores,
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        custom_config=custom_config,
    )
    snakemake_obj = SnakemakeFlags(
        dry_run,
        forceall,
        configfile,
        unlock,
        cores,
        jobs,
        latency_wait,
        rerun_incomplete,
        rerun_triggers,
        sdm,
        executor,
        quiet,
        snakemake_args,
    )

    if humann_params:
        module_obj.humann_params = humann_params

    vomix_actions_instance = vomix_actions()
    vomix_actions_instance.run_module("prok-annotate", module_obj, snakemake_obj)
    logging.info(END_MODULE_RUN_LOG)


# Viral End to End Module


@cli.command(
    "viral-end-to-end",
    context_settings={"ignore_unknown_options": True},
    short_help="Run the Viral End-To-End module",
)
@common_options
@snakemake_options
def run_viral_end_to_end(
    workdir,
    outdir,
    datadir,
    samplelist,
    fasta,
    fastadir,
    sample_name,
    assembly_ids,
    latest_run,
    splits,
    keep_intermediates,
    setup_database,
    max_cores,
    ncbi_email,
    ncbi_api_key,
    custom_config,
    dry_run,
    forceall,
    configfile,
    unlock,
    cores,
    jobs,
    latency_wait,
    rerun_incomplete,
    rerun_triggers,
    sdm,
    executor,
    quiet,
    snakemake_args,
):
    logging.info("Running module: viral-end-to-end")
    logging.info(f"outdir: {outdir}, datadir: {datadir}, samplelist: {samplelist}")

    module_obj = ViralEndToEndModule()
    module_obj.name = "viral-end-to-end"

    # Set the attributes of the module object
    module_obj = setOptions(
        module_obj=module_obj,
        workdir=workdir,
        outdir=outdir,
        datadir=datadir,
        samplelist=samplelist,
        fasta=fasta,
        fastadir=fastadir,
        sample_name=sample_name,
        assembly_ids=assembly_ids,
        latest_run=latest_run,
        splits=splits,
        keep_intermediates=keep_intermediates,
        setup_database=setup_database,
        max_cores=max_cores,
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        custom_config=custom_config,
    )
    snakemake_obj = SnakemakeFlags(
        dry_run,
        forceall,
        configfile,
        unlock,
        cores,
        jobs,
        latency_wait,
        rerun_incomplete,
        rerun_triggers,
        sdm,
        executor,
        quiet,
        snakemake_args,
    )

    vomix_actions_instance = vomix_actions()
    vomix_actions_instance.run_module("viral-end-to-end", module_obj, snakemake_obj)
    logging.info(END_MODULE_RUN_LOG)


# Cluster Fast Module


@cli.command(
    "cluster-fast",
    context_settings={"ignore_unknown_options": True},
    short_help="Run the Cluster Fast module",
)
@common_options
@click.option(
    "--clustering-fast",
    required=False,
    default=None,
    help="Flag to run fast clustering using CheckV's MEGABLAST approach. If set to False, CD-HIT will be used. Proceed with caution as it can be extremely slow at large sequence numbers. || default: True [True or False]",
)
@click.option(
    "--cdhit-params",
    required=False,
    default=None,
    help='Additional parameters to pass on to CD-HIT if clustering-fast is set to False. Read more at https://github.com/weizhongli/cdhit/blob/master/doc/cdhit-user-guide.wiki || default: "-c 0.95 -aS 0.85 -d 400 -M 0 -n 5" [STR]',
)
@click.option(
    "--vOTU-ani",
    required=False,
    default=None,
    help="Minimum average nucleotide identity for fast clustering algorithm of viral contigs || default: 95 [INT]",
)
@click.option(
    "--vOTU-targetcov",
    required=False,
    default=None,
    help="Minimum target coverage for fast clustering algorithm of viral contigs || default: 85 [NUM]",
)
@click.option(
    "--vOTU-querycov",
    required=False,
    default=None,
    help="Minimum query coverage for fast clustering algorithm of viral contigs || default: 0 [NUM]",
)
@snakemake_options
def run_cluster_fast(
    workdir,
    outdir,
    datadir,
    samplelist,
    custom_config,
    fasta,
    fastadir,
    sample_name,
    assembly_ids,
    latest_run,
    splits,
    keep_intermediates,
    setup_database,
    max_cores,
    ncbi_email,
    ncbi_api_key,
    clustering_fast,
    cdhit_params,
    votu_ani,
    votu_targetcov,
    votu_querycov,
    dry_run,
    forceall,
    configfile,
    unlock,
    cores,
    jobs,
    latency_wait,
    rerun_incomplete,
    rerun_triggers,
    sdm,
    executor,
    quiet,
    snakemake_args,
):
    logging.info("Running module: cluster-fast")
    logging.info(f"fasta: {fasta}, outdir: {outdir}")

    module_obj = ClusterFastModule()
    module_obj.name = "cluster-fast"
    # Set the attributes of the module object
    module_obj = setOptions(
        module_obj=module_obj,
        workdir=workdir,
        outdir=outdir,
        datadir=datadir,
        samplelist=samplelist,
        fasta=fasta,
        fastadir=fastadir,
        sample_name=sample_name,
        assembly_ids=assembly_ids,
        latest_run=latest_run,
        splits=splits,
        keep_intermediates=keep_intermediates,
        setup_database=setup_database,
        max_cores=max_cores,
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        custom_config=custom_config,
    )

    if clustering_fast:
        module_obj.clustering_fast = clustering_fast
        module_obj.hasOptions = True
    if cdhit_params:
        module_obj.cdhit_params = cdhit_params
        module_obj.hasOptions = True
    if votu_ani:
        module_obj.vOTU_ani = votu_ani
        module_obj.hasOptions = True
    if votu_targetcov:
        module_obj.vOTU_targetcov = votu_targetcov
        module_obj.hasOptions = True
    if votu_querycov:
        module_obj.vOTU_querycov = votu_querycov
        module_obj.hasOptions = True
    snakemake_obj = SnakemakeFlags(
        dry_run,
        forceall,
        configfile,
        unlock,
        cores,
        jobs,
        latency_wait,
        rerun_incomplete,
        rerun_triggers,
        sdm,
        executor,
        quiet,
        snakemake_args,
    )

    vomix_actions_instance = vomix_actions()
    vomix_actions_instance.run_module("cluster-fast", module_obj, snakemake_obj)
    logging.info(END_MODULE_RUN_LOG)


# CheckV-PyHMMER Module


@cli.command(
    "checkv-pyhmmer",
    context_settings={"ignore_unknown_options": True},
    short_help="Run the CheckV PyHMMER module",
)
@common_options
@click.option(
    "--checkv-original",
    required=False,
    default=None,
    help="Flag to use CheckV original instead of the much faster version in vOMIX-MEGA, CheckV-PyHMMER. || default: False [True or False]",
)
@click.option(
    "--checkv-params",
    required=False,
    default=None,
    help='Additional parameters to pass on to CheckV. Read more at https://bitbucket.org/berkeleylab/CheckV/src || default: "" [STR]',
)
@click.option(
    "--checkv-database",
    required=False,
    default=None,
    help='Path to CheckV\'s database || default: "workflow/database/checkv" [STR]',
)
@snakemake_options
def run_checkv_pyhmmer(
    workdir,
    outdir,
    datadir,
    samplelist,
    custom_config,
    fasta,
    fastadir,
    sample_name,
    assembly_ids,
    latest_run,
    splits,
    keep_intermediates,
    setup_database,
    max_cores,
    ncbi_email,
    ncbi_api_key,
    checkv_original,
    checkv_params,
    checkv_database,
    dry_run,
    forceall,
    configfile,
    unlock,
    cores,
    jobs,
    latency_wait,
    rerun_incomplete,
    rerun_triggers,
    sdm,
    executor,
    quiet,
    snakemake_args,
):
    logging.info("Running module: checkv-pyhmmer")
    logging.info(f"fasta: {fasta}, outdir: {outdir}")

    module_obj = CheckVPyHMMERModule()
    module_obj.name = "checkv-pyhmmer"
    # Set the attributes of the module object
    module_obj = setOptions(
        module_obj=module_obj,
        workdir=workdir,
        outdir=outdir,
        datadir=datadir,
        samplelist=samplelist,
        fasta=fasta,
        fastadir=fastadir,
        sample_name=sample_name,
        assembly_ids=assembly_ids,
        latest_run=latest_run,
        splits=splits,
        keep_intermediates=keep_intermediates,
        setup_database=setup_database,
        max_cores=max_cores,
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        custom_config=custom_config,
    )

    if checkv_original:
        module_obj.checkv_original = checkv_original
        module_obj.hasOptions = True
    if checkv_params:
        module_obj.checkv_params = checkv_params
        module_obj.hasOptions = True
    if checkv_database:
        module_obj.checkv_database = checkv_database
        module_obj.hasOptions = True
    snakemake_obj = SnakemakeFlags(
        dry_run,
        forceall,
        configfile,
        unlock,
        cores,
        jobs,
        latency_wait,
        rerun_incomplete,
        rerun_triggers,
        sdm,
        executor,
        quiet,
        snakemake_args,
    )

    vomix_actions_instance = vomix_actions()
    vomix_actions_instance.run_module("checkv-pyhmmer", module_obj, snakemake_obj)
    logging.info(END_MODULE_RUN_LOG)


# Setup-Database Module


@cli.command(
    "setup-database",
    context_settings={"ignore_unknown_options": True},
    short_help="Run the Setup Database module",
)
@common_options
@click.option("--hostile-index-db", required=False, default=None, help="")
@click.option(
    "--PhaBox2-db",
    required=False,
    default=None,
    help='Path to PhaBox2 database for download || default: "workflow/database/phabox_db_v2" [STR]',
)
@click.option(
    "--phabox2-db-name",
    required=False,
    default=None,
    help='The designated database name or identifier file package required for the PhaBox2 classification tool execution. (default: "phabox_db_v2")',
)
@click.option(
    "--phabox2-db-baselink",
    required=False,
    default=None,
    help='The primary remote server base link URL used to fetch and download resource updates for the PhaBox2 database. (default: "https://github.com/KennthShang/PhaBOX/releases/download/v2")',
)
@click.option(
    "--genomad-db",
    required=False,
    default=None,
    help='Path to geNomad database for download || default: "workflow/database/genomad" [STR]',
)
@click.option(
    "--virsorter2-db",
    required=False,
    default=None,
    help='The directory path where the VirSorter2 database is installed or will be downloaded. Defaults to the Snakemake base directory under workflow/databases. (default: "database/virsorter2")',
)
@click.option(
    "--vibrant-db",
    required=False,
    default=None,
    help='The directory path where the VIBRANT database is installed or will be downloaded. Defaults to the Snakemake base directory under workflow/databases. (default: "database/virsorter2")',
)
@click.option(
    "--checkv-db",
    required=False,
    default=None,
    help='Path to CheckV database for download || default: "workflow/database/phabox_db_v2" [STR]',
)
@click.option(
    "--eggNOG-db",
    required=False,
    default=None,
    help='Path to eggNOG v2 database for download || default: "workflow/database/eggNOGv2" [STR]',
)
@click.option(
    "--eggNOG-db-params",
    required=False,
    default=None,
    help='Parameters for downloading eggNOG v2 database || default: "" [STR]',
)
@click.option(
    "--checkm2-db",
    required=False,
    default=None,
    help="The directory path where the CheckM2 database is installed or will be downloaded. Defaults to the Snakemake base directory under workflow/databases.",
)
@click.option(
    "--GTDBTk-db",
    required=False,
    default=None,
    help='The directory path where the GTDB-Tk database is installed or will be downloaded. Defaults to the Snakemake base directory under workflow/databases. (default: "database/GTDB-Tk")',
)
@click.option(
    "--GTDBTk-db-version",
    required=False,
    default=232,
    help="The reference version of the GTDB-Tk database. Ensure that the database version corresponds with the local GTDB-Tk installation environment as detailed at https://ecogenomics.github.io/GTDBTk/installing/index.html#gtdb-tk-reference-data. (default: 232)",
)
@click.option(
    "--iphop-db",
    required=False,
    default=None,
    help='Path to iPHoP database for download || default: "workflow/database/iphop/Aug_2023_pub_rw" [STR]',
)
@click.option(
    "--iphop-db-version",
    required=False,
    default=None,
    help='The version identifier for the iPHoP database. Verify database compatibility mapping at https://bitbucket.org/srouxjgi/iphop/src/main/#markdown-header-host-databases-and-versions. (default: "iPHoP_db_Aug23_rw")',
)
@click.option(
    "--iphop-db-basename",
    required=False,
    default=None,
    help='The primary base name of the iPHoP database. Verify compatibility criteria at https://bitbucket.org/srouxjgi/iphop/src/main/#markdown-header-host-databases-and-versions. (default: "Aug_2023_pub_rw")',
)
@click.option(
    "--humann-db",
    required=False,
    default=None,
    help='Path to HUMAnN3 databases for download || default: "workflow/database/humann" [STR]',
)
@click.option(
    "--metacerberus-db",
    required=False,
    default=None,
    help='The directory path where the MetaCerberus database is installed or will be downloaded. Defaults to the Snakemake base directory under workflow/databases. (default: "database/metacerberus")',
)
@click.option(
    "--metacerberus-setup-params",
    required=False,
    default=None,
    help='Operational configurations supplied to initialize build or index the MetaCerberus database environment. (default: "")',
)
@click.option(
    "--pharokka-db",
    required=False,
    default=None,
    help='The directory path where the pharokka database is installed or will be downloaded. Defaults to the Snakemake base directory under workflow/databases. (default: "database/pharokka")',
)
@snakemake_options
def run_setup_database(
    workdir,
    outdir,
    datadir,
    samplelist,
    custom_config,
    fasta,
    fastadir,
    sample_name,
    assembly_ids,
    latest_run,
    splits,
    keep_intermediates,
    setup_database,
    max_cores,
    ncbi_email,
    ncbi_api_key,
    phabox2_db,
    phabox2_db_name,
    phabox2_db_baselink,
    genomad_db,
    virsorter2_db,
    vibrant_db,
    checkv_db,
    eggnog_db,
    eggnog_db_params,
    checkm2_db,
    gtdbtk_db,
    gtdbtk_db_version,
    iphop_db,
    iphop_db_version,
    iphop_db_basename,
    humann_db,
    metacerberus_db,
    metacerberus_setup_params,
    pharokka_db,
    dry_run,
    forceall,
    configfile,
    unlock,
    cores,
    jobs,
    latency_wait,
    rerun_incomplete,
    rerun_triggers,
    sdm,
    executor,
    quiet,
    snakemake_args,
    hostile_index_db,
):
    logging.info("Running module: setup-database")
    logging.info(f"fasta: {fasta}, outdir: {outdir}")

    module_obj = SetupDatabaseModule()
    module_obj.name = "setup-database"
    # Set the attributes of the module object
    module_obj = setOptions(
        module_obj=module_obj,
        workdir=workdir,
        outdir=outdir,
        datadir=datadir,
        samplelist=samplelist,
        fasta=fasta,
        fastadir=fastadir,
        sample_name=sample_name,
        assembly_ids=assembly_ids,
        latest_run=latest_run,
        splits=splits,
        keep_intermediates=keep_intermediates,
        setup_database=setup_database,
        max_cores=max_cores,
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        custom_config=custom_config,
    )

    if phabox2_db:
        module_obj.PhaBox2_db = phabox2_db
        module_obj.hasOptions = True
    if phabox2_db_name:
        module_obj.phabox2_db_name = phabox2_db_name
        module_obj.hasOptions = True
    if phabox2_db_baselink:
        module_obj.phabox2_db_baselink = phabox2_db_baselink
        module_obj.hasOptions = True
    if genomad_db:
        module_obj.genomad_db = genomad_db
        module_obj.hasOptions = True
    if checkv_db:
        module_obj.checkv_db = checkv_db
        module_obj.hasOptions = True
    if eggnog_db:
        module_obj.eggNOG_db = eggnog_db
        module_obj.hasOptions = True
    if eggnog_db_params:
        module_obj.eggNOG_db_params = eggnog_db_params
        module_obj.hasOptions = True
    if virsorter2_db:
        module_obj.virsorter2_db = virsorter2_db
        module_obj.hasOptions = True
    if iphop_db:
        module_obj.iphop_db = iphop_db
        module_obj.hasOptions = True
    if humann_db:
        module_obj.humann_db = humann_db
        module_obj.hasOptions = True
    if hostile_index_db:
        module_obj.hostile_index_db = hostile_index_db
        module_obj.hasOptions = True
    if metacerberus_db:
        module_obj.metacerberus_db = metacerberus_db
        module_obj.hasOptions = True
    if metacerberus_setup_params:
        module_obj.metacerberus_setup_params = metacerberus_setup_params
        module_obj.hasOptions = True
    if pharokka_db:
        module_obj.pharokka_db = pharokka_db
        module_obj.hasOptions = True
    if checkm2_db:
        module_obj.checkm2_db = checkm2_db
        module_obj.hasOptions = True
    if gtdbtk_db:
        module_obj.GTDBTk_db = gtdbtk_db
        module_obj.hasOptions = True
    if gtdbtk_db_version:
        module_obj.GTDBTk_db_version = gtdbtk_db_version
        module_obj.hasOptions = True
    if iphop_db_version:
        module_obj.iphop_db_version = iphop_db_version
        module_obj.hasOptions = True
    if iphop_db_basename:
        module_obj.iphop_db_basename = iphop_db_basename
        module_obj.hasOptions = True
    if vibrant_db:
        module_obj.vibrant_db = vibrant_db
        module_obj.hasOptions = True

    snakemake_obj = SnakemakeFlags(
        dry_run,
        forceall,
        configfile,
        unlock,
        cores,
        jobs,
        latency_wait,
        rerun_incomplete,
        rerun_triggers,
        sdm,
        executor,
        quiet,
        snakemake_args,
    )

    vomix_actions_instance = vomix_actions()
    vomix_actions_instance.run_module("setup-database", module_obj, snakemake_obj)
    logging.info(END_MODULE_RUN_LOG)
