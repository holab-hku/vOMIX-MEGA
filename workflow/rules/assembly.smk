import sys
import os
from rich.console import Console
from rich.panel import Panel

console = Console()

# ----------------------------------------------------------------------
# Configuration & setup
# ----------------------------------------------------------------------
short_assembler = config.get("short-read-assembler", "megahit")
logdir = relpath(os.path.join("assembly", "logs"))
tmpd = relpath(os.path.join("assembly", "tmp"))
benchmarks = relpath(os.path.join("assembly", "benchmarks"))

os.makedirs(logdir, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)

email = config["NCBI-email"]
api_key = config["NCBI-API-key"]
nowstr = config["latest-run"]
outdir = config["outdir"]
datadir = config["datadir"]

parse_quiet = config.get("module") in ["viral-end-to-end", "run-all"]
parse_verbose = config.get("verbose", False)
samples, assemblies = parse_sample_list(
    config["samplelist"],
    datadir,
    outdir,
    email,
    api_key,
    nowstr,
    quiet=parse_quiet,
    verbose=parse_verbose
)

# ----------------------------------------------------------------------
# Mapping: read_type -> assembler name
# ----------------------------------------------------------------------
assembler_map = {
    "paired": short_assembler,
    "single": short_assembler,
    "pacbio": "metamdbg",
    "nanopore": "nanomdbg",
}

# Check co‑assembly with SPAdes
if short_assembler == "spades" and len(assemblies.keys()) != len(samples.keys()):
    console.print(
        Panel.fit(
            "[bold red]ERROR:[/] SPAdes does not support co‑assembly.\n"
            "Use assembler='megahit' or split assemblies.",
            title="SPAdes Co‑assembly Error",
            border_style="red"
        )
    )
    sys.exit(1)

# ----------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------
def get_assembler(assembly_id):
    """Return the assembler name for a given assembly."""
    rt = assemblies[assembly_id]["read_type"]
    return assembler_map.get(rt, short_assembler)

def get_assembler_specific_path(assembly_id):
    """Return the assembler-specific final contigs FASTA path."""
    assembler = get_assembler(assembly_id)
    return relpath(os.path.join("assembly", assembler, "samples", assembly_id, "output", "final.contigs.fa"))

def get_common_output_dir(assembly_id):
    """Return the common output directory (final copy location)."""
    return relpath(os.path.join("assembly", "samples", assembly_id, "output"))

def get_common_fasta(assembly_id):
    """Return the final contigs FASTA path in the common directory."""
    return os.path.join(get_common_output_dir(assembly_id), "final.contigs.fa")

def get_sample_ids(assembly_id):
    """Return a list of sample IDs belonging to the assembly."""
    ids = assemblies[assembly_id]["sample_id"]
    if isinstance(ids, str):
        return [ids]
    return list(ids)

# ----------------------------------------------------------------------
# MASTER RULE
# ----------------------------------------------------------------------
rule done:
    name: "assembly.smk Done. removing tmp files"
    localrule: True
    input:
        expand(
            relpath("assembly/samples/{assembly_id}/output/final.contigs.fa"),
            assembly_id=assemblies.keys()
        ),
        expand(
            relpath(os.path.join("assembly", "reports", "{summary_type}.tsv")),
            summary_type=["assemblystats", "assembly_size_dist"],
        ),
    output:
        os.path.join(logdir, "done.log")
    shell:
        """
        touch {output}
        """

# ----------------------------------------------------------------------
# ASSEMBLY RULES
# ----------------------------------------------------------------------

# ---- 1. MEGAHIT (short reads, paired or single) ----
rule megahit:
    name: "assembly.smk MEGAHIT assembly"
    input:
        R1s=lambda wildcards: expand(
            relpath("preprocess/samples/{sample_id}/output/{sample_id}_R1_cut.trim.filt.fastq.gz"),
            sample_id=get_sample_ids(wildcards.assembly_id)
        ),
        R2s=lambda wildcards: expand(
            relpath("preprocess/samples/{sample_id}/output/{sample_id}_R2_cut.trim.filt.fastq.gz"),
            sample_id=get_sample_ids(wildcards.assembly_id)
        ),
    output:
        fasta=relpath("assembly/megahit/samples/{assembly_id}/output/final.contigs.fa"),
    params:
        parameters=config["megahit-params"],
        minlen=config["megahit-min-len"],
        outdir=relpath("assembly/megahit/samples/{assembly_id}/output"),
        tmpdir=os.path.join(tmpd, "megahit"),
        read_type=lambda wildcards: assemblies[wildcards.assembly_id]["read_type"],
    log: os.path.join(logdir, "megahit_{assembly_id}.log")
    benchmark: os.path.join(benchmarks, "megahit_{assembly_id}.log")
    conda: "../envs/megahit.yml"
    threads: 24
    resources:
        mem_mb=lambda wildcards, attempt, threads, input: 8000 * attempt
    shell:
        """
        rm -rf {params.tmpdir}/{wildcards.assembly_id} {params.outdir}/*
        mkdir -p {params.outdir} {params.tmpdir}

        # Determine if paired or single
        if [ "{params.read_type}" == "paired" ]; then
            megahit \
                -1 $(echo "{input.R1s}" | tr ' ' ',') \
                -2 $(echo "{input.R2s}" | tr ' ' ',') \
                --min-contig-len {params.minlen} \
                -o {params.tmpdir}/{wildcards.assembly_id} \
                -t {threads} \
                {params.parameters} &> {log}
        else
            # single-end
            megahit \
                -r $(echo "{input.R1s}" | tr ' ' ',') \
                --min-contig-len {params.minlen} \
                -o {params.tmpdir}/{wildcards.assembly_id} \
                -t {threads} \
                {params.parameters} &> {log}
        fi

        mv {params.tmpdir}/{wildcards.assembly_id}/final.contigs.fa {output.fasta}
        mv {params.tmpdir}/{wildcards.assembly_id}/* {params.outdir}
        """

# ---- 2. SPAdes (short reads, paired or single) ----
rule spades:
    name: "assembly.smk SPAdes assembly"
    input:
        R1s=lambda wildcards: expand(
            relpath("preprocess/samples/{sample_id}/output/{sample_id}_R1_cut.trim.filt.fastq.gz"),
            sample_id=get_sample_ids(wildcards.assembly_id)
        ),
        R2s=lambda wildcards: expand(
            relpath("preprocess/samples/{sample_id}/output/{sample_id}_R2_cut.trim.filt.fastq.gz"),
            sample_id=get_sample_ids(wildcards.assembly_id)
        ),
    output:
        fasta=relpath("assembly/spades/samples/{assembly_id}/output/final.contigs.fa"),
    params:
        parameters=config["spades-params"],
        memory=config["spades-memory"],
        outdir=relpath("assembly/spades/samples/{assembly_id}/output"),
        tmpdir=os.path.join(tmpd, "spades"),
        read_type=lambda wildcards: assemblies[wildcards.assembly_id]["read_type"],
    log: os.path.join(logdir, "spades_{assembly_id}.log")
    benchmark: os.path.join(benchmarks, "spades_{assembly_id}.log")
    conda: "../envs/spades.yml"
    threads: 24
    resources:
        mem_mb=config["spades-memory"] * 1024
    shell:
        """
        rm -rf {params.tmpdir}/{wildcards.assembly_id} {params.outdir}/*
        mkdir -p {params.outdir} {params.tmpdir}

        # SPAdes supports both paired and single with -1/-2 or -s
        if [ "{params.read_type}" == "paired" ]; then
            spades.py \
                -1 $(echo "{input.R1s}" | tr ' ' ',') \
                -2 $(echo "{input.R2s}" | tr ' ' ',') \
                -o {params.tmpdir}/{wildcards.assembly_id} \
                -m {params.memory} \
                -t {threads} \
                {params.parameters} &> {log}
        else
            spades.py \
                -s $(echo "{input.R1s}" | tr ' ' ',') \
                -o {params.tmpdir}/{wildcards.assembly_id} \
                -m {params.memory} \
                -t {threads} \
                {params.parameters} &> {log}
        fi

        mv {params.tmpdir}/{wildcards.assembly_id}/scaffolds.fasta {output.fasta}
        mv {params.tmpdir}/{wildcards.assembly_id}/* {params.outdir}
        """

# ---- 3. metaMDBG (PacBio) ----
rule metamdbg:
    name: "assembly.smk metaMDBG assembly (PacBio)"
    input:
        reads=lambda wildcards: expand(
            relpath("preprocess/samples/{sample_id}/output/{sample_id}_R1_cut.trim.filt.fastq.gz"),
            sample_id=get_sample_ids(wildcards.assembly_id)
        ),
    output:
        fasta=relpath("assembly/metamdbg/samples/{assembly_id}/output/final.contigs.fa"),
    params:
        parameters=config.get("metamdbg-params", ""),
        outdir=relpath("assembly/metamdbg/samples/{assembly_id}/output"),
        tmpdir=os.path.join(tmpd, "metamdbg"),
    log: os.path.join(logdir, "metamdbg_{assembly_id}.log")
    benchmark: os.path.join(benchmarks, "metamdbg_{assembly_id}.log")
    conda: "../envs/metamdbg.yml"
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt, threads, input: 16000 * attempt
    shell:
        """
        rm -rf {params.tmpdir} {params.outdir}/*
        mkdir -p {params.outdir} {params.tmpdir}

        # metaMDBG expects a single input file (or list). We'll assume all reads are in one file.
        # If multiple, we can concatenate.
        if [ $(echo "{input.reads}" | wc -w) -gt 1 ]; then
            cat {input.reads} > {params.tmpdir}/combined.fastq.gz
            reads={params.tmpdir}/combined.fastq.gz
        else
            reads={input.reads}
        fi

        metamdbg assemble \
            -i $reads \
            -o {params.tmpdir}/assembly \
            -t {threads} \
            {params.parameters} &> {log}

        mv {params.tmpdir}/assembly/contigs.fasta {output.fasta}
        mv {params.tmpdir}/assembly/* {params.outdir}
        """

# ---- 4. nanoMDBG (ONT) ----
rule nanomdbg:
    name: "assembly.smk nanoMDBG assembly (Nanopore)"
    input:
        reads=lambda wildcards: expand(
            relpath("preprocess/samples/{sample_id}/output/{sample_id}_R1_cut.trim.filt.fastq.gz"),
            sample_id=get_sample_ids(wildcards.assembly_id)
        ),
    output:
        fasta=relpath("assembly/nanomdbg/samples/{assembly_id}/output/final.contigs.fa"),
    params:
        parameters=config.get("nanomdbg-params", ""),
        outdir=relpath("assembly/nanomdbg/samples/{assembly_id}/output"),
        tmpdir=os.path.join(tmpd, "nanomdbg"),
    log: os.path.join(logdir, "nanomdbg_{assembly_id}.log")
    benchmark: os.path.join(benchmarks, "nanomdbg_{assembly_id}.log")
    conda: "../envs/nanomdbg.yml"
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt, threads, input: 16000 * attempt
    shell:
        """
        rm -rf {params.tmpdir} {params.outdir}/*
        mkdir -p {params.outdir} {params.tmpdir}

        if [ $(echo "{input.reads}" | wc -w) -gt 1 ]; then
            cat {input.reads} > {params.tmpdir}/combined.fastq.gz
            reads={params.tmpdir}/combined.fastq.gz
        else
            reads={input.reads}
        fi

        nanomdbg assemble \
            -i $reads \
            -o {params.tmpdir}/assembly \
            -t {threads} \
            {params.parameters} &> {log}

        mv {params.tmpdir}/assembly/contigs.fasta {output.fasta}
        mv {params.tmpdir}/assembly/* {params.outdir}
        """

# ---- 5. Finalize: copy correct assembler output to common location ----
rule finalize_assembly:
    name: "assembly.smk finalize assembly output"
    localrule: True
    input:
        lambda wildcards: get_assembler_specific_path(wildcards.assembly_id)
    output:
        get_common_fasta("{assembly_id}")
    params:
        outdir = lambda wildcards: get_common_output_dir(wildcards.assembly_id)
    shell:
        """
        mkdir -p {params.outdir}
        cp {input} {output}
        """

# ----------------------------------------------------------------------
# ASSEMBLY STATS (aggregate reports)
# ----------------------------------------------------------------------
rule assembly_stats:
    name: "assembly.smk aggregate assembly statistics"
    localrule: False
    input:
        expand(
            relpath("assembly/samples/{assembly_id}/output/final.contigs.fa"),
            assembly_id=assemblies.keys()
        )
    output:
        stats=relpath(os.path.join("assembly", "reports", "assemblystats.tsv")),
        sizedist=relpath(os.path.join("assembly", "reports", "assembly_size_dist.tsv"))
    params:
        script="workflow/scripts/assembly_stats.py",
        outdir=relpath("assembly/reports"),
        tmpdir=os.path.join(tmpd, "reports")
    log: os.path.join(logdir, "stats.log")
    conda: "../envs/seqkit-biopython.yml"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt, threads, input: 8000 * attempt
    shell:
        """
        rm -rf {params.tmpdir} {params.outdir}/*
        mkdir -p {params.tmpdir} {params.outdir}

        python {params.script} -i {input} --size-dist-file {params.tmpdir}/tmp1.tsv > {params.tmpdir}/tmp2.tsv 2> {log}

        mv {params.tmpdir}/tmp1.tsv {output.sizedist}
        mv {params.tmpdir}/tmp2.tsv {output.stats}
        """