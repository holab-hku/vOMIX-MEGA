import sys
from rich.console import Console
from rich.panel import Panel

console = Console()

logdir = relpath("community/metaphlan/logs")
tmpd = relpath("community/metaphlan/tmp")
benchmarks = relpath("community/metaphlan/benchmarks")

email = config["NCBI-email"]
api_key = config["NCBI-API-key"]
nowstr = config["latest-run"]
outdir = config["outdir"]
datadir = config["datadir"]

parse_quiet = config.get("module") == "viral-end-to-end"
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

long_read_samples = [
    s for s, info in samples.items()
    if info.get("read_type") in ["pacbio", "nanopore"]
]
if long_read_samples:
    console.print(
        Panel.fit(
            f"[bold red]ERROR:[/] MetaPhlAn does not support long reads.\n"
            f"Found long‑read samples: {long_read_samples}\n"
            f"Please remove these samples from the sample list for this module.",
            title="Long‑Read Not Supported",
            border_style="red"
        )
    )
    sys.exit(1)

os.makedirs(logdir, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)

# ----------------------------------------------------------------------
# MASTER RULE
# ----------------------------------------------------------------------
rule done_log:
    name: "prok-community.smk done. Removing tmp files"
    localrule: True
    input:
        expand(relpath("community/metaphlan/samples/{sample_id}/{sample_id}.txt"), sample_id=samples.keys()),
        relpath("community/metaphlan/output/metaphlan_out.txt"),
        expand(relpath("community/metaphlan/output/metaphlan_out_{level}.txt"),
               level=['phylum', 'class', 'order', 'family', 'genus', 'species', 'SGB'])
    output:
        os.path.join(logdir, "done.log")
    shell:
        """
        touch {output}
        """

# ----------------------------------------------------------------------
# METAPHLAN RUN (supports paired and single-end)
# ----------------------------------------------------------------------
rule metaphlan:
    name: "prok-community.smk MetaPhlAn with Bowtie2"
    input:
        R1=relpath("preprocess/samples/{sample_id}/output/{sample_id}_R1_cut.trim.filt.fastq.gz"),
        R2=relpath("preprocess/samples/{sample_id}/output/{sample_id}_R2_cut.trim.filt.fastq.gz"),
    output:
        rpkm=relpath("community/metaphlan/samples/{sample_id}/{sample_id}.txt"),
        sam=relpath("community/metaphlan/samples/{sample_id}/{sample_id}_bowtie_out.txt")
    params:
        outdir=relpath("community/metaphlan/samples/{sample_id}"),
        bowtiedir=relpath("community/metaphlan/samples/{sample_id}/bowtie"),
        db=os.path.join("workflow/database/metaphlan"),
        index_v=config["mpa-indexv"],
        parameters=config["mpa-params"],
        tmpdir=os.path.join(tmpd, "metaphlan/{sample_id}"),
        read_type=lambda wildcards: samples[wildcards.sample_id]["read_type"],
    log: os.path.join(logdir, "metaphlan_{sample_id}.log")
    benchmark: os.path.join(benchmarks, "metaphlan_{sample_id}.log")
    threads: 16
    resources:
        mem_mb=lambda wildcards, input, attempt: attempt * 20 * 10**3
    conda: "../envs/metaphlan.yml"
    shell:
        """
        rm -rf {params.tmpdir} {params.outdir}
        mkdir -p {params.tmpdir} {params.outdir}

        # Build input string: for paired-end, use both; for single, only R1
        if [ "{params.read_type}" == "single" ]; then
            INPUT="{input.R1}"
        else
            INPUT="{input.R1},{input.R2}"
        fi

        metaphlan $INPUT \
            --bowtie2db {params.db} \
            --bowtie2out {params.tmpdir}/tmp.sam \
            --input_type fastq \
            --index {params.index_v} \
            --sample_id {wildcards.sample_id} \
            -t rel_ab_w_read_stats \
            -o {params.tmpdir}/tmp.txt \
            --tmp_dir {params.tmpdir} \
            --offline {params.parameters} &> {log}

        mv {params.tmpdir}/tmp.sam {output.sam}
        mv {params.tmpdir}/tmp.txt {output.rpkm}
        """

# ----------------------------------------------------------------------
# MERGE METAPHLAN OUTPUTS
# ----------------------------------------------------------------------
rule metphlan_merge:
    name: "prok-community.smk merging MetaPhlAn outputs"
    localrule: True
    input:
        expand(relpath("community/metaphlan/samples/{sample_id}/{sample_id}.txt"), sample_id=samples.keys())
    output:
        relpath("community/metaphlan/output/metaphlan_out.txt")
    params:
        outdir=relpath("community/metaphlan/output"),
        tmpdir=os.path.join(tmpd, "metaphlan")
    log: os.path.join(logdir, "merge_metaphlan.log")
    threads: 1
    conda: "../envs/metaphlan.yml"
    shell:
        """
        rm -rf {params.tmpdir} {params.outdir}
        mkdir -p {params.tmpdir} {params.outdir}

        merge_metaphlan_tables.py {input} > {params.tmpdir}/tmp.txt

        mv {params.tmpdir}/tmp.txt {output}
        """

# ----------------------------------------------------------------------
# SEPARATE BY TAXONOMY LEVEL
# ----------------------------------------------------------------------
rule metaphlan_separate:
    name: "prok-community.smk separate MetaPhlAn outputs on taxonomy level"
    localrule: True
    input:
        relpath("community/metaphlan/output/metaphlan_out.txt")
    output:
        plev=relpath("community/metaphlan/output/metaphlan_out_phylum.txt"),
        clev=relpath("community/metaphlan/output/metaphlan_out_class.txt"),
        olev=relpath("community/metaphlan/output/metaphlan_out_order.txt"),
        flev=relpath("community/metaphlan/output/metaphlan_out_family.txt"),
        glev=relpath("community/metaphlan/output/metaphlan_out_genus.txt"),
        slev=relpath("community/metaphlan/output/metaphlan_out_species.txt"),
        tlev=relpath("community/metaphlan/output/metaphlan_out_SGB.txt")
    log: os.path.join(logdir, "process_metaphlan.log")
    threads: 1
    shell:
        """
        head {input} -n 2 > {output.tlev}
        head {input} -n 2 > {output.slev}
        head {input} -n 2 > {output.glev}
        head {input} -n 2 > {output.flev}
        head {input} -n 2 > {output.olev}
        head {input} -n 2 > {output.clev}
        head {input} -n 2 > {output.plev}

        grep -E "t__" {input} >> {output.tlev}
        grep -E "s__" {input} | grep -v "t__" >> {output.slev}
        grep -E "g__" {input} | grep -v "s__" >> {output.glev}
        grep -E "f__" {input} | grep -v "g__" >> {output.flev}
        grep -E "o__" {input} | grep -v "f__" >> {output.olev}
        grep -E "c__" {input} | grep -v "o_" >> {output.clev}
        grep -E "p__" {input} | grep -v "c__" >> {output.plev}
        """