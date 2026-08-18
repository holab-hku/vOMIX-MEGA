# ----------------------------------------------------------------------
# Configuration & setup
# ----------------------------------------------------------------------
console = vomix_utils.console
Panel = vomix_utils.Panel

vomix_module = vomix_utils.current_module
logdir = vomix_module.logdir
benchmarks = vomix_module.benchmarks
tmpd = vomix_module.tmpd
samples = vomix_module.samples
assemblies = vomix_module.assemblies
fastap = vomix_module.fastap
sample_id = vomix_module.sample_id
assembly_ids = vomix_module.assembly_ids


# ----------------------------------------------------------------------
# Give error that longread is not supported
# ----------------------------------------------------------------------
long_read_samples = [
    s for s, info in samples.items()
    if info.get("read_type") in ["pacbio", "nanopore"]
]
if long_read_samples:
    console.print(
        Panel.fit(
            f"[bold red]ERROR:[/] Prokaryotic annotation (HUMAnN3) does not support long reads.\n"
            f"Found long‑read samples: {long_read_samples}\n"
            f"Please remove these samples from the sample list for this module, or use a different workflow for long‑read data.",
            title="Long‑Read Not Supported",
            border_style="red"
        )
    )
    sys.exit(1)

# ----------------------------------------------------------------------
# MASTER RULE
# ----------------------------------------------------------------------
rule prok_community_done:
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