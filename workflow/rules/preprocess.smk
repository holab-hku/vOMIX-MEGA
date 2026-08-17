# ----------------------------------------------------------------------
# Configuration & setup
# ----------------------------------------------------------------------
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
# Validation: Bowtie2 cannot handle long reads
# ----------------------------------------------------------------------
if config.get("decontam-host", False):
    aligner = config.get("hostile-aligner", "minimap2")
    if aligner == "bowtie2":
        long_read_samples = [
            s for s, info in samples.items()
            if info.get("read_type") in ["pacbio", "nanopore"]
        ]
        if long_read_samples:
            console.print(
                Panel.fit(
                    f"[bold red]ERROR:[/] Bowtie2 is not suitable for long reads (PacBio/Nanopore).\n"
                    f"Found long-read samples: {long_read_samples}\n"
                    f"Please set [cyan]hostile-aligner: minimap2[/] in your config for long-read data.\n"
                    f"Alternatively, set [cyan]decontam-host: false[/] to skip host removal.",
                    title="Hostile Aligner Incompatibility",
                    border_style="red",
                )
            )
            sys.exit(1)

os.makedirs(logdir, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)


def retrieve_accessions(wildcards):
    try:
        acc = samples[wildcards.sample_id]["accession"]
    except KeyError:
        acc = wildcards.sample_id
    return acc


# ----------------------------------------------------------------------
# MASTER RULE
# ----------------------------------------------------------------------
if config["dwnld-only"]:
    rule preprocess_done:
        name: "preprocessing.py download SRA only Done."
        localrule: True
        input:
            expand(os.path.join(datadir, "{sample_id}_{i}.fastq.gz"), sample_id=samples.keys(),i=[1, 2])
        output:
            os.path.join(logdir, "done.log")
        shell: "touch {output}"

elif config["keep-intermediates"]:
    rule preprocess_done:
        name: "preprocessing.py Done. deleting all tmp files"
        localrule: True
        input:
            expand(relpath("preprocess/samples/{sample_id}/{sample_id}_R{i}.fastq.gz"), sample_id=samples.keys(), i=[1, 2]),
            expand(os.path.join(datadir, "{sample_id}_{i}.fastq.gz"), sample_id=samples.keys(), i=[1, 2]),
            expand(relpath("preprocess/samples/{sample_id}/output/{sample_id}_R{i}_cut.trim.filt.fastq.gz"), sample_id=samples.keys(), i=[1, 2]),
            relpath("preprocess/reports/preprocess_report.html"),
            relpath("preprocess/reports/library_size_stats.csv"),
        output:
            os.path.join(logdir, "done.log")
        shell: "touch {output}"

else:
    rule preprocess_done:
        name: "preprocessing.py Done. deleting all tmp and intermediate files."
        localrule: True
        input:
            expand(relpath("preprocess/samples/{sample_id}/{sample_id}_R{i}.fastq.gz"), sample_id=samples.keys(), i=[1, 2],),
            expand(os.path.join(datadir, "{sample_id}_{i}.fastq.gz"),sample_id=samples.keys(),i=[1, 2],),
            expand(relpath("preprocess/samples/{sample_id}/output/{sample_id}_R{i}_cut.trim.filt.fastq.gz"),sample_id=samples.keys(),i=[1, 2],),
            relpath("preprocess/reports/preprocess_report.html"),
            relpath("preprocess/reports/library_size_stats.csv"),
        output:
            os.path.join(logdir, "done.log")
        params:
            intermediate=expand(
                relpath(
                    "preprocess/samples/{sample_id}/output/{sample_id}_R{i}_cut.trim.filt.nodecontam.fastq.gz"
                ),
                sample_id=samples.keys(),
                i=[1, 2],
            )
        shell:
            """
            rm -f {params.intermediate}
            touch {output}
            """


# ----------------------------------------------------------------------
# DOWNLOAD
# ----------------------------------------------------------------------
rule download_fastq:
    name: "preprocessing.py download fastq from SRA"
    output:
        R1=os.path.join(datadir, "{sample_id}_1.fastq.gz"),
        R2=os.path.join(datadir, "{sample_id}_2.fastq.gz")
    params:
        download=config["dwnld-params"],
        pigz=config["pigz-params"],
        logdir=os.path.join(datadir, ".log"),
        accessions=lambda wildcards: retrieve_accessions(wildcards),
        tmpdir=os.path.join(datadir, ".tmp/{sample_id}"),
        read_type=lambda wildcards: samples[wildcards.sample_id]["read_type"],
    log: os.path.join(datadir, ".log/{sample_id}.log")
    benchmark: os.path.join(benchmarks, "download_{sample_id}.log")
    conda: "../envs/sratools-pigz.yml"
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4 * 10**3
    shell:
        """
        mkdir -p {params.tmpdir} {params.logdir}

        fasterq-dump {params.accessions} \
            {params.download} \
            --split-3 \
            --skip-technical \
            --outdir {params.tmpdir} \
            --temp {params.tmpdir} \
            --threads {threads} &> {log}

        pigz -p {threads} -c {params.pigz} {params.tmpdir}/*_1.fastq > {output.R1} 2>> {log}

        if [ "{params.read_type}" == "paired" ]; then
            pigz -p {threads} -c {params.pigz} {params.tmpdir}/*_2.fastq > {output.R2} 2>> {log}
        else
            touch {output.R2}
        fi

        rm -rf {params.tmpdir}
        """


# ----------------------------------------------------------------------
# QC & TRIMMING – short reads use fastp, long reads use fastplong
# ----------------------------------------------------------------------
rule fastp:
    name: "preprocessing.py preprocess reads"
    input:
        R1=os.path.join(datadir, "{sample_id}_1.fastq.gz"),
        R2=os.path.join(datadir, "{sample_id}_2.fastq.gz")
    output:
        R1=relpath("preprocess/samples/{sample_id}/output/{sample_id}_R1_cut.trim.filt.fastq.gz"),
        R2=relpath("preprocess/samples/{sample_id}/output/{sample_id}_R2_cut.trim.filt.fastq.gz"),
        html=relpath("preprocess/samples/{sample_id}/report.fastp.html"),
        json=relpath("preprocess/samples/{sample_id}/report.fastp.json")
    params:
        fastp_paired=config["fastp-params"],
        fastp_single=config.get("fastp-single-params", ""),
        fastplong_params=config.get("fastplong-params", "--min_length 1000"),
        outdir=relpath("preprocess/samples/{sample_id}/output"),
        tmpdir=os.path.join(tmpd, "fastp/{sample_id}"),
        read_type=lambda wildcards: samples[wildcards.sample_id]["read_type"],
    log: os.path.join(logdir, "fastp_{sample_id}.log")
    benchmark: os.path.join(benchmarks, "fastp_{sample_id}.log")
    threads: 12
    resources:
        mem_mb=lambda wildcards, input, attempt: attempt * max(5 * input.size_mb, 4000)
    conda: "../envs/fastp.yml"   # assumes fastplong is in the same env
    shell:
        """
        rm -rf {params.tmpdir}/*
        mkdir -p {params.tmpdir} {params.outdir}

        if [ "{params.read_type}" == "paired" ]; then
            fastp -i {input.R1} -I {input.R2} \
                -o {params.tmpdir}/R1.fastq.gz \
                -O {params.tmpdir}/R2.fastq.gz \
                --thread {threads} \
                --html {params.tmpdir}/tmp.html \
                --json {params.tmpdir}/tmp.json \
                {params.fastp_paired} &> {log}
        elif [ "{params.read_type}" == "single" ]; then
            fastp -i {input.R1} \
                -o {params.tmpdir}/R1.fastq.gz \
                --thread {threads} \
                --html {params.tmpdir}/tmp.html \
                --json {params.tmpdir}/tmp.json \
                {params.fastp_single} &> {log}
            touch {params.tmpdir}/R2.fastq.gz
        elif [[ "{params.read_type}" == "pacbio" ]] || [[ "{params.read_type}" == "nanopore" ]]; then
            # Use fastplong for long-read QC
            fastplong -i {input.R1} \
                -o {params.tmpdir}/R1.fastq.gz \
                --thread {threads} \
                --html {params.tmpdir}/tmp.html \
                --json {params.tmpdir}/tmp.json \
                {params.fastplong_params} &> {log}
            touch {params.tmpdir}/R2.fastq.gz
        else
            echo "ERROR: Unknown read_type '{params.read_type}'" &> {log}
            exit 1
        fi

        mv {params.tmpdir}/R1.fastq.gz {output.R1}
        mv {params.tmpdir}/R2.fastq.gz {output.R2}
        mv {params.tmpdir}/tmp.html {output.html}
        mv {params.tmpdir}/tmp.json {output.json}

        rm -rf {params.tmpdir}
        """


# ----------------------------------------------------------------------
# HOST DECONTAMINATION – runs for all read types (hostile supports long reads)
# ----------------------------------------------------------------------
if config["decontam-host"]:
    # Index paths
    if config["hostile-aligner"] == "minimap2":
        index_path = os.path.join(config["hostile-index-db"], config["hostile-index-name"] + ".fa.gz")
        index_check = os.path.join(config["hostile-index-db"], config["hostile-index-name"] + ".fa.gz")
    elif config["hostile-aligner"] == "bowtie2":
        index_path = os.path.join(config["hostile-index-db"], config["hostile-index-name"])
        index_check = expand(
            os.path.join(config["hostile-index-db"], config["hostile-index-name"]) + ".{i}.bt2",
            i=[1, 2, 3, 4],
        )

    rule decontam:
        name: "preprocess.py Hostile host decontamination"
        input:
            R1=relpath("preprocess/samples/{sample_id}/output/{sample_id}_R1_cut.trim.filt.fastq.gz"),
            R2=relpath("preprocess/samples/{sample_id}/output/{sample_id}_R2_cut.trim.filt.fastq.gz"),
            indexdb=index_check,
        output:
            R1=relpath("preprocess/samples/{sample_id}/output/{sample_id}_R1_cut.trim.filt.fastq.gz"),
            R2=relpath("preprocess/samples/{sample_id}/output/{sample_id}_R2_cut.trim.filt.fastq.gz"),
        params:
            parameters=config["hostile-params"],
            aligner=config["hostile-aligner"],
            alignerp=config["hostile-aligner-params"],
            indexpath=index_path,
            outdir=relpath("preprocess/samples/{sample_id}/output"),
            tmpdir=os.path.join(tmpd, "hostile/{sample_id}"),
            read_type=lambda wildcards: samples[wildcards.sample_id]["read_type"],
        log: os.path.join(logdir, "hostile_{sample_id}.log")
        benchmark: os.path.join(benchmarks, "hostile_{sample_id}.log")
        threads: 8
        resources:
            mem_mb=lambda wildcards, input, attempt: attempt * 16 * 10**3
        conda: "../envs/hostile.yml"
        shell:
            """
            rm -rf {params.tmpdir}
            mkdir -p {params.tmpdir} {params.outdir}

            # hostile accepts both paired and single; long-reads use minimap2 (single-end)
            if [ "{params.read_type}" == "paired" ]; then
                hostile clean \
                    --fastq1 {input.R1} \
                    --fastq2 {input.R2} \
                    --aligner {params.aligner} \
                    --aligner-args "{params.alignerp}" \
                    --index {params.indexpath} \
                    --threads {threads} \
                    --output {params.tmpdir} &> {log}
                mv {params.tmpdir}/{wildcards.sample_id}_R1_cut.trim.filt.nodecontam.clean_1.fastq.gz {output.R1}
                mv {params.tmpdir}/{wildcards.sample_id}_R2_cut.trim.filt.nodecontam.clean_2.fastq.gz {output.R2}
            else
                # single / pacbio / nanopore: only one input
                hostile clean \
                    --fastq1 {input.R1} \
                    --aligner {params.aligner} \
                    --aligner-args "{params.alignerp}" \
                    --index {params.indexpath} \
                    --threads {threads} \
                    --output {params.tmpdir} &> {log}
                mv {params.tmpdir}/{wildcards.sample_id}_R1_cut.trim.filt.nodecontam.clean_1.fastq.gz {output.R1}
                touch {output.R2}
            fi

            rm -rf {params.tmpdir}
            """

else:
    # No decontam: we just pass the fastp output through; the done rule will see the final files.
    pass


# ----------------------------------------------------------------------
# SYMLINK
# ----------------------------------------------------------------------
rule symlink:
    name: "preprocessing.py creating symbolic links"
    localrule: True
    input:
        R1=relpath("preprocess/samples/{sample_id}/output/{sample_id}_R1_cut.trim.filt.fastq.gz"),
        R2=relpath("preprocess/samples/{sample_id}/output/{sample_id}_R2_cut.trim.filt.fastq.gz"),
    output:
        R1=relpath("preprocess/samples/{sample_id}/{sample_id}_R1.fastq.gz"),
        R2=relpath("preprocess/samples/{sample_id}/{sample_id}_R2.fastq.gz"),
    shell:
        """
        ln -s {input.R1} {output.R1}
        ln -s {input.R2} {output.R2}
        """


# ----------------------------------------------------------------------
# AGGREGATE FASTP STATS (includes all samples, updated parser handles long-read JSON)
# ----------------------------------------------------------------------
rule aggregate_fastp:
    name: "preprocess.py summarise fastp stats"
    localrule: True
    input:
        jsons=expand(
            relpath("preprocess/samples/{sample_id}/report.fastp.json"),
            sample_id=samples.keys(),
        )
    output:
        relpath("preprocess/reports/library_size_stats.csv")
    params:
        script="workflow/scripts/fastp_parse.py",
        names=list(samples.keys()),
        outdir=relpath("preprocess/reports"),
        tmpdir=os.path.join(tmpd, "fastp/summary"),
    log: os.path.join(logdir, "fastp_summary_stats.log")
    benchmark: os.path.join(benchmarks, "fastp_summary_stats.log")
    conda: "../envs/seqkit-biopython.yml"
    threads: 1
    shell:
        """
        rm -rf {params.tmpdir} {params.outdir}
        mkdir -p {params.tmpdir} {params.outdir}

        echo "{params.names}" > {params.tmpdir}/tmp.names
        echo "{input.jsons}" > {params.tmpdir}/tmp.jsons

        python {params.script} \
            --names {params.tmpdir}/tmp.names \
            --jsons {params.tmpdir}/tmp.jsons > {params.tmpdir}/tmp.csv 2> {log}

        mv {params.tmpdir}/tmp.csv {output}
        rm -rf {params.tmpdir}/*
        """


# ----------------------------------------------------------------------
# MULTIQC REPORT – scans everything; long-read JSONs are harmless
# ----------------------------------------------------------------------
rule multiqc:
    name: "preprocessing.py MultiQC preprocess report"
    input:
        R1s=expand(
            relpath("preprocess/samples/{sample_id}/output/{sample_id}_R1_cut.trim.filt.fastq.gz"),
            sample_id=samples.keys(),
        ),
        R2s=expand(
            relpath("preprocess/samples/{sample_id}/output/{sample_id}_R2_cut.trim.filt.fastq.gz"),
            sample_id=samples.keys(),
        ),
        logs=expand(
            relpath("preprocess/samples/{sample_id}/report.fastp.json"),
            sample_id=samples.keys(),
        ),
    output:
        relpath("preprocess/reports/preprocess_report.html"),
        relpath("preprocess/reports/preprocess_report_data/multiqc.log"),
    params:
        searchdir=relpath("preprocess/"),
        outdir=relpath("preprocess/reports"),
        tmpdir=os.path.join(tmpd, "multiqc"),
    log: os.path.join(logdir, "multiqc.log")
    benchmark: os.path.join(benchmarks, "multiqc.log")
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: attempt * 8 * 10**3
    conda: "../envs/multiqc.yml"
    shell:
        """
        rm -rf {params.tmpdir}
        mkdir -p {params.tmpdir} {params.outdir}

        multiqc {params.searchdir} -f -o {params.tmpdir} -n preprocess_report.html 2> {log}
        mv {params.tmpdir}/*.html {params.outdir}
        mv {params.tmpdir}/preprocess* {params.outdir}

        rm -rf {params.tmpdir}/*
        """