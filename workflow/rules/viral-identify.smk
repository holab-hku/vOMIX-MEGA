import os
import sys
from rich.console import Console
from rich.panel import Panel

console = Console()

logdir = relpath("identify/viral/logs")
benchmarks = relpath("identify/viral/benchmarks")
tmpd = relpath("identify/viral/tmp")

email = config["NCBI-email"]
api_key = config["NCBI-API-key"]
nowstr = config["latest-run"]
outdir = config["outdir"]
datadir = config["datadir"]

os.makedirs(logdir, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)

short_read_assembler = config.get('short-read-assembler', 'megahit')

# ------------------------------------------------------------
# Read input (fasta, fastadir, or sample list)
# ------------------------------------------------------------
if config['fasta'] != "":
    fastap = readfasta(config['fasta'])
    sample_id = config["sample-name"]
    assembly_ids = [sample_id]
    # fastap is a string (single file)
    def get_fasta_input(wildcards):
        return fastap
    fastap_func = get_fasta_input

elif config['fastadir'] != "":
    fasta_files = readfastadir(config['fastadir'])
    assembly_ids = config["assembly-ids"]
    # Build mapping from sample_id (basename) to full path
    # We assume the order matches
    fasta_dict = dict(zip(assembly_ids, fasta_files))
    def get_fasta_input(wildcards):
        return fasta_dict[wildcards.sample_id]
    fastap_func = get_fasta_input

else:
    if config.get("module") == "viral-end-to-end":
        parse_quiet = True
    else:
        parse_quiet = False
    samples, assemblies = parse_sample_list(
        config["samplelist"], datadir, outdir, email, api_key, nowstr, parse_quiet
    )
    assembly_ids = list(assemblies.keys())

    # Function to get assembly path based on read_type
    def get_assembly_path(sample_id):
        rt = assemblies[sample_id]['read_type']
        if rt in ["paired", "single"]:
            assembler = short_read_assembler
        elif rt == "pacbio":
            assembler = "metamdbg"
        elif rt == "nanopore":
            assembler = "nanomdbg"
        else:
            # fallback
            assembler = short_read_assembler
        return relpath(os.path.join("assembly", assembler, "samples", sample_id, "output", "final.contigs.fa"))

    def get_fasta_input(wildcards):
        return get_assembly_path(wildcards.sample_id)

    fastap_func = get_fasta_input

# ------------------------------------------------------------
# MASTER RULE
# ------------------------------------------------------------
rule done_log:
    name: "viral-identify.smk Done. removing tmp files"
    localrule: True
    input:
        relpath("identify/viral/output/combined.final.vOTUs.fa"),
    output:
        os.path.join(logdir, "done.log")
    params:
        tmpdir = tmpd
    log: os.path.join(logdir, "done.log")
    shell:
        """
        rm -rf {params.tmpdir}/*
        touch {output}
        """

# ------------------------------------------------------------
# RULES
# ------------------------------------------------------------
rule filter_contigs:
    name: "viral-identify.smk filter contigs [length]"
    localrule: True
    input:
        fastap_func
    output:
        relpath("identify/viral/samples/{sample_id}/tmp/final.contigs.filtered.fa")
    params:
        minlen = config['contig-min-len'],
        outdir = relpath("identify/viral/samples/{sample_id}/tmp"),
        tmpdir = os.path.join(tmpd, "contigs/{sample_id}")
    log: os.path.join(logdir, "filtercontig_{sample_id}.log")
    conda: "../envs/seqkit-biopython.yml"
    threads: 1
    shell:
        """
        rm -rf {params.tmpdir}/* {params.outdir}
        mkdir -p {params.tmpdir} {params.outdir}

        seqkit seq {input} --min-len {params.minlen} > {params.tmpdir}/tmp.fa

        mv {params.tmpdir}/tmp.fa {output}
        """



rule genomad_classify:
  name: "viral-identify.smk geNomad classify" 
  input: 
    fna=relpath("identify/viral/samples/{sample_id}/tmp/final.contigs.filtered.fa"),
    db=os.path.join(config['genomad-db'], "genomad_db.source")
  output: 
    relpath("identify/viral/samples/{sample_id}/intermediate/genomad/final.contigs.filtered_summary/final.contigs.filtered_virus_summary.tsv")
  params:
    genomadparams=config['genomad-params'],
    dbdir=config['genomad-db'],
    outdir=relpath("identify/viral/samples/{sample_id}/intermediate/genomad/"),
    splits=config['contig-splits'],
    tmpdir=os.path.join(tmpd, "genomad/{sample_id}")
  log: os.path.join(logdir, "genomad_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "genomad_{sample_id}.log")
  conda: "../envs/genomad.yml"
  threads: 32
  resources:
    mem_mb=lambda wildcards, attempt, input: 24 * 10**3 * attempt
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir} 2> {log}
    mkdir -p {params.tmpdir} {params.outdir} 2> {log}

    genomad end-to-end \
        {input.fna} \
        {params.tmpdir} \
        {params.dbdir} \
        --threads {threads} \
        --splits {params.splits} \
        --cleanup \
        {params.genomadparams} &> {log}

    mv {params.tmpdir}/* {params.outdir}
    rm -rf {params.tmpdir}
    """


rule genomad_filter:
  name : "viral-identify.smk filter geNomad output"
  input:
    fna=relpath("identify/viral/samples/{sample_id}/tmp/final.contigs.filtered.fa"),
    tsv=relpath("identify/viral/samples/{sample_id}/intermediate/genomad/final.contigs.filtered_summary/final.contigs.filtered_virus_summary.tsv"),
  output:
    fna=relpath("identify/viral/samples/{sample_id}/output/viral.contigs.fa"),
    scrs=relpath("identify/viral/samples/{sample_id}/output/merged_scores_filtered.csv"),
    hits=relpath("identify/viral/samples/{sample_id}/output/viralhits_list")
  params:
    script="workflow/scripts/genomad_filter.py", 
    minlen=config['genomad-min-len'],
    cutoff=config['genomad-cutoff'],
    outdir=relpath("identify/viral/samples/{sample_id}/output/"),
    tmpdir=os.path.join(tmpd, "filter/{sample_id}")
  log: os.path.join(logdir, "genomad_filter_{sample_id}.log")
  conda: "../envs/seqkit-biopython.yml"
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
  shell:
    """
    rm -rf {params.tmpdir}/*
    mkdir -p {params.tmpdir} {params.outdir}
    
    python {params.script} \
        --genomad_out {input.tsv} \
        --genomad_min_score {params.cutoff} \
        --genomad_min_len {params.minlen} \
        --output_path {params.tmpdir}/tmp.csv \
        --hitlist_path {params.tmpdir}/tmplist &> {log}

    mv {params.tmpdir}/tmp.csv {output.scrs}
    cat {params.tmpdir}/tmplist | uniq > {output.hits}
    
    seqkit grep {input.fna} -f {output.hits} | seqkit replace -p  "\\s.*" -r "" | seqkit replace -p $ -r _{wildcards.sample_id}  > {params.tmpdir}/tmp.fa 2> {log}
    mv {params.tmpdir}/tmp.fa {output.fna}

    rm -rf {params.tmpdir}/*
    """
    
 

rule cat_contigs:
  name : "viral-identify.smk combine viral contigs"
  input:
    fna=expand(relpath("identify/viral/samples/{sample_id}/output/viral.contigs.fa"), sample_id=assembly_ids),
    scrs=expand(relpath("identify/viral/samples/{sample_id}/output/merged_scores_filtered.csv"), sample_id=assembly_ids)
  output: 
    fna=relpath("identify/viral/intermediate/scores/combined.viralcontigs.fa"),
    scrs=relpath("identify/viral/intermediate/scores/combined_viral_scores.csv")
  params:
    script="workflow/scripts/viral_merge_scores.py", 
    names=list(assembly_ids),
    tmpdir=tmpd
  log: os.path.join(logdir, "catcontigs.log")
  conda: "../envs/seqkit-biopython.yml"
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir}
    
    echo "{params.names}" > {params.tmpdir}/tmp.names
    echo "{input.scrs}" > {params.tmpdir}/tmp.csv.paths
    
    python {params.script} \
        --names {params.tmpdir}/tmp.names \
        --csvs {params.tmpdir}/tmp.csv.paths > {params.tmpdir}/tmp.csv 2> {log}
    cat {input.fna} > {params.tmpdir}/tmp.fa 2> {log}

    mv {params.tmpdir}/tmp.csv {output.scrs}
    mv {params.tmpdir}/tmp.fa {output.fna}

    rm -rf {params.tmpdir}
    """


# THEN GOES INTO 
# 1) VIRAL BINNING [optional]
# 2) CLUSTERING [sensitive or fast]
# 3) CHECKV-PYHMMER 
# THEN COMES BACK HERE TO GET vCONTIGS


rule combine_classifications:
  name: "viral-identify.smk combine derepped classification results"
  input:
    checkv_out=relpath("identify/viral/output/checkv/quality_summary.tsv"),
    classify_out=relpath("identify/viral/intermediate/scores/combined_viral_scores.csv")
  output:
    relpath("identify/viral/output/checkv/combined_classification_results.csv")
  params:
    script="workflow/scripts/viral_combine_class.py",
    tmpdir=tmpd
  log: os.path.join(logdir, "combine_classification.log")
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
  conda : "../envs/seqkit-biopython.yml"
  shell:
    """
    rm -rf {params.tmpdir}/*
    mkdir -p {params.tmpdir}

    python {params.script} \
        --mergedclassify {input.classify_out} \
        --checkvsummary {input.checkv_out} \
        --output {params.tmpdir}/tmp.csv 2> {log}

    mv {params.tmpdir}/tmp.csv {output}
    """


rule consensus_filtering:
  name: "viral-identify.smk consensus vOTU filtering"
  input:
    relpath("identify/viral/output/checkv/combined_classification_results.csv")
  output:
    summary=relpath("identify/viral/output/classification_summary_vOTUs.csv"),
    proviruslist=relpath("identify/viral/output/provirus.list.txt"),
    viruslist=relpath("identify/viral/output/virus.list.txt")
  params:
    script="workflow/scripts/genomad_consensus_filtering.py",
    genomad=config['genomad-cutoff-s'],
    tmpdir=tmpd
  log: os.path.join(logdir, "consensus_filtering.log")
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
  conda: "../envs/seqkit-biopython.yml"
  shell:
    """
    rm -rf {params.tmpdir}/*
    mkdir -p {params.tmpdir}

    python {params.script} \
        --classification_results {input} \
        --genomad_min_score {params.genomad} \
        --summary_out {params.tmpdir}/tmp.csv \
        --provirus_list {params.tmpdir}/tmp.list.1 \
        --virus_list {params.tmpdir}/tmp.list.2  2> {log}

    mv {params.tmpdir}/tmp.csv {output.summary}
    mv {params.tmpdir}/tmp.list.1 {output.proviruslist}
    mv {params.tmpdir}/tmp.list.2 {output.viruslist}
    """


rule votu:
  name: "viral-identify.smk generate final vContigs"
  input:
    provirusfasta=relpath("identify/viral/output/checkv/proviruses.fna"),
    virusfasta=relpath("identify/viral/output/checkv/viruses.fna"), 
    provirushits=relpath("identify/viral/output/provirus.list.txt"),
    virushits=relpath("identify/viral/output/virus.list.txt")
  output:
    combined=relpath("identify/viral/output/combined.final.vOTUs.fa"),
    provirus=relpath("identify/viral/output/provirus.final.vOTUs.fa"),
    virus=relpath("identify/viral/output/virus.final.vOTUs.fa"),
    tsv=relpath("identify/viral/output/GC_content_vOTUs.tsv")
  params:
    outdir=relpath("identify/viral/output/"),
    tmpdir=tmpd
  log: os.path.join(logdir, "vOTUs.log")
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
  conda: "../envs/seqkit-biopython.yml"
  shell:
    """
    rm -rf {params.tmpdir}/*
    mkdir -p {params.tmpdir} {params.outdir}

    seqkit replace {input.provirusfasta} --f-use-regexp -p "(.+)_\\d\\s.+$" -r '$1' | \
        seqkit grep -f {input.provirushits} > {params.tmpdir}/tmp1.fa 2> {log}
    seqkit replace {input.virusfasta} --f-use-regexp -p "(.+)_\\d\\s.+$" -r '$1' | \
        seqkit grep -f {input.provirushits} >> {params.tmpdir}/tmp1.fa 2> {log}

    seqkit grep {input.virusfasta} -f {input.virushits} > {params.tmpdir}/tmp2.fa 2> {log}
    seqkit grep {input.provirusfasta} -f {input.virushits} >> {params.tmpdir}/tmp2.fa 2> {log}

    cat {params.tmpdir}/tmp1.fa {params.tmpdir}/tmp2.fa  > {params.tmpdir}/tmp3.fa 2> {log}

    seqkit rmdup {params.tmpdir}/tmp1.fa > {output.provirus} 2> {log}
    seqkit rmdup {params.tmpdir}/tmp2.fa > {output.virus} 2> {log}
    seqkit rmdup {params.tmpdir}/tmp3.fa > {output.combined} 2> {log}
    
    echo -e "seq_id\tlength\tGC_percent" > {params.tmpdir}/tmp.tsv
    seqkit fx2tab -g -l -n {output.combined} >> {params.tmpdir}/tmp.tsv
    mv {params.tmpdir}/tmp.tsv {output.tsv}

    rm -rf {params.tmpdir}/*
    """


rule benchmark_summary:
  name: "viral-identify.smk summaries performance benchmarks"
  localrule: True
  input:
    os.path.join(benchmarks, "checkv.log")
  output:
    os.path.join(benchmarks, "summary.tsv")
  params:
    indir=benchmarks
  threads: 1
  shell:
    """
    cat {params.indir}/* | head -n 1 > {output}
    for file in $(find {params.indir}/*.log -type f); do echo -e "$(tail -n +2 $file)\t  $(basename $file)" >> {output}; done
    """
