import os 
import json 

from rich.console import Console
from rich.progress import Progress
from rich.layout import Layout
from rich.panel import Panel
console = Console()

configdict = config['viral-identify']
logdir=relpath("identify/viral/logs")
benchmarks=relpath("identify/viral/benchmarks")
tmpd=relpath("identify/viral/tmp")

os.makedirs(logdir, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)

n_cores = config['cores'] 
assembler = config['assembler']

### Read fasta or fastadir input 
if config['fasta'] != "":
  fastap = readfasta(config['fasta'])
  sample_id = config["sample-name"]
  assembly_ids = [sample_id]
elif config['fastadir'] != "":
  fastap = readfastadir(config['fastadir'])
  assembly_ids = config["assembly-ids"]
else:
  samples, assemblies = parse_sample_list(config["samplelist"], datadir, outdir, email, nowstr)
  fastap = relpath(os.path.join("assembly", assembler, "samples/{sample_id}/output/final.contigs.fa"))
  assembly_ids = assemblies.keys()


### MASTER RULE 

rule done_log:
  name: "viral-benchmark.smk Done. removing tmp files"
  localrule: True
  input:
    expand(relpath("identify/viral/samples/{sample_id}/intermediate/genomad/final.contigs.filtered_summary/final.contigs.filtered_virus_summary.tsv"), sample_id=assembly_ids), 
    expand(relpath("identify/viral/samples/{sample_id}/intermediate/dvf/final_score.txt"), sample_id=assembly_ids),
    expand(relpath("identify/viral/samples/{sample_id}/intermediate/phamer/final_prediction/phamer_prediction.tsv"), sample_id=assembly_ids),
    expand(relpath("identify/viral/{sample_id}/intermediate/virsorter2/final-viral-score.tsv"), sample_id=assembly_ids),
    expand(relpath("identify/viral/{sample_id}/intermediate/virfinder/output.tsv"), sample_id=assembly_ids),
  output:
    os.path.join(logdir, "done_benchmarks.log")
  params:
    filteredcontigs=expand(relpath("identify/viral/samples/{sample_id}/tmp"), sample_id=assembly_ids),
    tmpdir=tmpd
  log: os.path.join(logdir, "done_benchmarks.log")
  shell:
    """
    rm -rf {params.tmpdir}/*
    touch {output}
    """


### RULES

rule filter_contigs:
  name: "viral-benchmark.smk filter contigs [length]"
  localrule: True
  input:
    fastap
  output:
    relpath("identify/viral/samples/{sample_id}/tmp/final.contigs.filtered.fa")
  params:
    minlen=configdict['contigminlen'],
    outdir=relpath("identify/viral/samples/{sample_id}/tmp"),
    tmpdir=os.path.join(tmpd, "contigs/{sample_id}")
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
  name: "viral-benchmark.smk geNomad classify" 
  input:
    fna=relpath("identify/viral/samples/{sample_id}/tmp/final.contigs.filtered.fa"),
    db=os.path.join(configdict['genomaddb'], "genomad_db")
  output:
    relpath("identify/viral/samples/{sample_id}/intermediate/genomad/final.contigs.filtered_summary/final.contigs.filtered_virus_summary.tsv")
  params:
    genomadparams=configdict['genomadparams'],
    dbdir=configdict['genomaddb'],
    outdir=relpath("identify/viral/samples/{sample_id}/intermediate/genomad/"),
    tmpdir=os.path.join(tmpd, "genomad/{sample_id}")
  log: os.path.join(logdir, "genomad_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "genomad_{sample_id}.log")
  conda: "../envs/genomad.yml"
  threads: min(32, n_cores//3)
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
        --cleanup \
        {params.genomadparams} &> {log}

    mv {params.tmpdir}/* {params.outdir}
    rm -rf {params.tmpdir}
    """



rule dvf_classify:
  name : "viral-benchmark.smk DeepVirFinder classify"
  input:
    fna=relpath("identify/viral/samples/{sample_id}/tmp/final.contigs.filtered.fa")
  output:
    relpath("identify/viral/samples/{sample_id}/intermediate/dvf/final_score.txt")
  params:
    script="workflow/software/DeepVirFinder/dvf.py",
    parameters=configdict['dvfparams'], 
    modeldir="workflow/software/DeepVirFinder/models/",
    outdir=relpath("identify/viral/samples/{sample_id}/intermediate/dvf/"),
    tmpdir=os.path.join(tmpd, "dvf/{sample_id}")
  log: os.path.join(logdir, "dvf_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "dvf_{sample_id}.log")
  conda: "../envs/dvf.yml"
  threads: min(32, n_cores//3)
  resources:
    mem_mb=lambda wildcards, attempt, input, threads: max(6 * threads * 10**3 * attempt, 8000)
  shell:
    """
    rm -rf {params.tmpdir}/* 
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} \
        -i {input.fna} \
        -l 0 \
        -m {params.modeldir} \
        -c {threads} \
        -o {params.tmpdir} \
        {params.parameters} &> {log}

    mv {params.tmpdir}/* {output}
    rm -rf {params.tmpdir} 
    """


rule phamer_classify:
  name: "viral-benchmark.smk PhaMer classify"
  input:
    fna=relpath("identify/viral/samples/{sample_id}/tmp/final.contigs.filtered.fa")
  output:
    relpath("identify/viral/samples/{sample_id}/intermediate/phamer/final_prediction/phamer_prediction.tsv")
  params:
    parameters=configdict['phamerparams'],
    dbdir=configdict['PhaBox2db'],
    outdir=relpath("identify/viral/samples/{sample_id}/intermediate/phamer/"),
    tmpdir=os.path.join(tmpd, "phamer/{sample_id}")
  log: os.path.join(logdir, "phamer_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "phamer_{sample_id}.log")
  conda: "../envs/phabox2.yml"
  threads: min(32, n_cores//3)
  resources:
    mem_mb=lambda wildcards, attempt, input, threads: max(1 * threads * 10**3 * attempt, 8000)
  shell:
    """
    rm -rf {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    phabox2 --task phamer \
        --contigs {input.fna} \
        --len 0 \
        --threads {threads} \
        --outpth {params.tmpdir} \
        --dbdir {params.dbdir} \
        {params.parameters} &> {log}

    mv -f {params.tmpdir}/* {params.outdir}
    rm -rf {params.tmpdir}
    """



rule virsorter2:
  name: "viral-benchmark.smk VirSorter2 classify"
  input: 
    fna=relpath("identify/viral/samples/{sample_id}/tmp/final.contigs.filtered.fa"), 
    db=os.path.join(configdict['virsorter2db'], "db.tgz")
  output: relpath("identify/viral/{sample_id}/intermediate/virsorter2/final-viral-score.tsv")
  params: 
    parameters=configdict['virsorter2params'],
    dbdir=configdict['virsorter2db'],
    outdir=relpath("identify/viral/{sample_id}/intermediate/virsorter2/"),
    tmpdir=os.path.join(tmpd, "virsorter2/{sample_id}")
  log: os.path.join(logdir, "virsorter2_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "virsorter2_{sample_id}.log")
  conda: "../envs/virsorter2.yml"
  threads: min(32, n_cores//8)
  resources:
    mem_mb=lambda wildcards, attempt, input, threads: max(1 * threads * 10**3 * attempt, 8000)
  shell:
    """
    rm -rf {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    virsorter run \
        -i {input.fna} \
        -w {params.tmpdir} \
        --db-dir {params.dbdir} \
        -j {threads} \
        {params.parameters} \
        all &> {log}

    mv {params.tmpdir}/* {params.outdir}
    """

rule virfinder_parallel:
  name: "viral-benchmark.smk VirFinder Parallel run"
  input: relpath("identify/viral/samples/{sample_id}/tmp/final.contigs.filtered.fa")
  output: relpath("identify/viral/{sample_id}/intermediate/virfinder/output.tsv")
  params: 
    parameters=configdict['vfparams'],
    outdir=relpath("identify/viral/{sample_id}/intermediate/virfinder/"),
    tmpdir=os.path.join(tmpd, "virfinder/{sample_id}")
  log: os.path.join(logdir, "virfinder_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "virfinder_{sample_id}.log")
  conda: "../envs/parallel-virfinder.yml"
  threads: min(32, n_cores//8)
  resources:
    mem_mb=lambda wildcards, attempt, input, threads: max(1 * threads * 10**3 * attempt, 8000)
  shell:
    """
    rm -rf {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    parallel-virfinder.py \
        -i {input} \
        -t {params.tmpdir}/tmp
        -o {params.tmpdir}/tmp.csv \
        -n {threads} \ 
        {params.parameters} 2> {log}

    mv {params.tmpdir}/tmp.csv {output}
    rm -rf {params.tmpdir}
    """

#rule ppr-meta:
#rule seeker:
#rule metaphinder:
#rule VIBRANT: