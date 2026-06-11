logdir = relpath("annotate/viral/logs")
tmpd = relpath("annotate/viral/tmp")
benchmarks = relpath("annotate/viral/benchmarks")

os.makedirs(logdir, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)

n_cores = config['max-cores']

### Read single fasta file if input
if config['fasta'] != "":
  fastap = readfasta(config['fasta'])
  sample_id = config["sample-name"]
  assembly_ids = [sample_id]
else:
  fastap = relpath("identify/viral/output/combined.final.vOTUs.fa")


### MASTER RULE 
rule done_log:
  name: "viral-annotate.smk Done. removing tmp files"
  localrule: True
  input:
    os.path.join(config['eggNOG-db'], "eggnog.db"),
    os.path.join(config['eggNOG-db'], "eggnog_proteins.dmnd"), 
    relpath("annotate/viral/output/proteins.vOTUs.faa"), 
    relpath("annotate/viral/output/PhaVIP/final_prediction/phavip_prediction.tsv"), 
    relpath("annotate/viral/output/eggNOGv2/out.emapper.emapper.annotations"), 
    relpath("annotate/viral/output/MetaCerberus/time.tsv"),
    relpath("annotate/viral/output/Pharokka/pharokka_cds_final_merged_output.tsv"),
    #relpath("annotate/viral/output/VirSorter2/final-viral-score.tsv"), 
    #relpath("annotate/viral/output/VirSorter2/for-dramv/final-viral-combined-for-dramv.fa"),
    #relpath("annotate/viral/output/VirSorter2/for-dramv/viral-affi-contigs-for-dramv.tab"),
    #relpath("annotate/viral/output/DRAMv/annotations.tsv")
  output:
    os.path.join(logdir, "done.log")
  params:
    tmpdir=tmpd
  log: os.path.join(logdir, "done.log")
  shell:
    """
    rm -rf {params.tmpdir}/*
    #touch {output}
    """


### RULES
rule prodigalgv_taxonomy:
  name: "viral-annotate.smk prodigal-gv vTOUs [parallelized]"
  input:
    fastap
  output:
    relpath("annotate/viral/output/proteins.vOTUs.faa")
  params:
    script="workflow/scripts/parallel_prodigal_gv.py",
    outdir=relpath("annotate/viral/output"),
    tmpdir=os.path.join(tmpd, "prodigal-gv")
  conda: "../envs/prodigal-gv.yml"
  log: os.path.join(logdir, "prodigal-gv.log")
  benchmark: os.path.join(benchmarks, "prodigal-gv.log")
  threads: 64
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 72 * 10**3
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} \
        -i {input} \
        -a {params.tmpdir}/tmp.faa \
        -t {threads} &> {log}

    mv {params.tmpdir}/tmp.faa {output}
    rm -rf {params.tmpdir}
    """


rule eggNOGmapper:
  name: "viral-annotate.smk eggNOG-mapper v2 Run"
  input:
    faa=relpath("annotate/viral/output/proteins.vOTUs.faa"),
    db=os.path.join(config['eggNOG-db'], "eggnog.db"),
    diamond=os.path.join(config['eggNOG-db'], "eggnog_proteins.dmnd")
  output:
    relpath("annotate/viral/output/eggNOGv2/out.emapper.emapper.annotations")
  params:
    parameters=config["eggNOG-params"],
    outdir=relpath("annotate/viral/output/eggNOGv2"),
    dbdir=config["eggNOG-db"], 
    tmpdir=os.path.join(tmpd, "eggNOGv2-mapper")
  conda: "../envs/eggnog-mapper.yml"
  log: os.path.join(logdir, "eggNOGv2-mapper.log")
  benchmark: os.path.join(benchmarks, "eggNOGv2-mapper.log")
  threads: 64
  resources: 
    mem_mb=lambda wildcards, attempt: attempt * 64 * 10**3
  shell:
    """
    rm -rf {params.outdir}
    mkdir -p {params.tmpdir}/tmp {params.outdir}

    emapper.py \
        -i {input.faa} \
        --output out.emapper \
        --temp_dir {params.tmpdir}/tmp \
        --output_dir {params.tmpdir} \
        --data_dir {params.dbdir} \
        -m diamond \
        --cpu {threads} \
        {params.parameters} 2> {log}
    
    mv {params.tmpdir}/* {params.outdir}
    """


rule PhaVIP:
  name: "viral-annotate.smk PhaVIP protein annotation"
  input:
    fna=fastap
  output:
    relpath("annotate/viral/output/PhaVIP/final_prediction/phavip_prediction.tsv")
  params:
    parameters=config['PhaVIP-params'],
    dbdir=config['PhaBox2-db'],
    outdir=relpath("annotate/viral/output/PhaVIP"),
    tmpdir=os.path.join(tmpd, "PhaVIP")
  conda: "../envs/phabox2.yml"
  log: os.path.join(logdir, "PhaVIP.log")
  benchmark: os.path.join(benchmarks, "PhaVIP.log")
  threads: 32
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 16 * 10**3
  shell:
    """
    mkdir -p {params.tmpdir} {params.outdir}

    phabox2 --task phavip \
        --contigs {input.fna} \
        --threads {threads} \
        --outpth {params.tmpdir} \
        --dbdir {params.dbdir} \
        {params.parameters} &> {log}

    mv {params.tmpdir}/* {params.outdir}/
    rm -rf {params.tmpdir}
    """


rule virsorter2:
  name: "viral-annotate.smk VirSorter2 annotate for DRAM-v"
  input:
    fna=fastap,
    db=os.path.join(config['virsorter2-db'], "Done_all_setup")
  output: 
    fna=relpath("annotate/viral/output/VirSorter2/for-dramv/final-viral-combined-for-dramv.fa"), 
    tab=relpath("annotate/viral/output/VirSorter2/for-dramv/viral-affi-contigs-for-dramv.tab")
  params:
    parameters=config['virsorter2-annotate-params'],
    dbdir=config['virsorter2-db'],
    outdir=relpath("annotate/viral/output/VirSorter2/"),
    tmpdir=os.path.join(tmpd, "VirSorter2")
  log: os.path.join(logdir, "virsorter2.log")
  benchmark: os.path.join(benchmarks, "virsorter2.log")
  conda: "../envs/virsorter2.yml"
  threads: 64
  resources:
    mem_mb=lambda wildcards, attempt, input, threads: 8 * 10**3 * attempt
  shell:
    """
    rm -rf {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    virsorter run \
        -i {input.fna} \
        -w {params.tmpdir} \
        --prep-for-dramv \
        --db-dir {params.dbdir} \
        -j {threads} \
        {params.parameters} \
        all &> {log}

    mv {params.tmpdir}/* {params.outdir}
    """

rule dramv_annotate:
  name: "viral-annotate.smk DRAM-v annotate"
  input:
    fna=relpath("annotate/viral/output/VirSorter2/for-dramv/final-viral-combined-for-dramv.fa"),
    tab=relpath("annotate/viral/output/VirSorter2/for-dramv/viral-affi-contigs-for-dramv.tab"),
    db=os.path.join(config["dram-db"], "database_processing.log"), 
    json=os.path.join(config["dram-db"], "dram_config.json")
  output: relpath("annotate/viral/output/DRAMv/annotations.tsv")
  params:
    parameters=config["dram-v-annotate-params"],
    dbdir=config["dram-db"],
    outdir=relpath("annotate/viral/output/DRAMv"),
    tmpdir=os.path.join(tmpd, "DRAM-v/annotate")
  log: os.path.join(logdir, "DRAM-v-annotate.log")
  benchmark: os.path.join(benchmarks, "DRAM-v-annotate.log")
  conda: "../envs/dram.yml"
  threads: 64
  resources: 
    mem_mb=lambda wildcards, attempt, input, threads: 64 * 10**3 * attempt
  shell:
    """
    rm -rf {params.outdir} {params.tmpdir}
    mkdir -p {params.outdir}

    DRAM-setup.py import_config --config_loc {input.json}

    DRAM-v.py annotate \
        -i {input.fna} \
        -v {input.tab} \
        -o {params.tmpdir} \
        --threads {threads} \
        {params.parameters} 2> {log}

    mv {params.tmpdir}/* {params.outdir}
    touch {output}
    """

rule dramv_distill:
  name: "viral-annotate.smk DRAM-v distill"
  input: relpath("annotate/viral/output/DRAMv/annotations.tsv")
  output: relpath("annotate/viral/output/DRAMv/distilled/tmpout.log")
  params:
    parameters=config["dram-v-distill-params"],
    outdir=relpath("annotate/viral/output/DRAMv/distilled"),
    tmpdir=os.path.join(tmpd, "DRAM-v/distill")
  log: os.path.join(logdir, "DRAM-v-distill.log")
  benchmark: os.path.join(benchmarks, "DRAM-v-distill.log")
  conda: "../envs/dram.yml"
  threads: 8
  resources:
    mem_mb=lambda wildcards, attempt, input, threads: 8 * 10**3 * attempt
  shell:
    """
    rm -rf {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    DRAM-v.py distill \
        -i {input} \
        -o {params.tmpdir} \
        {params.parameters} 2> {log}

    mv {params.tmpdir}/* {params.outdir}
    touch {output}
    """

rule MetaCerberus:
  name: "viral-annotate.smk MetaCerberus annotate"
  input:
    db=expand(os.path.join(config["metacerberus-db"], "{db}.{extension}"), db = ["PVOG", "VOG"], extension = ["tsv", "hmm.gz"]), 
    fna=fastap
  output:
    relpath("annotate/viral/output/MetaCerberus/time.tsv")
  params:
    parameters=config['metacerberus-params'],
    dbdir=config['metacerberus-db'],
    outdir=relpath("annotate/viral/output/MetaCerberus"),
    tmpdir=os.path.join(tmpd, "MetaCerberus")
  conda: "../envs/metacerberus.yml"
  log: os.path.join(logdir, "metacerberus.log")
  benchmark: os.path.join(benchmarks, "metacerberus.log")
  threads: 64
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 16 * 10**3
  shell:
    """
    rm -rf {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    metacerberus.py \
        --prodigalgv {input.fna} \
        --db-path {params.dbdir} \
        --dir_out {params.tmpdir} \
        --cpus {threads} \
        {params.parameters} 2> {log}

    mv {params.tmpdir}/* {params.outdir}/
    rm -rf {params.tmpdir}
    """

rule Pharokka:
  name: "viral-annotate.smk Pharokka annotate"
  input:
    db=os.path.join(config["pharokka-db"], "all_phrogs.h3m"),
    fna=fastap
  output:
    relpath("annotate/viral/output/Pharokka/pharokka_cds_final_merged_output.tsv")
  params:
    parameters=config['pharokka-params'],
    dbdir=config['pharokka-db'],
    outdir=relpath("annotate/viral/output/Pharokka"),
    tmpdir=os.path.join(tmpd, "Pharokka")
  conda: "../envs/pharokka.yml"
  log: os.path.join(logdir, "pharokka.log")
  benchmark: os.path.join(benchmarks, "pharokka.log")
  threads: 64
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 16 * 10**3
  shell:
    """
    rm -rf {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    pharokka.py \
        -i {input.fna} \
        -d {params.dbdir} \
        -o {params.tmpdir} \
        -f \
        -t {threads} \
        {params.parameters} 2> {log}

    mv {params.tmpdir}/* {params.outdir}/
    rm -rf {params.tmpdir}
    """
