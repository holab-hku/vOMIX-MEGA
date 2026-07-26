import os

logdir=relpath("identify/viral/logs")
benchmarks=relpath("identify/viral/benchmarks")
tmpd=relpath("identify/viral/tmp")

email=config["NCBI-email"]
api_key=config["NCBI-API-key"]
nowstr=config["latest-run"]
outdir=config["outdir"]
datadir=config["datadir"]

os.makedirs(logdir, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)

n_cores = config['max-cores']
assembler = config['assembler']
split_part = list(range(1, config["contig-splits"] + 2))
split_part_ids = [f"{i:03d}" for i in range(1, config["contig-splits"] + 2)] # matches seqkit


### Read fasta or fastadir input
if config['fasta'] != "":
  fastap = readfasta(config['fasta'])
  sample_id = config["sample-name"]
  assembly_ids = [sample_id]
elif config['fastadir'] != "":
  fastap = readfastadir(config['fastadir'])
  assembly_ids = config["assembly-ids"]
else:
  samples, assemblies = parse_sample_list(config["samplelist"], datadir, outdir, email, api_key, nowstr)
  fastap = relpath(os.path.join("assembly", assembler, "samples/{sample_id}/output/final.contigs.fa"))
  sample_id = "final.contigs"
  assembly_ids = assemblies.keys()


### MASTER RULE 
rule done_log:
  name: "viral-benchmark.smk Done. removing tmp files"
  localrule: True
  input:
    expand(relpath("identify/viral/samples/{sample_id}/intermediate/genomad/{sample_id}_filtered_summary/{sample_id}_filtered_virus_summary.tsv"), sample_id=assembly_ids), 
    expand(relpath("identify/viral/samples/{sample_id}/tmp/splits/{sample_id}_filtered.part_{part}.fa"), sample_id=assembly_ids, part=split_part_ids),
    expand(relpath("identify/viral/samples/{sample_id}/intermediate/dvf/splits/split-{part}/final_score.txt"), sample_id=assembly_ids, part=split_part_ids),
    expand(relpath("identify/viral/samples/{sample_id}/intermediate/dvf/final_score.txt"), sample_id=assembly_ids),
    expand(relpath("identify/viral/samples/{sample_id}/intermediate/phamer/splits/split-{part}/final_prediction/phamer_prediction.tsv"), sample_id=assembly_ids, part=split_part_ids),
    expand(relpath("identify/viral/samples/{sample_id}/intermediate/phamer/final_prediction/phamer_prediction.tsv"), sample_id=assembly_ids),
    expand(relpath("identify/viral/samples/{sample_id}/intermediate/virsorter2/splits/split-{part}/final-viral-score.tsv"), sample_id=assembly_ids, part=split_part_ids),
    expand(relpath("identify/viral/samples/{sample_id}/intermediate/virsorter2/final-viral-score.tsv"), sample_id=assembly_ids),
    expand(relpath("identify/viral/samples/{sample_id}/intermediate/virfinder/splits/split-{part}/output.tsv"), sample_id=assembly_ids, part=split_part_ids),
    expand(relpath("identify/viral/samples/{sample_id}/intermediate/virfinder/output.tsv"), sample_id=assembly_ids),
    expand(relpath("identify/viral/samples/{sample_id}/intermediate/vibrant/splits/split-{part}/VIBRANT_{sample_id}_filtered.part_{part}/VIBRANT_phages_{sample_id}_filtered.part_{part}/{sample_id}_filtered.part_{part}.phages_combined.txt"), sample_id=assembly_ids, part=split_part_ids),
    expand(relpath("identify/viral/samples/{sample_id}/intermediate/vibrant/VIBRANT_{sample_id}_filtered/VIBRANT_phages_{sample_id}_filtered/{sample_id}_filtered.phages_combined.txt"), sample_id=assembly_ids),
    expand(relpath("identify/viral/samples/{sample_id}/output/viral_benchmark_summary.tsv"), sample_id=assembly_ids), 
    relpath("identify/viral/output/viral_benchmark_merged.csv")
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
  name: "viral-benchmark.smk filter short contigs"
  localrule: True
  input:
    fastap
  output:
    relpath("identify/viral/samples/{sample_id}/tmp/{sample_id}_filtered.fa")
  params:
    minlen=config['contig-min-len'],
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

rule split_contigs:
  name: "viral-benchmark.smk split filtered contigs"
  input:
    relpath("identify/viral/samples/{sample_id}/tmp/{sample_id}_filtered.fa")
  output:
    expand(relpath("identify/viral/samples/{{sample_id}}/tmp/splits/{{sample_id}}_filtered.part_{part}.fa"), part=split_part_ids)
  params:
    parts=config["contig-splits"] + 1,
    tmpdir=os.path.join(tmpd, "contigs/{sample_id}/splits"),
    outdir=relpath("identify/viral/samples/{sample_id}/tmp/splits")
  log: os.path.join(logdir, "splitcontig_{sample_id}.log")
  conda: "../envs/seqkit-biopython.yml"
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input: 4 * 10**3 * attempt
  shell:
    """
    rm -f {output}
    mkdir -p {params.outdir} {params.tmpdir}

    seqkit split2 \
      --by-part {params.parts} \
      --threads {threads} {input} \
      --out-dir {params.tmpdir} 2> {log}

    mv {params.tmpdir}/* {params.outdir}
    """
    

# already has split functionality built in
rule genomad_classify:
  name: "viral-benchmark.smk geNomad classify"
  input:
    fna=relpath("identify/viral/samples/{sample_id}/tmp/{sample_id}_filtered.fa"),
    db=os.path.join(config['genomad-db'], "genomad_db.source")
  output:
    relpath("identify/viral/samples/{sample_id}/intermediate/genomad/{sample_id}_filtered_summary/{sample_id}_filtered_virus_summary.tsv")
  params:
    genomadparams=config['genomad-params'],
    dbdir=config['genomad-db'],
    outdir=relpath("identify/viral/samples/{sample_id}/intermediate/genomad/"),
    splits=config['contig-splits'],
    tmpdir=os.path.join(tmpd, "genomad/{sample_id}")
  log: os.path.join(logdir, "genomad_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "genomad_{sample_id}.log")
  conda: "../envs/genomad.yml"
  threads: 64
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


rule dvf_classify:
  name : "viral-benchmark.smk DeepVirFinder classify"
  input:
    fna=relpath("identify/viral/samples/{sample_id}/tmp/splits/{sample_id}_filtered.part_{part}.fa")
  output:
    relpath("identify/viral/samples/{sample_id}/intermediate/dvf/splits/split-{part}/final_score.txt")
  params:
    script="workflow/software/DeepVirFinder/dvf.py",
    parameters=config['dvf-params'], 
    modeldir="workflow/software/DeepVirFinder/models/",
    outdir=relpath("identify/viral/samples/{sample_id}/intermediate/dvf/splits/split-{part}/"),
    tmpdir=os.path.join(tmpd, "dvf/{sample_id}/split-{part}")
  log: os.path.join(logdir, "dvf_{sample_id}_{part}.log")
  benchmark: os.path.join(benchmarks, "dvf_{sample_id}_{part}.log")
  conda: "../envs/dvf.yml"
  threads: 8
  resources:
    mem_mb=lambda wildcards, attempt, input, threads: max(1 * threads * 10**3 * attempt, 8000)
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

rule dvf_classify_merge:
  name: "viral-benchmark.smk DeepVirFinder merge results"
  localrule: True
  input: 
    expand(relpath("identify/viral/samples/{{sample_id}}/intermediate/dvf/splits/split-{part}/final_score.txt"), part=split_part_ids),
  output: 
    relpath("identify/viral/samples/{sample_id}/intermediate/dvf/final_score.txt")
  params:
    script="workflow/scripts/tables_row_bind.py",
  log: os.path.join(logdir, "dvf_merge_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "dvf_merge_{sample_id}.log")
  conda: "../envs/seqkit-biopython.yml"
  threads: 1
  shell:
    """
    python {params.script} --inputs {input} --output {output}
    """
    

rule phamer_classify:
  name: "viral-benchmark.smk PhaMer classify"
  input:
    fna=relpath("identify/viral/samples/{sample_id}/tmp/splits/{sample_id}_filtered.part_{part}.fa"),
    db=os.path.join(config['PhaBox2-db'], "genus2hostlineage.pkl")
  output:
    phamer=relpath("identify/viral/samples/{sample_id}/intermediate/phamer/splits/split-{part}/final_prediction/phamer_prediction.tsv"), 
    phavip=relpath("identify/viral/samples/{sample_id}/intermediate/phamer/splits/split-{part}/final_prediction/phavip_prediction.tsv")
  params:
    parameters=config['phamer-params'],
    dbdir=config['PhaBox2-db'],
    outdir=relpath("identify/viral/samples/{sample_id}/intermediate/phamer/splits/split-{part}/"),
    tmpdir=os.path.join(tmpd, "phamer/{sample_id}/split-{part}")
  log: os.path.join(logdir, "phamer_{sample_id}_{part}.log")
  benchmark: os.path.join(benchmarks, "phamer_{sample_id}_{part}.log")
  conda: "../envs/phabox2.yml"
  threads: 8
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


rule phamer_classify_merge:
  name: "viral-benchmark.smk PhaMer merge results"
  localrule: True
  input: 
    phamer=expand(relpath("identify/viral/samples/{{sample_id}}/intermediate/phamer/splits/split-{part}/final_prediction/phamer_prediction.tsv"), part=split_part_ids),
    phavip=expand(relpath("identify/viral/samples/{{sample_id}}/intermediate/phamer/splits/split-{part}/final_prediction/phavip_prediction.tsv"), part=split_part_ids),
  output: 
    phamer=relpath("identify/viral/samples/{sample_id}/intermediate/phamer/final_prediction/phamer_prediction.tsv"),
    phavip=relpath("identify/viral/samples/{sample_id}/intermediate/phamer/final_prediction/phavip_prediction.tsv")
  params:
    script="workflow/scripts/tables_row_bind.py",
    outdir=relpath("identify/viral/samples/{sample_id}/intermediate/phamer/final_prediction/"), 
    tmpdir=os.path.join(tmpd, "phamer_merge/{sample_id}")
  log: os.path.join(logdir, "phamer_merge_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "phamer_merge_{sample_id}.log")
  conda: "../envs/seqkit-biopython.yml"
  threads: 1
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} --inputs {input.phamer} --output {params.tmpdir}/phamer.tsv
    python {params.script} --inputs {input.phavip} --output {params.tmpdir}/phavip.tsv

    mv {params.tmpdir}/phamer.tsv {output.phamer}
    mv {params.tmpdir}/phavip.tsv {output.phavip}
    rm -rf {params.tmpdir}
    """


rule virsorter2:
  name: "viral-benchmark.smk VirSorter2 classify"
  input: 
    fna=relpath("identify/viral/samples/{sample_id}/tmp/splits/{sample_id}_filtered.part_{part}.fa"),
    db=os.path.join(config['virsorter2-db'], "Done_all_setup")
  output: 
    relpath("identify/viral/samples/{sample_id}/intermediate/virsorter2/splits/split-{part}/final-viral-score.tsv")
  params: 
    parameters=config['virsorter2-params'],
    dbdir=config['virsorter2-db'],
    outdir=relpath("identify/viral/samples/{sample_id}/intermediate/virsorter2/splits/split-{part}/"),
    tmpdir=os.path.join(tmpd, "virsorter2/{sample_id}/split-{part}")
  log: os.path.join(logdir, "virsorter2_{sample_id}_{part}.log")
  benchmark: os.path.join(benchmarks, "virsorter2_{sample_id}_{part}.log")
  conda: "../envs/virsorter2.yml"
  threads: 8
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

rule virsorter2_merge:
  name: "viral-benchmark.smk VirSorter2 merge results"
  localrule: True
  input: 
    expand(relpath("identify/viral/samples/{{sample_id}}/intermediate/virsorter2/splits/split-{part}/final-viral-score.tsv"), part=split_part_ids),
  output: 
    relpath("identify/viral/samples/{sample_id}/intermediate/virsorter2/final-viral-score.tsv")
  params:
    script="workflow/scripts/tables_row_bind.py",
  log: os.path.join(logdir, "virsorter2_merge_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "virsorter2_merge_{sample_id}.log")
  conda: "../envs/seqkit-biopython.yml"
  threads: 1
  shell:
    """
    python {params.script} --inputs {input} --output {output}
    """

rule virfinder_parallel:
  name: "viral-benchmark.smk VirFinder Parallel run"
  input: relpath("identify/viral/samples/{sample_id}/tmp/splits/{sample_id}_filtered.part_{part}.fa")
  output: relpath("identify/viral/samples/{sample_id}/intermediate/virfinder/splits/split-{part}/output.tsv")
  params: 
    parameters=config['vf-params'],
    outdir=relpath("identify/viral/samples/{sample_id}/intermediate/virfinder/splits/split-{part}/"),
    tmpdir=os.path.join(tmpd, "virfinder/{sample_id}/split-{part}")
  log: os.path.join(logdir, "virfinder_{sample_id}_{part}.log")
  benchmark: os.path.join(benchmarks, "virfinder_{sample_id}_{part}.log")
  conda: "../envs/parallel-virfinder.yml"
  threads: 8
  resources:
    mem_mb=lambda wildcards, attempt, input, threads: max(1 * threads * 10**3 * attempt, 8000)
  shell:
    """
    rm -rf {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    parallel-virfinder.py \
        -i {input} \
        -o {params.tmpdir}/tmp.csv \
        -n {threads} \
        {params.parameters} 2> {log}

    mv {params.tmpdir}/tmp.csv {output}
    rm -rf {params.tmpdir}
    """

rule virfinder_merge:
  name: "viral-benchmark.smk VirFinder merge results"
  localrule: True
  input: 
    expand(relpath("identify/viral/samples/{{sample_id}}/intermediate/virfinder/splits/split-{part}/output.tsv"), part=split_part_ids),
  output: 
    relpath("identify/viral/samples/{sample_id}/intermediate/virfinder/output.tsv")
  params:
    script="workflow/scripts/tables_row_bind.py",
  log: os.path.join(logdir, "virfinder_merge_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "virfinder_merge_{sample_id}.log")
  conda: "../envs/seqkit-biopython.yml"
  threads: 1
  shell:
    """
    python {params.script} --inputs {input} --output {output}
    """

rule VIBRANT:
  name: "viral-benchmark.smk VIBRANT classify"
  input:
    fna=relpath("identify/viral/samples/{sample_id}/tmp/splits/{sample_id}_filtered.part_{part}.fa"),
    db=os.path.join(config['vibrant-db'], "files/VIBRANT_machine_model.sav")
  output: 
    txt=relpath("identify/viral/samples/{sample_id}/intermediate/vibrant/splits/split-{part}/VIBRANT_{sample_id}_filtered.part_{part}/VIBRANT_phages_{sample_id}_filtered.part_{part}/{sample_id}_filtered.part_{part}.phages_combined.txt"),
    tsv=relpath("identify/viral/samples/{sample_id}/intermediate/vibrant/splits/split-{part}/VIBRANT_{sample_id}_filtered.part_{part}/VIBRANT_results_{sample_id}_filtered.part_{part}/VIBRANT_summary_results_{sample_id}_filtered.part_{part}.tsv")
  params:
    parameters=config['vibrant-params'],
    dbdir=config['vibrant-db'],
    outdir=relpath("identify/viral/samples/{sample_id}/intermediate/vibrant/splits/split-{part}/"),
    tmpdir=os.path.join(tmpd, "vibrant/{sample_id}/split-{part}")
  log: os.path.join(logdir, "vibrant_{sample_id}_{part}.log")
  benchmark: os.path.join(benchmarks, "vibrant_{sample_id}_{part}.log")
  conda: "../envs/vibrant.yml"
  threads: 8
  resources:
    mem_mb=lambda wildcards, attempt, input, threads: max(1 * threads * 10**3 * attempt, 8000)
  shell:
    """
    rm -rf {params.outdir} {params.tmpdir}
    mkdir -p {params.outdir} {params.tmpdir}

    VIBRANT_run.py \
        -i {input.fna} \
        -folder {params.tmpdir} \
        -d {params.dbdir}/databases/ \
        -m {params.dbdir}/files/ \
        -f nucl \
        -t {threads} \
        {params.parameters} 2> {log} || true

    # Move results if the directory isn't empty
    if [ -d {params.tmpdir} ] && [ "$(ls -A {params.tmpdir})" ]; then
        mv {params.tmpdir}/* {params.outdir}
    fi

    # Ensure output directories exist
    mkdir -p "$(dirname {output.txt})"
    mkdir -p "$(dirname {output.tsv})"

    # Handle empty/missing TXT file
    if [ ! -s "{output.txt}" ]; then
        touch "{output.txt}"
    fi

    # Handle empty/missing TSV file (write headers)
    if [ ! -s "{output.tsv}" ]; then
        printf "scaffold\ttotal genes\tall KEGG\tKEGG v-score\tall Pfam\tPfam v-score\tall VOG\tVOG v-score\tKEGG int-rep\tKEGG zero\tPfam int-rep\tPfam zero\tVOG redoxin\tVOG rec-tran\tVOG int\tVOG RnR\tVOG DNA\tKEGG restriction check\tKEGG toxin check\tVOG special\tannotation check\tp_v check\tp_k check\tk_v check\tk check\tp check\tv check\th check\n" > "{output.tsv}"
    fi

    rm -rf {params.tmpdir}    
    """

rule merge_VIBRANT:
  name: "viral-benchmark.smk VIBRANT merge results"
  localrule: True
  input: 
    txt=expand(relpath("identify/viral/samples/{{sample_id}}/intermediate/vibrant/splits/split-{part}/VIBRANT_{{sample_id}}_filtered.part_{part}/VIBRANT_phages_{{sample_id}}_filtered.part_{part}/{{sample_id}}_filtered.part_{part}.phages_combined.txt"), part=split_part_ids),
    tsv=expand(relpath("identify/viral/samples/{{sample_id}}/intermediate/vibrant/splits/split-{part}/VIBRANT_{{sample_id}}_filtered.part_{part}/VIBRANT_results_{{sample_id}}_filtered.part_{part}/VIBRANT_summary_results_{{sample_id}}_filtered.part_{part}.tsv"), part=split_part_ids),
  output: 
    txt=relpath("identify/viral/samples/{sample_id}/intermediate/vibrant/VIBRANT_{sample_id}_filtered/VIBRANT_phages_{sample_id}_filtered/{sample_id}_filtered.phages_combined.txt"),
    tsv=relpath("identify/viral/samples/{sample_id}/intermediate/vibrant/VIBRANT_{sample_id}_filtered/VIBRANT_results_{sample_id}_filtered/VIBRANT_summary_results_{sample_id}_filtered.tsv")
  params:
    script="workflow/scripts/tables_row_bind.py",
    outdir_res=relpath("identify/viral/samples/{sample_id}/intermediate/vibrant/VIBRANT_{sample_id}_filtered/VIBRANT_results_{sample_id}_filtered/"),
    outdir_phages=relpath("identify/viral/samples/{sample_id}/intermediate/vibrant/VIBRANT_{sample_id}_filtered/VIBRANT_phages_{sample_id}_filtered/"),
  log: os.path.join(logdir, "vibrant_merge_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "vibrant_merge_{sample_id}.log")
  conda: "../envs/seqkit-biopython.yml"
  threads: 1
  shell:
    """
    rm -rf {params.outdir_res} {params.outdir_phages}
    mkdir -p {params.outdir_res} {params.outdir_phages}
    cat {input.txt} > {output.txt}
    python {params.script} --inputs {input.tsv} --output {output.tsv}
    """

rule merge_viral_summaries:
  name: "viral-benchmark.smk merge all tools"
  input:
    genomad=relpath("identify/viral/samples/{sample_id}/intermediate/genomad/{sample_id}_filtered_summary/{sample_id}_filtered_virus_summary.tsv"),
    dvf=relpath("identify/viral/samples/{sample_id}/intermediate/dvf/final_score.txt"),
    phamer=relpath("identify/viral/samples/{sample_id}/intermediate/phamer/final_prediction/phamer_prediction.tsv"),
    phavip=relpath("identify/viral/samples/{sample_id}/intermediate/phamer/final_prediction/phavip_prediction.tsv"),
    virsorter2=relpath("identify/viral/samples/{sample_id}/intermediate/virsorter2/final-viral-score.tsv"),
    virfinder=relpath("identify/viral/samples/{sample_id}/intermediate/virfinder/output.tsv"),
    vibrant=relpath("identify/viral/samples/{sample_id}/intermediate/vibrant/VIBRANT_{sample_id}_filtered/VIBRANT_results_{sample_id}_filtered/VIBRANT_summary_results_{sample_id}_filtered.tsv"),
  output:
    relpath("identify/viral/samples/{sample_id}/output/viral_benchmark_summary.tsv")
  params:
    script="workflow/scripts/viral_benchmark_merge.py",
    outdir=relpath("identify/viral/samples/{sample_id}/output/"),
    tmpdir=os.path.join(tmpd, "merge_viral_tools/{sample_id}")
  log: os.path.join(logdir, "merge_viral_tools_{sample_id}.log")
  benchmark: os.path.join(benchmarks, "merge_viral_tools_{sample_id}.log")
  conda: "../envs/seqkit-biopython.yml" 
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input: 4 * 10**3 * attempt
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} \
      --genomad {input.genomad} \
      --dvf {input.dvf} \
      --phamer {input.phamer} \
      --phavip {input.phavip} \
      --virsorter2 {input.virsorter2} \
      --virfinder {input.virfinder} \
      --vibrant {input.vibrant} \
      --output {params.tmpdir}/tmp.tsv 2> {log}
    
    mv {params.tmpdir}/tmp.tsv {output}
    rm -rf {params.tmpdir}
    """

rule merge_all_viral_benchmarks:
  name: "viral-benchmark.smk merge all samples to final csv"
  localrule: True
  input:
    expand(relpath("identify/viral/samples/{sample_id}/output/viral_benchmark_summary.tsv"), sample_id=assembly_ids)
  output:
    relpath("identify/viral/output/viral_benchmark_merged.csv")
  params:
    script="workflow/scripts/tables_row_bind.py", 
    sample_ids=list(assembly_ids),
    outdir=relpath("identify/viral/output/")
  log: os.path.join(logdir, "merge_all_viral_benchmarks.log")
  conda: "../envs/seqkit-biopython.yml" 
  threads: 1
  shell:
    """
    mkdir -p {params.outdir}

    python {params.script} \
      --inputs {input} \
      --sample-ids {params.sample_ids} \
      --output {output} 2> {log}
    """