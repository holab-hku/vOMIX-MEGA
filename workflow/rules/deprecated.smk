##############################
# VIRSORTER 2 CLASSIFICATION #
##############################


#rule vs2_classify:
#  input:
#    os.path.join(config['contigdir'], "{sample_id}/final.contigs.fa")
#  output:
#    "output/viralcontigident/{sample_id}/intermediate/vs2/final-viral-score.tsv"
#  params:
#    vs2params = config['vs2params'],
#    db_dir = config['vs2db'],
#    output_dir = "output/viralcontigident/{sample_id}/intermediate/vs2/",
#    tmp_dir = "$TMPDIR/{sample_id}"
#  log: "logs/viralcontigident_{sample_id}_vs2.log"
#  benchmark: "benchmarks/viralcontigident_{sample_id}_vs2.log"
#  conda: "pipeline/envs/vs2.yml"
#  threads: 8
#  resources:
#    runtime = lambda wildcards, attempt: attempt*attempt*60
#  shell:
#    """
#    rm -rf {params.output_dir}
#    mkdir -p {params.tmp_dir}
#
#    virsorter run \
#        -i {input} \
#        -w {params.tmp_dir} \
#        --db-dir {params.db_dir} \
#        -j {threads} \
#        {params.vs2params} \
#        all &> {log}
#
#    mkdir -p {params.output_dir}
#    mv {params.tmp_dir}/* {params.output_dir}
#    """

#########################
# VIRBOT CLASSIFICATION #
#########################


#rule virbot_classify:
#  input:
#    os.path.join(config['contigdir'], "{sample_id}/final.contigs.fa")
#  output:
#    "output/viralcontigident/{sample_id}/intermediate/virbot/pos_contig_score.csv"
#  params:
#    virbotparams = config['virbotparams'],
#    output_dir = "output/viralcontigident/{sample_id}/intermediate/virbot/",
#    tmp_dir = "$TMPDIR/_virbot/{sample_id}",
#    tmp_dir_parent = "$TMPDIR/_virbot"
#  log: "logs/viralcontigident_{sample_id}_virbot.log"
#  benchmark: "benchmarks/viralcontigident_{sample_id}_virbot.log"
#  conda: "pipeline/envs/virbot.yml"
#  threads: 8
#  shell:
#    """
#    rm -rf {params.tmp_dir} # remove any old directories
#    rm -rf {params.output_dir}
#    mkdir -p {params.tmp_dir_parent}l
#    mkdir -p {params.output_dir}
#
#    python pipeline/bin/VirBot/VirBot.py \
#        --input {input} \
#        --output {params.tmp_dir} \
#        --threads {threads} \
#        {params.virbotparams} &> {log}
#
#    #mv {params.tmp_dir}/tmp {params.tmp_dir}/intermediate
#    mv -f {params.tmp_dir}/* {params.output_dir}
#    rm -r {params.tmp_dir}
#    """

#rule phabox_classify:
#  input:
#    os.path.join(config['contigdir'], "{sample_id}/final.contigs.fa")
#  output:
#    "output/phaboxout/{sample_id}/phamer/out/phamer_prediction.csv"
#  params:
#    phamerparams = config['phamerparams'],
#    db_dir = config['phamerdb'],
#    params_dir = "pipeline/params/phabox/",
#    output_dir = "output/phaboxout/{sample_id}/",
#    script_dir = "pipeline/bin/PhaBOX/scripts/",
#    tmp_dir = "$TMPDIR/{sample_id}"
#  log: "logs/viralcontigident_{sample_id}_phabox.log"
#  conda: "pipeline/envs/phabox.yml"
#  threads: 32
#  shell:
#    """
#    rm -rf {params.output_dir}
#    mkdir -p {params.output_dir}
#    mkdir -p {params.tmp_dir}
#
#    python pipeline/bin/PhaBOX/main.py \
#        --contigs {input} \
#        --threads {threads} \
#        --rootpth {params.tmp_dir} \
#        --dbdir {params.db_dir} \
#        --parampth {params.params_dir} \
#        --scriptpth {params.script_dir} \
#        {params.phamerparams} &> {log}
#
#    mv -f {params.tmp_dir}/* {params.output_dir}
#    """


#rule pseudophabox:
#  input:
#    expand("output/phaboxout/{sample_id}/phamer/out/phamer_prediction.csv", sample_id = samples.keys())
#  output:
#    "test.txt"
#  shell: "touch {output}"


rule bowtie2_build:
  name: "viral-binning.py bowtie2 build viral contig index"
  input:
    relpath("viralcontigident/intermediate/scores/combined.viralcontigs.fa")
  output:
    expand(relpath("binning/viral/intermediate/combined.viralcontigs.fa.{i}.bt2l"), i=range(1,5))
  params:
    parameters=config["bowtiebuildparams"],
    outdir=relpath("binning/viral/intermediate"), 
    tmpdir=os.path.join(tmpd, "build")
  log: os.path.join(logdir, "bowtie2build.log")
  conda: "../envs/bowtie2.yml"
  threads: 32 
  shell:
    """
    rm -rf {params.tmpdir}/* {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    bowtie2-build \
        --large-index \
        --threads {threads} \
        {params.parameters} \
        {input} \
        {params.tmpdir}/combined.viralcontigs.fa &> {log}

    mv {params.tmpdir}/combined.viralcontigs.fa.*.bt2l {params.outdir}
    
    rm -rf {params.tmpdir}/*
    """

rule bowtie2:
  name: "viral-binning.py bowtie2 map viral contigs"
  input:
    R1=relpath("preprocess/samples/{sample_id}/output/{sample_id}_R1_cut.trim.filt.fastq.gz"),
    R2=relpath("preprocess/samples/{sample_id}/output/{sample_id}_R2_cut.trim.filt.fastq.gz"),
    bowtie=expand(relpath("binning/viral/intermediate/combined.viralcontigs.fa.{i}.bt2l"), i=range(1,5))
  output:
    relpath("binning/viral/samples/{sample_id}/intermediate/{sample_id}.bam")
  params:
    bowtie2params=config["bowtie2params"],
    prefix=relpath("binning/viral/intermediate/combined.viralcontigs.fa"),
    outdir=relpath("binning/viral/samples/{sample_id}/intermediate"),
    tmpdir=os.path.join(tmpd, "bowtie2/{sample_id}")
  log: os.path.join(logdir, "bowtie2_{sample_id}.log")
  conda: "../envs/bowtie2.yml"
  threads: 8
  resources:
  shell:
    """
    rm -rf {params.tmpdir} {output}
    mkdir -p {params.tmpdir} {params.outdir}

    bowtie2 \
        -1 {input.R1} -2 {input.R2} \
        -p {threads} \
        -x {params.prefix} \
        {params.bowtie2params}  2> {log} | samtools view -hb - 2> {log} | samtools sort - -o {params.tmpdir}/tmp.bam 2> {log}
    
    mv {params.tmpdir}/tmp.bam {output}
    """


rule bowtie2_build:
  name: "prokaryote.smk Bowtie2 build insex"
  input:
    relpath("output/assembly/samples/{assembly_id}/output/final.contigs.fa")
  output:
    expand(relpath("binning/prokaryotic/samples/{assembly_id}/{assembly_id}.fa.{i}.bt2l"), i=range(1,5))
  params:
    parameters=configdict["bowtiebuildparams"],
    outdir=relpath("binning/prokaryotic/samples/{assembly_id}"),
    tmpdir=os.path.join(tmpd, "bowtie2_build")
  log: os.path.join(logdir, "{assembly_id}_bowtie2build.log")
  conda: "../envs/bowtie2.yml"
  threads: 8
  shell:
    """
    rm -rf {params.tmpdir}/* {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    bowtie2-build \
        --large-index \
        --threads {threads} \
        {params.parameters} \
        {input} \
        {params.tmpdir}/combined.viralcontigs.fa &> {log}

    mv {params.tmpdir}/combined.viralcontigs.fa.*.bt2l {params.outdir}
    
    rm -rf {params.tmpdir}/*
    """

rule mapping_bowtie2:
  name: "prokaryote.smk Bowtie2 mapping"
  input:
    R1s=lambda wildcards: expand(relpath("preprocess/samples/{sample_id}/{sample_id}_R1.fastq.gz"),
        sample_id = assemblies[wildcards.assembly_id]["sample_id"]),
    R2s=lambda wildcards: expand(relpath("preprocess/samples/{sample_id}/{sample_id}_R1.fastq.gz"),
        sample_id = assemblies[wildcards.assembly_id]["sample_id"])
    fasta=relpath("assembly/samples/{assembly_id}/output/final.contigs.fa")
    index=expand(relpath("binning/prokaryotic/samples/{assembly_id}/{assembly_id}.fa.{i}.bt2l"), i=range(1,5))
  output:
    bam=relpath("results/prokaryote/samples/{assembly_id}/{assembly_id}.bam")
  params:
    bowtie2params=configdict["bowtie2params"],
    prefix=relpath("assembly/samples/{assembly_id}/output/final.contigs.fa"),
    tmp_dir=os.path.join(tmpd, "bowtie2"),
  log: os.path.join(logdir, "bowtie2_{assembly_id}.log")
  conda: "../envs/bowtie2.yml"
  threads: 8
  shell:
    """
    rm -rf {params.tmp_dir}
    mkdir -p {params.tmp_dir}

    bowtie2 \
        -1 $(echo "{input.R1s}" | tr ' ' ',') \
        -2 $(echo "{input.R2s}" | tr ' ' ',') \
        -p {threads} \
        -x {params.prefix} \
        {params.bowtie2params}  2>{log} | samtools view -bSu - | samtools sort - -o {params.tmp_out} 2> {log}
    
    mv {params.tmp_out} {output}
    """


rule comebin:
  name: "prokaryote.smk COMEBin binning"
  input:
    bams=lambda wildcards: expand(relpath("binning/prokaryotic/samples/{sample_id}/strobealign/{sample_id}.sorted.bam"),
        sample_id = assemblies[wildcards.assembly_id]["sample_id"]),
    fasta=relpath("assembly/samples/{assembly_id}/output/final.contigs.fa")
  output:
    relpath("binning/prokaryotic/assemblies/{assembly_id}/COMEBin/comebin_res/comebin_res.tsv")
  params:
    parameters=configdict["COMEBinparams"],
    outdir=relpath("binning/prokaryotic/assemblies/{assembly_id}/COMEBin"),
    tmpdir=os.path.join(tmpd, "COMEBin/{assembly_id}")
  log: os.path.join(logdir, "COMEBin_{assembly_id}.log")
  conda: "../envs/comebin.yml"
  threads: 8
  resources:
    mem_mb=lambda wildcards, attempt, input: 20 * 10**3 * attempt
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir}/bams 
  
    cp {input.bams} {params.tmpdir}/bams

    run_comebin.sh \
        -a {input.fasta} \
        -p {params.tmpdir}/bams \
        -o {params.tmpdir} \
        -t {threads} \
        {params.parameters} 2> {log}
    
    rm -rf {params.tmpdir}/bams
    mv {params.tmpdir}/* {params.outdir}
    """

rule vamb:
  name: "prokaryote.smk VAMB binning"
  input:
    jgi=relpath("binning/prokaryotic/assemblies/{assembly_id}/MetaBAT2/depthfile.txt"),
    fasta=relpath("assembly/samples/{assembly_id}/output/final.contigs.fa")
  output:
    relpath("binning/prokaryotic/assemblies/{assembly_id}/VAMB/vae_clusters_unsplit.tsv")
  params:
    parameters=configdict["VAMBparams"],
    outdir=relpath("binning/prokaryotic/samples/{assembly_id}/VAMB/bins"),
    tmpdir=os.path.join(tmpd, "VAMB/{assembly_id}")
  log: os.path.join(logdir, "VAMB_{assembly_id}.log")
  benchmark: os.path.join(benchmarks, "MaxBin2_{assembly_id}.log")
  conda: "../envs/vamb.yml"
  threads: min(8, n_cores)
  resources:
    mem_mb=lambda wildcards, attempt, input: 20 * 10**3 * attempt
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir}

    vamb \
        --outdir {params.outdir} \
        --fasta {input.fasta} \
        --jgi {input.jgi} \
        {params.parameters} 2> {log}

    mv {params.tmpdir}/* {params.outdir}
    """
