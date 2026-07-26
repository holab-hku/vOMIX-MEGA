logdir=relpath("identify/viral/logs")
benchmarks=relpath("identify/viral/benchmarks")
tmpd=relpath("identify/viral/tmp")

os.makedirs(logdir, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)

n_cores = config['max-cores']
split_part = list(range(1, config["checkv-splits"] + 2))
split_part_ids = [f"{i:03d}" for i in range(1, config["checkv-splits"] + 2)] # matches seqkit
parts=config["checkv-splits"] + 1

### Read single fasta file if input
if config['fasta'] != "" and config["module"] == "checkv-pyhmmer":
  fastap = readfasta(config['fasta'])
  sample_id = config["sample-name"]
  assembly_ids = [sample_id]
else:
  fastap = relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa")
  sample_id = "combined.viralcontigs.derep"

### MASTER RULE
rule done_log:
  name: "checkv-pyhmmer.py Done. removing tmp files"
  localrule: True
  input:
    expand(relpath(f"identify/viral/tmp/splits/{sample_id}.part_{{part}}.fa"), part=split_part_ids),
    relpath("identify/viral/output/checkv/viruses.fna"),
    relpath("identify/viral/output/checkv/proviruses.fna"),
    relpath("identify/viral/output/checkv/quality_summary.tsv")
  output:
    os.path.join(logdir, "checkv-done.log")
  shell: "touch {output}"


### RULES
rule split_contigs:
  name: "checkv-pyhmmer.smk split dereplicated contigs"
  input: fastap
  output:
    expand(relpath(f"identify/viral/tmp/splits/{sample_id}.part_{{part}}.fa"), part=split_part_ids)
  params:
    parts=parts,
    sample_id=sample_id, # Pass sample_id explicitly to params so the shell can read it
    tmpdir=os.path.join(tmpd, "checkv/splits/"),
    outdir=relpath("identify/viral/tmp/splits/"),
  log: os.path.join(logdir, "checkv_splitcontig.log")
  conda: "../envs/seqkit-biopython.yml"
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input: 4 * 10**3 * attempt
  shell:
    """
    rm -rf {params.tmpdir} {params.outdir}
    mkdir -p {params.tmpdir} {params.outdir}

    # count how many sequences (lines starting with '>') are in the input file
    seq_count=$(grep -c "^>" "{input}")

    # check if we have fewer sequences than requested parts
    if [ "$seq_count" -lt "{params.parts}" ]; then
        echo "ERROR: Input file has only $seq_count sequences, but {params.parts} parts were requested!" >> {log}
        exit 1
    fi

    # Run seqkit split2 (overwrite/start the log)
    seqkit split2 \
      --by-part {params.parts} \
      --threads {threads} {input} \
      --out-dir {params.tmpdir} 2> {log}

    # Clean and rename files to force a strict '.fa' extension
    for file in "{params.tmpdir}"/*; do
        if [ -f "$file" ]; then
            base_name=$(basename "$file")
            part_suffix=$(echo "$base_name" | grep -o "part_[0-9]\\+")
            mv "$file" "{params.outdir}/{params.sample_id}.${{part_suffix}}.fa"
        fi
    done

    rm -rf {params.tmpdir}    
    """    
    
# A) Original CheckV
if config["checkv-original"]:
  rule checkv:
    name: "checkv.smk CheckV split contigs"
    input:
      fna=relpath(f"identify/viral/tmp/splits/{sample_id}.part_{{part}}.fa"), 
      db=expand(os.path.join(config['checkv-db'], "hmm_db/checkv_hmms/{index}.hmm"), index=range(1, 81))
    output:
      relpath("identify/viral/output/splits/split-{part}/checkv/viruses.fna"),
      relpath("identify/viral/output/splits/split-{part}/checkv/proviruses.fna"),
      relpath("identify/viral/output/splits/split-{part}/checkv/quality_summary.tsv")
    params:
      parameters= config['checkv-params'],
      outdir=relpath("identify/viral/output/splits/split-{part}/checkv"),
      dbdir=config["checkv-db"], 
      tmpdir=os.path.join(tmpd, "checkv/splits/split-{part}"),
    log: os.path.join(logdir, "checkv_split_{part}.log")
    benchmark: os.path.join(benchmarks, "checkv_split_{part}.log")
    threads: max(1, round(64 / parts))
    resources:
      mem_mb=lambda wildcards, attempt, input: attempt * max(1, round(72 / parts)) * 10**3
    conda: "../envs/checkv.yml"
    shell:
      """
      rm -rf {params.tmpdir}
      mkdir -p {params.tmpdir} {params.outdir}

      checkv end_to_end \
          {input.fna} \
          {params.outdir} \
          -d {params.dbdir} \
          -t {threads} \
          {params.parameters} 2> {log}

      rm -rf {params.tmpdir}
      """

  rule checkv_merge:
    name: "checkv.smk CheckV merge split results"
    localrule: True
    input:
      virus=expand(relpath("identify/viral/output/splits/split-{part}/checkv/viruses.fna"), part=split_part_ids),
      provirus=expand(relpath("identify/viral/output/splits/split-{part}/checkv/proviruses.fna"), part=split_part_ids),
      summary=expand(relpath("identify/viral/output/splits/split-{part}/checkv/quality_summary.tsv"), part=split_part_ids)
    output:
      virus=relpath("identify/viral/output/checkv/viruses.fna"),
      provirus=relpath("identify/viral/output/checkv/proviruses.fna"),
      summary=relpath("identify/viral/output/checkv/quality_summary.tsv")
    params:
      script="workflow/scripts/tables_row_bind.py",
      outdir=relpath("identify/viral/output/checkv"), 
      tmpdir=os.path.join(tmpd, "checkv/merge")
    log: os.path.join(logdir, "checkv_merge.log")
    benchmark: os.path.join(benchmarks, "checkv_merge.log")
    conda: "../envs/seqkit-biopython.yml"
    shell:
      """
      rm -rf {params.tmpdir} {params.outdir}
      mkdir -p {params.tmpdir} {params.outdir}

      cat {input.virus} > {params.tmpdir}/tmp.viruses.fna
      cat {input.provirus} > {params.tmpdir}/tmp.proviruses.fna
      python {params.script} --inputs {input.summary} --output {params.tmpdir}/tmp.quality_summary.tsv

      mv {params.tmpdir}/tmp.viruses.fna {output.virus}
      mv {params.tmpdir}/tmp.proviruses.fna {output.provirus}
      mv {params.tmpdir}/tmp.quality_summary.tsv {output.summary}
      rm -rf {params.tmpdir}
      """

# B) CheckV-PyHMMER
else:
  rule checkv_prodigalgv:
    name: "checkv-pyhmmer.smk CheckV run prodigal-gv"
    input: 
      relpath(f"identify/viral/tmp/splits/{sample_id}.part_{{part}}.fa"), 
    output:
      relpath("identify/viral/output/splits/split-{part}/checkv/tmp/proteins.faa")
    params:
      script="workflow/scripts/parallel_prodigal_gv.py",
      outdir=relpath("identify/viral/output/splits/split-{part}/checkv/tmp"),
      tmpdir=os.path.join(tmpd, "checkv/splits/split-{part}/prodigal-gv")
    log: os.path.join(logdir, "checkv_split_{part}_prodigal-gv.log")
    benchmark: os.path.join(benchmarks, "checkv_split_{part}_prodigal-gv.log")
    conda: "../envs/prodigal-gv.yml"
    threads: max(1, round(128 / parts))
    resources:
      mem_mb=lambda wildcards, attempt: attempt * max(1, round(72 / parts)) * 10**3
    shell:
      """
      rm -rf {params.tmpdir} {params.outdir}
      mkdir -p {params.tmpdir} {params.outdir}

      python {params.script} \
          -i {input} \
          -a {params.tmpdir}/tmp.faa \
          -t {threads} &> {log}

      mv {params.tmpdir}/tmp.faa {output}
      rm -rf {params.tmpdir}/*
      """

  rule checkv_pyhmmer:
    name: "checkv-pyhmmer.smk CheckV PyHMMER hmmsearch"
    localrule: False
    input:
      faa=relpath("identify/viral/output/splits/split-{part}/checkv/tmp/proteins.faa"), 
      db=os.path.join(config["checkv-db"], "hmm_db/checkv_hmms/{index}.hmm"), 
      dbfull=expand(os.path.join(config["checkv-db"], "hmm_db/checkv_hmms/{index}.hmm"), index=range(1, 81))
    output:
      relpath("identify/viral/output/splits/split-{part}/checkv/tmp/hmmsearch/{index}.hmmout")
    params:
      script="workflow/scripts/pyhmmer_wrapper.py",
      outdir=relpath("identify/viral/output/splits/split-{part}/checkv/tmp/hmmsearch"),
      tmpdir=os.path.join(tmpd, "checkv/splits/split-{part}/hmmsearch/{index}"), 
      ecutoff=10.0
    log : os.path.join(logdir, "checkv_hmmsearch_split_{part}_{index}.log")
    benchmark: os.path.join(benchmarks, "checkv_hmmsearch_split_{part}_{index}.log")
    conda: "../envs/pyhmmer.yml"
    threads: 1
    resources:
      mem_mb=lambda wildcards, attempt: attempt * 4 * 10**3
    shell:
      """
      rm -rf {params.tmpdir}
      mkdir -p {params.tmpdir} {params.outdir}

      python {params.script} \
           --proteins {input.faa} \
           --hmmdb {input.db} \
           --cores {threads} \
           --e_value {params.ecutoff} \
           --tblout {params.tmpdir}/tmp.hmmout &> {log}

      mv {params.tmpdir}/tmp.hmmout {output}
      rm -rf {params.tmpdir}
      """

  rule checkv_hmm_merge:
    name: "checkv-pyhmmer.smk CheckV hmmsearch merge"
    localrule: True
    input:
      expand(relpath("identify/viral/output/splits/split-{{part}}/checkv/tmp/hmmsearch/{index}.hmmout"), index = range(1, 81))
    output:
      relpath("identify/viral/output/splits/split-{part}/checkv/tmp/hmmsearch.txt")
    shell:
      """
      cat {input} > {output}
      """

  rule checkv_hmmer_checkpoint:
    name: "checkv-pyhmmer.smk CheckV hmmsearch checkpoint"
    localrule: True
    input:
      relpath("identify/viral/output/splits/split-{part}/checkv/tmp/hmmsearch.txt")
    output:
      relpath("identify/viral/output/splits/split-{part}/checkv/tmp/hmmsearch_checkpoint")
    shell:
      """
      touch {output}
      """

  # This rule currently does not operate in tmpdir so that it matches the checkpoint
  # Fix this later please thanks!
  rule checkv:
    name: "checkv-pyhmmer.smk CheckV dereplicated contigs"
    input:
      checkpoint=relpath("identify/viral/output/splits/split-{part}/checkv/tmp/hmmsearch_checkpoint"),
      fna=relpath(f"identify/viral/tmp/splits/{sample_id}.part_{{part}}.fa"),
      db=expand(os.path.join(config['checkv-db'], "hmm_db/checkv_hmms/{index}.hmm"), index=range(1, 81))
    output:
      relpath("identify/viral/output/splits/split-{part}/checkv/viruses.fna"),
      relpath("identify/viral/output/splits/split-{part}/checkv/proviruses.fna"),
      relpath("identify/viral/output/splits/split-{part}/checkv/quality_summary.tsv")
    params:
      params= config['checkv-params'],
      dbdir=config["checkv-db"],
      outdir=relpath("identify/viral/output/splits/split-{part}/checkv"),
      tmpdir=os.path.join(tmpd, "checkv/splits/split-{part}"),
    log: os.path.join(logdir, "checkv_split-{part}.log")
    benchmark: os.path.join(benchmarks, "checkv_split-{part}.log")
    threads: 1
    resources:
      mem_mb=lambda wildcards, attempt, input: attempt * 72/parts * 10**3
    conda: "../envs/checkv.yml"
    shell:
      """
      rm -rf {params.tmpdir}
      mkdir -p {params.tmpdir} {params.outdir}

      checkv end_to_end \
          {input.fna} \
          {params.outdir} \
          -d {params.dbdir} \
          -t {threads} \
          {params.params} 2> {log}

      rm -rf {params.tmpdir}
      """

  rule checkv_merge:
    name: "checkv-pyhmmer.smk CheckV merge split results"
    localrule: True
    input:
      virus=expand(relpath("identify/viral/output/splits/split-{part}/checkv/viruses.fna"), part=split_part_ids),
      provirus=expand(relpath("identify/viral/output/splits/split-{part}/checkv/proviruses.fna"), part=split_part_ids),
      summary=expand(relpath("identify/viral/output/splits/split-{part}/checkv/quality_summary.tsv"), part=split_part_ids)
    output:
      virus=relpath("identify/viral/output/checkv/viruses.fna"),
      provirus=relpath("identify/viral/output/checkv/proviruses.fna"),
      summary=relpath("identify/viral/output/checkv/quality_summary.tsv")
    params:
      script="workflow/scripts/tables_row_bind.py",
      outdir=relpath("identify/viral/output/checkv/"), 
      tmpdir=os.path.join(tmpd, "checkv/merge")
    log: os.path.join(logdir, "checkv_merge.log")
    benchmark: os.path.join(benchmarks, "checkv_merge.log")
    threads: 1
    resources:
      mem_mb=lambda wildcards, attempt, input: attempt * 4 * 10**3
    conda: "../envs/seqkit-biopython.yml"
    shell:
      """
      rm -rf {params.tmpdir} {params.outdir}
      mkdir -p {params.tmpdir} {params.outdir}

      cat {input.virus} > {params.tmpdir}/tmp.viruses.fna
      cat {input.provirus} > {params.tmpdir}/tmp.proviruses.fna
      python {params.script} --inputs {input.summary} --output {params.tmpdir}/tmp.quality_summary.tsv

      mv {params.tmpdir}/tmp.viruses.fna {output.virus}
      mv {params.tmpdir}/tmp.proviruses.fna {output.provirus}
      mv {params.tmpdir}/tmp.quality_summary.tsv {output.summary}
      rm -rf {params.tmpdir}
      """
