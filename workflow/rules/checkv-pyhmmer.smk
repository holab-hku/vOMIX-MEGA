configdict = config['viral-contigident']['checkv-pyhmmer']
logdir=relpath("viralcontigident/logs")
benchmarks=relpath("viralcontigident/benchmarks")
tmpd=relpath("viralcontigident/tmp")

os.makedirs(logdir, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)


if isinstance(config['cores'], int):
 n_cores = config['cores']
else:
  console.print(Panel.fit(f"config['cores'] is not an integer: {config['cores']}, you can change the parameter in config/config.yml file", title="Error", subtitle="config['cores'] not integer"))
  sys.exit(1)

############################
# Single-Sample Processing #
############################

if config['fasta']!="":

  fastap = config['fasta']
  _, extension = os.path.splitext(fastap)

  console.print(f"\n[dim]The config['fasta'] parameter is not empty, using '{fastap}' as input.")

  if extension.lower() not in ['.fa', '.fasta', '.fna']:
    console.print(Panel.fit("File path does not end with .fa, .fasta, or .fna", title = "Error", subtitle="Input not fasta file"))
    sys.exit(1)

  cwd = os.getcwd()
  fasta_path = os.path.join(cwd, fastap)

  if not os.path.exists(fastap):
    console.print(Panel.fit("The fasta file path provided does not exist.", title="Error", subtitle="Contig File Path"))
    sys.exit(1)

  outdir_p = os.path.join(cwd, relpath("viralcontigident/output/checkv/"))
  console.print(f"[dim]Output file will be written to the '{outdir_p}' directory.\n")

  try:
    if len(os.listdir(outdir_p)) > 0:
      console.print(Panel.fit(f"Output directory '{outdir_p}' already exists and is not empty.", title = "Warning", subtitle="Output Directory Not Empty"))
  except Exception:
    pass

  sample_id = os.path.splitext(os.path.basename(fastap))[0]

else:
  fasta_path = relpath("viralcontigident/output/derep/combined.viralcontigs.derep.fa")

### MASTER RULE

rule done_log:
  name: "checkv-pyhmmer.py Done. removing tmp files"
  localrule: True
  input:
    relpath("viralcontigident/output/checkv/viruses.fna"),
    relpath("viralcontigident/output/checkv/proviruses.fna"),
    relpath("viralcontigident/output/checkv/quality_summary.tsv")
  output:
    os.path.join(logdir, "checkv-done.log")
  shell: "touch {output}"


### RULES

rule checkv_prodigalgv:
  name: "checkv-pyhmmer.smk CheckV run prodigal-gv"
  input: fasta_path
  output:
    relpath("viralcontigident/output/checkv/tmp/proteins.faa")
  params:
    script="workflow/software/prodigal-gv/parallel-prodigal-gv.py",
    outdir=relpath("viralcontigident/output/checkv/tmp"),
    tmpdir=os.path.join(tmpd, "checkv/prodigal-gv")
  log: os.path.join(logdir, "checkv_prodigal-gv.log")
  conda: "../envs/prodigal-gv.yml"
  threads: min(64, n_cores)
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
    rm -rf {params.tmpdir}/*
    """

rule checkv_pyhmmer:
  name: "checkv-pyhmmer.smk CheckV PyHMMER hmmsearch"
  input:
    faa=relpath("viralcontigident/output/checkv/tmp/proteins.faa"), 
    db="workflow/database/checkv/hmm_db/checkv_hmms/{index}.hmm"
  output:
    relpath("viralcontigident/output/checkv/tmp/hmmsearch/{index}.hmmout")
  params:
    script="workflow/scripts/taxonomy/pyhmmer_wrapper.py",
    outdir=relpath("viralcontigident/output/checkv/tmp/hmmsearch"),
    tmpdir=os.path.join(tmpd, "checkv/hmmsearch/{index}"), 
    ecutoff=10.0
  log : os.path.join(logdir, "checkv_hmmsearch_{index}.log")
  benchmark: os.path.join(benchmarks, "checkv_hmmsearch_{index}.log")
  conda: "../envs/pyhmmer.yml"
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt: attempt * 16 * 10**3
  shell:
    """
    rm -rf {params.tmpdir}
    mkdir -p {params.tmpdir} {params.outdir}

    python {params.script} \
         --proteins {input.faa} \
         --hmmdb {input.db} \
         --cores {threads} \
         --e_value {params.ecutoff} \
         --tblout {params.tmpdir}/tmp.hmmout 2> {log}

    mv {params.tmpdir}/tmp.hmmout {output}
    rm -rf {params.tmpdir}
    """

rule checkv_hmm_merge:
  name: "checkv-pyhmmer.smk CheckV hmmsearch merge"
  localrule: True
  input:
    expand(relpath("viralcontigident/output/checkv/tmp/hmmsearch/{index}.hmmout"), index = range(1, 81))
  output:
    relpath("viralcontigident/output/checkv/tmp/hmmsearch.txt")
  shell:
    """
    cat {input} > {output}

    """

rule checkv_hmmer_checkpoint:
  name: "checkv-pyhmmer.smk CheckV hmmsearch checkpoint"
  localrule: True
  input:
    relpath("viralcontigident/output/checkv/tmp/hmmsearch.txt")
  output:
    relpath("viralcontigident/output/checkv/tmp/hmmsearch_checkpoint")
  shell:
    """
    touch {output}
    """

# This rule currently does not operate in tmpdir so that it matches the checkpoint
# Fix this later please thanks!

rule checkv:
  name: "checkv-pyhmmer.smk CheckV dereplicated contigs"
  input:
    checkpoint=relpath("viralcontigident/output/checkv/tmp/hmmsearch_checkpoint"),
    fna=fasta_path
  output:
    relpath("viralcontigident/output/checkv/viruses.fna"),
    relpath("viralcontigident/output/checkv/proviruses.fna"),
    relpath("viralcontigident/output/checkv/quality_summary.tsv")
  params:
    checkvparams= configdict['checkvparams'],
    outdir=relpath("viralcontigident/output/checkv"),
    tmpdir=os.path.join(tmpd, "checkv"),
    dbdir="workflow/database/checkv"
  log: os.path.join(logdir, "checkv.log")
  benchmark: os.path.join(benchmarks, "checkv.log")
  threads: 1
  resources:
    mem_mb=lambda wildcards, attempt, input: attempt * 72 * 10**3
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
        {params.checkvparams} 2> {log}

    rm -rf {params.tmpdir}
    """
