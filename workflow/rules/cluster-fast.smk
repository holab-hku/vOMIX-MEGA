import math
import os

wildcard_constraints:
    layer=r"\d+",
    chunk=r"\d+"

logdir = relpath("identify/viral/logs")
tmpd = relpath("identify/viral/tmp")
benchmarks=relpath("identify/viral/benchmarks")

os.makedirs(logdir, exist_ok=True)
os.makedirs(tmpd, exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)

clustering_iter = config["cluster-iter"] # nLayers
n_chunks_layer_1 =  2 ** (clustering_iter - 1) # nCluster Chunks for Layer 1

#### Helper function to choose inputs for recursive clustering
def get_iter_inputs(wildcards):
    """Calculates which files to merge based on the current layer and chunk."""
    layer = int(wildcards.layer)
    chunk = int(wildcards.chunk)
    
    if layer == 1:
        # layer 1 - just grabs the raw split pieces
        return os.path.join(relpath("identify/viral/output/derep/cluster-splits"), f"chunk_{chunk}.fa")
    else:
        # Layer > 1 - grab the two specific clustered chunks from the previous layer
        prev_layer = layer - 1
        child1 = chunk * 2
        child2 = chunk * 2 + 1
        return [
            os.path.join(relpath("identify/viral/output/derep/cluster-layers"), f"layer_{prev_layer}/chunk_{child1}.fa"), 
            os.path.join(relpath("identify/viral/output/derep/cluster-layers"), f"layer_{prev_layer}/chunk_{child2}.fa")
        ]
 
        
    
### Read single fasta file if input
if config['fasta'] != "" and config["module"] == "cluster-fast":
    fastap = readfasta(config["fasta"])
    sample_id = config["sample-name"]
    assembly_ids = [sample_id]
else:
    fastap = relpath("identify/viral/intermediate/scores/combined.viralcontigs.fa")
    sample_id = "combined.viralcontigs"
    assembly_ids = [sample_id]


### MASTER RULE
rule done_log:
    name: "clustering.smk Done. removing tmp files"
    localrule: True
    input:
        relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa.clstr"), 
        os.path.join(relpath("identify/viral/output/derep/cluster-layers"), f"layer_{clustering_iter}", "chunk_0.fa"),
    output:
        os.path.join(logdir, "clustering-done.log")
    shell: "touch {output}"


### RULES 
rule split_input:
    name: "clustering.smk split input fasta"
    localrule: True
    input:
        fastap
    output:
        expand(os.path.join(relpath("identify/viral/output/derep/cluster-splits"), f"chunk_{{chunk}}.fa"), chunk=range(n_chunks_layer_1))
    params:
        pieces = n_chunks_layer_1,
        outdir = relpath("identify/viral/output/derep/cluster-splits"),
        tmpdir = os.path.join(tmpd, "cluster-splits")
    log: os.path.join(logdir, "split_input.log")
    benchmark: os.path.join(benchmarks, "split_input.log")
    conda: "../envs/seqkit-biopython.yml"
    threads: 1
    shell:
        """
        rm -rf {params.outdir} {params.tmpdir}
        mkdir -p {params.outdir} {params.tmpdir}
            
        seqkit split2 {input} -p {params.pieces} -O {params.tmpdir}/
            
        counter=0
        shopt -s nullglob
        set -- {params.tmpdir}/*.fa {params.tmpdir}/*.fna {params.tmpdir}/*.fasta
        for file in "$@"; do
            mv "$file" "{params.outdir}/chunk_${{counter}}.fa"
            counter=$((counter+1))
        done       
        shopt -u nullglob 
        """


# MEGABLAST (clustering-fast=True)
if config["clustering-fast"]:
    rule megablast_prep_input:
        name: "clustering.smk prepare split input"
        localrule: True
        input:
            get_iter_inputs
        output:
            os.path.join(relpath("identify/viral/output/derep/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "input.fa")
        params:
            tmpdir = lambda wildcards, output: os.path.dirname(output[0])
        log: os.path.join(logdir, "megablast_prep_layer_{layer}_chunk_{chunk}.log")
        threads: 1
        shell:
            """
            mkdir -p {params.tmpdir}
            
            if [ "{wildcards.layer}" -eq "1" ]; then
              # layer = 1: input contains exactly one file so output is identical to input
              cp {input} {output} 2> {log}
            else
              # layer > 1: input automatically expands to multiple "chunk_0.fa chunk_1.fa ..." files and is pooled
              cat {input} > {output} 2> {log}
            fi
            """

    rule megablast_makeblastdb:
        name: "clustering.smk make blast db"
        input:
            os.path.join(relpath("identify/viral/output/derep/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "input.fa")
        output:
            expand(os.path.join(relpath("identify/viral/output/derep/cluster-layers"), "layer_{{layer}}", "chunk_{{chunk}}", "db.{suffix}"), suffix=["ntf", "ndb"])
        params:
            tmpdir = lambda wildcards, input: os.path.dirname(input[0]),
            dbtype = 'nucl'
        log: os.path.join(logdir, "megablast_makeblastdb_layer_{layer}_chunk_{chunk}.log")
        benchmark: os.path.join(benchmarks, "megablast_makeblastdb_layer_{layer}_chunk_{chunk}.log")
        conda: "../envs/checkv.yml"
        threads: 1
        resources:
            mem_mb = lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
        shell:
            """
            makeblastdb -in {input} -dbtype {params.dbtype} -out {params.tmpdir}/db &> {log}
            """

    rule megablast_run:
        name: "clustering.smk megablast run"
        input:
            fasta = os.path.join(relpath("identify/viral/output/derep/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "input.fa"),
            dbcheckpoints = expand(os.path.join(relpath("identify/viral/output/derep/cluster-layers"), "layer_{{layer}}", "chunk_{{chunk}}", "db.{suffix}"), suffix=["ntf", "ndb"])
        output:
            os.path.join(relpath("identify/viral/output/derep/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "blast_out.csv")
        params:
            db = lambda wildcards, input: os.path.join(os.path.dirname(input.fasta), "db"),
            outfmt = "'6 std qlen slen'",
            maxtargetseqs = 10000
        log: os.path.join(logdir, "megablast_run_layer_{layer}_chunk_{chunk}.log")
        benchmark: os.path.join(benchmarks, "megablast_run_layer_{layer}_chunk_{chunk}.log")
        conda: "../envs/checkv.yml"
        threads: 64
        resources:
          mem_mb = lambda wildcards, threads, attempt: int(attempt * (threads / 64) * 72 * 10**3)
        shell:
            """
            blastn -query {input.fasta} \
                -db {params.db} \
                -outfmt {params.outfmt} \
                -max_target_seqs {params.maxtargetseqs} \
                -out {output} \
                -num_threads {threads} &> {log}
            """
  
    rule megablast_anicalc:
        name: "clustering.smk calculate ani"
        input:
            os.path.join(relpath("identify/viral/output/derep/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "blast_out.csv")
        output:
            os.path.join(relpath("identify/viral/output/derep/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "ani.tsv")
        params:
            script = "workflow/scripts/clust_anicalc.py"
        log: os.path.join(logdir, "megablast_anicalc_layer_{layer}_chunk_{chunk}.log")
        benchmark: os.path.join(benchmarks, "megablast_anicalc_layer_{layer}_chunk_{chunk}.log")
        conda: "../envs/checkv.yml"
        threads: 1
        resources:
            mem_mb = lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
        shell:
            """
            python {params.script} \
                -i {input} \
                -o {output} &> {log}
            """

    rule megablast_aniclust:
        name: "clustering.smk cluster by ani"
        input:
            fa = os.path.join(relpath("identify/viral/output/derep/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "input.fa"),
            ani = os.path.join(relpath("identify/viral/output/derep/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "ani.tsv")
        output:
            clstr = os.path.join(relpath("identify/viral/output/derep/cluster-layers"), "layer_{layer}", "chunk_{chunk}.fa.clstr"),
            reps = os.path.join(relpath("identify/viral/output/derep/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "cluster_representatives.txt")
        params:
            script = "workflow/scripts/clust_ani.py",
            minani = config["vOTU-ani"],
            targetcov = config["vOTU-targetcov"],
            querycov = config["vOTU-querycov"],
            outdir = lambda wildcards, output: os.path.dirname(output.clstr)
        log: os.path.join(logdir, "megablast_aniclust_layer_{layer}_chunk_{chunk}.log")
        benchmark: os.path.join(benchmarks, "megablast_aniclust_layer_{layer}_chunk_{chunk}.log")
        conda: "../envs/checkv.yml"
        threads: 1
        resources:
            mem_mb = lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
        shell:
            """
            mkdir -p {params.outdir}
            python {params.script} \
                --fna {input.fa} \
                --ani {input.ani} \
                --out {output.clstr} \
                --min_ani {params.minani} \
                --min_tcov {params.targetcov} \
                --min_qcov {params.querycov} &> {log}
                
            cut -f1 {output.clstr} > {output.reps}
            """

    rule megablast_filtercontigs:
        name: "clustering.smk filter dereplicated viral contigs"
        input:
            fna = os.path.join(relpath("identify/viral/output/derep/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "input.fa"),
            reps = os.path.join(relpath("identify/viral/output/derep/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "cluster_representatives.txt")
        output:
            fa = os.path.join(relpath("identify/viral/output/derep/cluster-layers"), "layer_{layer}", "chunk_{chunk}.fa")
        params:
            tmpdir = lambda wildcards, input: os.path.dirname(input.fna),
            outdir = lambda wildcards, output: os.path.dirname(output.fa)
        log: os.path.join(logdir, "megablast_filtercontigs_layer_{layer}_chunk_{chunk}.log")
        conda: "../envs/seqkit-biopython.yml" 
        threads: 1
        resources:
            mem_mb = lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
        shell:
            """
            mkdir -p {params.outdir}
            seqkit grep {input.fna} -f {input.reps} > {output.fa} 2> {log}
            
            rm -rf {params.tmpdir}
            """


# CD-HIT (clustering-fast=False)
else:        
    rule cdhit_recursive_cluster:
        name: "clustering.smk CD-HIT recursive clustering"
        input:
            get_iter_inputs
        output:
            fa = os.path.join(relpath("identify/viral/output/derep/cluster-layers"), f"layer_{{layer}}", f"chunk_{{chunk}}.fa"),
            clstr = os.path.join(relpath("identify/viral/output/derep/cluster-layers"), f"layer_{{layer}}", f"chunk_{{chunk}}.fa.clstr")
        params:
            cdhitparams = config['cdhit-params'],
            outdir = os.path.join(relpath("identify/viral/output/derep/cluster-layers"), f"layer_{{layer}}"), 
            tmpdir = os.path.join(relpath("identify/viral/output/derep/cluster-layers"), f"layer_{{layer}}", f"chunk_{{chunk}}")
        threads: 64
        resources:
            mem_mb = lambda wildcards, threads, attempt: int(attempt * (threads / 64) * 72 * 10**3)
        conda: "../envs/cd-hit.yml"
        log: os.path.join(logdir, "cdhit_layer_{layer}_chunk_{chunk}.log")
        benchmark: os.path.join(benchmarks, "cdhit_layer_{layer}_chunk_{chunk}.log")
        shell:
            """
            rm -rf {params.tmpdir}
            mkdir -p {params.tmpdir} {params.outdir}
            
            chunk_name=$(basename "{output.fa}" | sed 's/\\.[^.]*$//')
            
            if [ "{wildcards.layer}" -eq "1" ]; then
    
                # layer = 1: input contains exactly one file
    
                echo "[$(date '+%Y-%m-%d %H:%M:%S')] PROCESSING: Chunk '${{chunk_name}}' | LAYER: {wildcards.layer} | CONDITION: Initial CD-HIT Clustering" >> {log}
                cd-hit -i {input} -o {params.tmpdir}/tmp_output -T {threads} {params.cdhitparams} >> {log} 2>&1
    
                mv {params.tmpdir}/tmp_output {output.fa}
                mv {params.tmpdir}/tmp_output.clstr {output.clstr}
                rm -rf {params.tmpdir}
    
            else
    
                # layer > 1: inputs automatically expands to "chunk_0.fa chunk_1.fa ..." 
    
                cat {input} > {params.tmpdir}/tmp_input
                echo "[$(date '+%Y-%m-%d %H:%M:%S')] PROCESSING: Chunk '${{chunk_name}}' | LAYER: {wildcards.layer} | CONDITION: Pooling multiple inputs & Clustered Processing" >> {log}
                cd-hit -i {params.tmpdir}/tmp_input -o {params.tmpdir}/tmp_output -T {threads} {params.cdhitparams} >> {log} 2>&1
                
                mv {params.tmpdir}/tmp_output {output.fa}
                mv {params.tmpdir}/tmp_output.clstr {output.clstr}
                rm -rf {params.tmpdir}
            fi        
            """

rule derep_finalize:
    name: "clustering.smk finalize dereplicated dataset"
    localrule: True
    input:
        fa = os.path.join(relpath("identify/viral/output/derep/cluster-layers"), f"layer_{clustering_iter}", "chunk_0.fa"),
        clstr = os.path.join(relpath("identify/viral/output/derep/cluster-layers"), f"layer_{clustering_iter}", "chunk_0.fa.clstr")
    output:
        fa = relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa"),
        clstr = relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa.clstr"), 
    log: os.path.join(logdir, "cdhit_final.log")
    benchmark: 
        os.path.join(benchmarks, "identify/viral_cdhit.log")
    shell:
        """
        cp {input.fa} {output.fa}
        cp {input.clstr} {output.clstr}
        """