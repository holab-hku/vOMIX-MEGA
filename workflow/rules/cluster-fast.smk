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

cluster_method_input = config.get("cluster-method")  # keep original config value
cluster_iter = config.get("cluster-iter")            # nLayers
n_chunks_layer_1 =  2 ** (cluster_iter - 1)          # nCluster Chunks for Layer 1

### Read single fasta file if input
if config.get("fasta", "") != "" and config.get("module", "") == "cluster-fast":
    fastap = readfasta(config.get("fasta", ""))
    sample_id = config.get("sample-name", "")
    assembly_ids = [sample_id]
else:
    fastap = relpath("identify/viral/intermediate/scores/combined.viralcontigs.fa")
    sample_id = "combined.viralcontigs"
    assembly_ids = [sample_id]


# ----- MASTER RULE -----
if cluster_method_input == "all":
    final_targets = [
        relpath("identify/viral/output/derep/checkv-megablast/combined.viralcontigs.derep.fa"),
        relpath("identify/viral/output/derep/checkv-megablast/combined.viralcontigs.derep.fa.clstr"),
        relpath("identify/viral/output/derep/cd-hit-est/combined.viralcontigs.derep.fa"),
        relpath("identify/viral/output/derep/cd-hit-est/combined.viralcontigs.derep.fa.clstr"),
        relpath("identify/viral/output/derep/vclust/combined.viralcontigs.derep.fa"),
        relpath("identify/viral/output/derep/vclust/combined.viralcontigs.derep.fa.clstr"),
        relpath("identify/viral/output/derep/linclust/combined.viralcontigs.derep.fa"),
        relpath("identify/viral/output/derep/linclust/combined.viralcontigs.derep.fa.clstr"),
        relpath("identify/viral/output/derep/dnaclust/combined.viralcontigs.derep.fa"),
        relpath("identify/viral/output/derep/dnaclust/combined.viralcontigs.derep.fa.clstr"),
        relpath("identify/viral/output/derep/vsearch/combined.viralcontigs.derep.fa"),
        relpath("identify/viral/output/derep/vsearch/combined.viralcontigs.derep.fa.clstr"),
        relpath("identify/viral/output/derep/viridic/combined.viralcontigs.derep.fa"),
        relpath("identify/viral/output/derep/viridic/combined.viralcontigs.derep.fa.clstr"),
    ]
else:
    final_targets = [
        relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa"),
        relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa.clstr"),
    ]

rule done_log:
    name: "clustering.smk Done. removing tmp files"
    localrule: True
    input:
        final_targets
    output:
        os.path.join(logdir, "clustering-done.log")
    shell: "touch {output}"


rule split_input:
    name: "clustering.smk split input fasta"
    localrule: True
    input:
        fastap
    output:
        expand(os.path.join(relpath("identify/viral/tmp/derep/cluster-splits"), f"chunk_{{chunk}}.fa"), chunk=range(n_chunks_layer_1))
    params:
        pieces = n_chunks_layer_1,
        outdir = relpath("identify/viral/tmp/derep/cluster-splits"),
        tmpdir = os.path.join(tmpd, "derep", "cluster-splits"), 
        seed = config.get("seed", 42)
    log: os.path.join(logdir, "split_input.log")
    benchmark: os.path.join(benchmarks, "split_input.log")
    conda: "../envs/seqkit-biopython.yml"
    threads: 1
    shell:
        """
        rm -rf {params.outdir} {params.tmpdir}
        mkdir -p {params.outdir} {params.tmpdir}
            
        seqkit split2 {input} -p {params.pieces} -O {params.tmpdir}/ -s {params.seed}
            
        counter=0
        shopt -s nullglob
        set -- {params.tmpdir}/*.fa {params.tmpdir}/*.fna {params.tmpdir}/*.fasta
        for file in "$@"; do
            mv "$file" "{params.outdir}/chunk_${{counter}}.fa"
            counter=$((counter+1))
        done       
        shopt -u nullglob 
        """


# ----- CheckV-MEGABLAST -----
if cluster_method_input in ["checkv-megablast", "all"]:
    cluster_method = "checkv-megablast"
    _method = cluster_method  # capture value for closure

    def get_iter_inputs(wildcards, method=_method):
        layer = int(wildcards.layer)
        chunk = int(wildcards.chunk)
        if layer == 1:
            return os.path.join(relpath("identify/viral/tmp/derep/cluster-splits"), f"chunk_{chunk}.fa")
        else:
            prev_layer = layer - 1
            child1 = chunk * 2
            child2 = chunk * 2 + 1
            return [
                os.path.join(relpath(f"identify/viral/tmp/derep/{method}/cluster-layers"), f"layer_{prev_layer}", f"chunk_{child1}.fa"),
                os.path.join(relpath(f"identify/viral/tmp/derep/{method}/cluster-layers"), f"layer_{prev_layer}", f"chunk_{child2}.fa")
            ]

    rule checkv_megablast_prep_input:
        name: "clustering.smk CheckV-MEGABLAST prepare split input"
        localrule: True
        input:
            get_iter_inputs
        output:
            os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "input.fa")
        params:
            tmpdir = lambda wildcards, output: os.path.dirname(output[0])
        log: os.path.join(logdir, "megablast_prep_layer_{layer}_chunk_{chunk}.log")
        threads: 1
        shell:
            """
            mkdir -p {params.tmpdir}
            
            if [ "{wildcards.layer}" -eq "1" ]; then
              cp {input} {output} 2> {log}
            else
              cat {input} > {output} 2> {log}
            fi
            """

    rule checkv_megablast_makeblastdb:
        name: "clustering.smk CheckV-MEGABLAST make blast db"
        input:
            os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "input.fa")
        output:
            expand(os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{{layer}}", "chunk_{{chunk}}", "db.{suffix}"), suffix=["ntf", "ndb"])
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

    rule checkv_megablast_run:
        name: "clustering.smk CheckV-MEGABLAST run"
        input:
            fasta = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "input.fa"),
            dbcheckpoints = expand(os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{{layer}}", "chunk_{{chunk}}", "db.{suffix}"), suffix=["ntf", "ndb"])
        output:
            os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "blast_out.csv")
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
  
    rule checkv_megablast_anicalc:
        name: "clustering.smk CheckV-MEGABLAST calculate ani"
        input:
            os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "blast_out.csv")
        output:
            os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "ani.tsv")
        params:
            script = "workflow/scripts/clust_anicalc.py"
        log: os.path.join(logdir, "checkv_megablast_anicalc_layer_{layer}_chunk_{chunk}.log")
        benchmark: os.path.join(benchmarks, "checkv_megablast_anicalc_layer_{layer}_chunk_{chunk}.log")
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

    rule checkv_megablast_aniclust:
        name: "clustering.smk CheckV-MEGABLAST cluster by ani"
        input:
            fa = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "input.fa"),
            ani = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "ani.tsv")
        output:
            clstr = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}.fa.clstr"),
            reps = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "cluster_representatives.txt")
        params:
            script = "workflow/scripts/clust_ani.py",
            minani = config.get("checkv-megablast-ani", ""),
            targetcov = config.get("checkv-megablast-targetcov", ""),
            querycov = config.get("checkv-megablast-querycov", ""),
            outdir = lambda wildcards, output: os.path.dirname(output.clstr)
        log: os.path.join(logdir, "checkv_megablast_aniclust_layer_{layer}_chunk_{chunk}.log")
        benchmark: os.path.join(benchmarks, "checkv_megablast_aniclust_layer_{layer}_chunk_{chunk}.log")
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

    rule checkv_megablast_filter_contigs:
        name: "clustering.smk CheckV-MEGABLAST filter dereplicated viral contigs"
        input:
            fna = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "input.fa"),
            reps = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "cluster_representatives.txt")
        output:
            fa = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}.fa")
        params:
            tmpdir = lambda wildcards, input: os.path.dirname(input.fna),
            outdir = lambda wildcards, output: os.path.dirname(output.fa)
        log: os.path.join(logdir, "checkv_megablast_filtercontigs_layer_{layer}_chunk_{chunk}.log")
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

    rule checkv_megablast_derep_finalize:
        name: "clustering.smk CheckV-MEGABLAST finalize dereplicated vOTUs"
        localrule: True
        input:
            fa = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), f"layer_{cluster_iter}", "chunk_0.fa"),
            clstr = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), f"layer_{cluster_iter}", "chunk_0.fa.clstr")
        output:
            fa = relpath(f"identify/viral/output/derep/{cluster_method}/combined.viralcontigs.derep.fa") if cluster_method_input == "all" else relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa"),
            clstr = relpath(f"identify/viral/output/derep/{cluster_method}/combined.viralcontigs.derep.fa.clstr") if cluster_method_input == "all" else relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa.clstr"),
        params: 
            outdir = lambda wildcards, output: os.path.dirname(output.fa)
        log: os.path.join(logdir, "check_megablast_finalize.log")
        benchmark: 
            os.path.join(benchmarks, "check_megablast_finalize.log")
        shell:
            """
            mkdir -p {params.outdir}

            cp {input.fa} {output.fa}
            cp {input.clstr} {output.clstr}
            """


# ----- CD-HIT-EST -----
if cluster_method_input in ["cd-hit-est", "all"]:
    cluster_method = "cd-hit-est"
    _method = cluster_method

    def get_iter_inputs(wildcards, method=_method):
        layer = int(wildcards.layer)
        chunk = int(wildcards.chunk)
        if layer == 1:
            return os.path.join(relpath("identify/viral/tmp/derep/cluster-splits"), f"chunk_{chunk}.fa")
        else:
            prev_layer = layer - 1
            child1 = chunk * 2
            child2 = chunk * 2 + 1
            return [
                os.path.join(relpath(f"identify/viral/tmp/derep/{method}/cluster-layers"), f"layer_{prev_layer}", f"chunk_{child1}.fa"),
                os.path.join(relpath(f"identify/viral/tmp/derep/{method}/cluster-layers"), f"layer_{prev_layer}", f"chunk_{child2}.fa")
            ]

    rule cdhit_recursive_cluster:
        name: "clustering.smk CD-HIT-EST recursive clustering"
        input:
            get_iter_inputs
        output:
            fa = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), f"layer_{{layer}}", f"chunk_{{chunk}}.fa"),
            clstr = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), f"layer_{{layer}}", f"chunk_{{chunk}}.fa.clstr")
        params:
            cdhitparams = config.get("cdhit-params", ""),
            outdir = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), f"layer_{{layer}}"), 
            tmpdir = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), f"layer_{{layer}}", f"chunk_{{chunk}}")
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

    rule cdhit_derep_finalize:
        name: "clustering.smk CD-HIT-EST finalize dereplicated vOTUs"
        localrule: True
        input:
            fa = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), f"layer_{cluster_iter}", "chunk_0.fa"),
            clstr = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), f"layer_{cluster_iter}", "chunk_0.fa.clstr")
        output:
            fa = relpath(f"identify/viral/output/derep/{cluster_method}/combined.viralcontigs.derep.fa") if cluster_method_input == "all" else relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa"),
            clstr = relpath(f"identify/viral/output/derep/{cluster_method}/combined.viralcontigs.derep.fa.clstr") if cluster_method_input == "all" else relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa.clstr"),
        params: 
            outdir = lambda wildcards, output: os.path.dirname(output.fa)
        log: os.path.join(logdir, "cdhit_finalize.log")
        benchmark: 
            os.path.join(benchmarks, "cdhit_finalize.log")
        shell:
            """
            mkdir -p {params.outdir}

            cp {input.fa} {output.fa}
            cp {input.clstr} {output.clstr}
            """


# ----- Vclust -----
if cluster_method_input in ["vclust", "all"]:
    cluster_method = "vclust"
    _method = cluster_method

    def get_iter_inputs(wildcards, method=_method):
        layer = int(wildcards.layer)
        chunk = int(wildcards.chunk)
        if layer == 1:
            return os.path.join(relpath("identify/viral/tmp/derep/cluster-splits"), f"chunk_{chunk}.fa")
        else:
            prev_layer = layer - 1
            child1 = chunk * 2
            child2 = chunk * 2 + 1
            return [
                os.path.join(relpath(f"identify/viral/tmp/derep/{method}/cluster-layers"), f"layer_{prev_layer}", f"chunk_{child1}.fa"),
                os.path.join(relpath(f"identify/viral/tmp/derep/{method}/cluster-layers"), f"layer_{prev_layer}", f"chunk_{child2}.fa")
            ]

    rule vclust_prep_input:
        name: "clustering.smk Vclust prepare splits"
        localrule: True
        input:
            get_iter_inputs
        output:
            os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "input.fa")
        params:
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            tmpdir = lambda wildcards, output: os.path.dirname(output[0])
        log: os.path.join(logdir, "vclust_prep_splits_{layer}_chunk_{chunk}.log")
        threads: 1
        shell:
            """
            mkdir -p {params.tmpdir} {params.outdir}
            
            if [ "{wildcards.layer}" -eq "1" ]; then
              cp {input} {output} 2> {log}
            else
              cat {input} > {output} 2> {log}
            fi
            """

    rule vclust_filter:
        name: "clustering.smk Vclust filter"
        input:
            os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "input.fa")
        output:
            os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{{layer}}", "chunk_{{chunk}}", "fltr.txt")
        params:
            parameters = config.get("vclust-filter-params"),
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            tmpdir = lambda wildcards, output: os.path.dirname(output[0])
        log: os.path.join(logdir, "vclust_filter_layer_{layer}_chunk_{chunk}.log")
        benchmark: os.path.join(benchmarks, "vclust_filter_layer_{layer}_chunk_{chunk}.log")
        conda: "../envs/vclust.yml"
        threads: 1
        resources:
            mem_mb = lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
        shell:
            """
            rm -rf {params.tmpdir}
            mkdir -p {params.tmpdir} {params.outdir}

            vclust filter \
                -i {input} \
                -o {params.tmpdir}/tmp.txt \
                {params.parameters} &> {log}

            mv {params.tmpdir}/tmp.txt {output}
            """

    rule vclust_align:
        name: "clustering.smk Vclust align"
        input:
            fna = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "input.fa"),
            txt = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{{layer}}", "chunk_{{chunk}}", "fltr.txt")
        output:
            ani = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{{layer}}", "chunk_{{chunk}}", "ani.tsv"), 
            ids = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{{layer}}", "chunk_{{chunk}}", "ani.ids.tsv"), 
        params:
            parameters = config.get("vclust-align-params"),
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            tmpdir = lambda wildcards, output: os.path.dirname(output[0])
        log: os.path.join(logdir, "vclust_align_layer_{layer}_chunk_{chunk}.log")
        benchmark: os.path.join(benchmarks, "vclust_align_layer_{layer}_chunk_{chunk}.log")
        conda: "../envs/vclust.yml"
        threads: 64
        resources:
            mem_mb = lambda wildcards, threads, attempt: int(attempt * (threads / 64) * 72 * 10**3)
        shell:
            """
            rm -rf {params.tmpdir}
            mkdir -p {params.tmpdir} {params.outdir}

            vclust align \
                -i {input.fna} \
                --filter {input.txt} \
                -o {params.tmpdir}/ani.tsv \
                {params.parameters} &> {log}

            mv {params.tmpdir}/ani.tsv {output.ani}
            mv {params.tmpdir}/ani.ids.tsv {output.ids}
            """

    rule vclust_cluster:
        name: "clustering.smk Vclust cluster"
        input:
            ani = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{{layer}}", "chunk_{{chunk}}", "ani.tsv"), 
            ids = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{{layer}}", "chunk_{{chunk}}", "ani.ids.tsv"), 
        output:
            clstr = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}.fa.clstr"),
            reps = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "cluster_representatives.txt")
        params:
            parameters = config.get("vclust-cluster-params"),
            algo = config.get("vclust-cluster-algorithm"),
            outdir = lambda wildcards, output: os.path.dirname(output[0]),
            tmpdir = lambda wildcards, output: os.path.dirname(output[0])
        log: os.path.join(logdir, "vclust_cluster_layer_{layer}_chunk_{chunk}.log")
        benchmark: os.path.join(benchmarks, "vclust_cluster_layer_{layer}_chunk_{chunk}.log")
        conda: "../envs/vclust.yml"
        threads: 64
        resources:
            mem_mb = lambda wildcards, threads, attempt: int(attempt * (threads / 64) * 72 * 10**3)
        shell:
            """
            rm -rf {params.tmpdir}
            mkdir -p {params.tmpdir} {params.outdir}

            vclust cluster \
                --algorithm {params.algo} \
                -i {input.ani} \
                --ids {input.ids} \
                -o {params.tmpdir}/tmp.tsv \
                --out-repr \
                {params.parameters} &> {log}
            
            cut -f1 {params.tmpdir}/tmp.tsv > {output.reps} 2>> {log}
            mv {params.tmpdir}/tmp.tsv {output.clstr}
            """

    rule vclust_filter_contigs:
        name: "clustering.smk Vclust filter dereplicated viral contigs"
        input:
            fna = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "input.fa"),
            reps = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "cluster_representatives.txt")
        output:
            fa = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}.fa")
        params:
            tmpdir = lambda wildcards, input: os.path.dirname(input.fna),
            outdir = lambda wildcards, output: os.path.dirname(output.fa)
        log: os.path.join(logdir, "vclust_filter_contigs_layer_{layer}_chunk_{chunk}.log")
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

    rule vclust_derep_finalize:
        name: "clustering.smk Vclust finalize dereplicated vOTUs"
        localrule: True
        input:
            fa = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), f"layer_{cluster_iter}", "chunk_0.fa"),
            clstr = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), f"layer_{cluster_iter}", "chunk_0.fa.clstr")
        output:
            fa = relpath(f"identify/viral/output/derep/{cluster_method}/combined.viralcontigs.derep.fa") if cluster_method_input == "all" else relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa"),
            clstr = relpath(f"identify/viral/output/derep/{cluster_method}/combined.viralcontigs.derep.fa.clstr") if cluster_method_input == "all" else relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa.clstr"),
        params: 
            outdir = lambda wildcards, output: os.path.dirname(output.fa)
        log: os.path.join(logdir, "cvclust_finalize.log")
        benchmark: 
            os.path.join(benchmarks, "vclust_finalize.log")
        shell:
            """
            mkdir -p {params.outdir}

            cp {input.fa} {output.fa}
            cp {input.clstr} {output.clstr}
            """


# ----- Linclust -----
if cluster_method_input in ["linclust", "all"]:
    if cluster_method_input == "all":
        cluster_method = "linclust"
    else:
        cluster_method = "linclust"
    _method = cluster_method

    def get_iter_inputs(wildcards, method=_method):
        layer = int(wildcards.layer)
        chunk = int(wildcards.chunk)
        if layer == 1:
            return os.path.join(relpath("identify/viral/tmp/derep/cluster-splits"), f"chunk_{chunk}.fa")
        else:
            prev_layer = layer - 1
            child1 = chunk * 2
            child2 = chunk * 2 + 1
            return [
                os.path.join(relpath(f"identify/viral/tmp/derep/{method}/cluster-layers"), f"layer_{prev_layer}", f"chunk_{child1}.fa"),
                os.path.join(relpath(f"identify/viral/tmp/derep/{method}/cluster-layers"), f"layer_{prev_layer}", f"chunk_{child2}.fa")
            ]

    rule linclust_prep_input:
        name: "clustering.smk MMSeqs2-Linclust prepare split input"
        localrule: True
        input:
            get_iter_inputs
        output:
            os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "input.fa")
        params:
            tmpdir = lambda wildcards, output: os.path.dirname(output[0])
        log: os.path.join(logdir, "linclust_prep_layer_{layer}_chunk_{chunk}.log")
        threads: 1
        shell:
            """
            mkdir -p {params.tmpdir}
            if [ "{wildcards.layer}" -eq "1" ]; then
                cp {input} {output} 2> {log}
            else
                cat {input} > {output} 2> {log}
            fi
            """

    rule linclust_makedb:
        name: "clustering.smk MMSeqs2-Linclust make db"
        input:
            os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "input.fa")
        output:
            db = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "linclust.db")
        params:
            outdir = lambda wildcards, input: os.path.dirname(input[0]),
            parameters = config.get("linclust-db-params", "")
        log: os.path.join(logdir, "linclust_makedb_layer_{layer}_chunk_{chunk}.log")
        benchmark: os.path.join(benchmarks, "linclust_makedb_layer_{layer}_chunk_{chunk}.log")
        conda: "../envs/mmseqs2.yml"
        threads: 1
        resources:
            mem_mb = lambda wildcards, attempt, input: max(2 * input.size_mb, 1000)
        shell:
            """
            mkdir -p {params.outdir}
            mmseqs createdb {input} {output.db} {params.parameters} &> {log}
            """

    rule linclust_run:
        name: "clustering.smk MMSeqs2-Linclust run"
        input:
            db = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "linclust.db")
        output:
            # MMseqs2 creates multiple files, so we use a directory
            clstr_db = directory(os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "cluster_db"))
        params:
            outdir = lambda wildcards, input: os.path.dirname(input.db),
            parameters = config.get("linclust-params", "")
        log: os.path.join(logdir, "linclust_run_layer_{layer}_chunk_{chunk}.log")
        benchmark: os.path.join(benchmarks, "linclust_run_layer_{layer}_chunk_{chunk}.log")
        conda: "../envs/mmseqs2.yml"
        threads: 64
        resources:
            mem_mb = lambda wildcards, threads, attempt: int(attempt * (threads / 64) * 72 * 10**3)
        shell:
            """
            mkdir -p {params.outdir}
            # The output prefix for mmseqs linclust is the directory path without trailing slash
            mmseqs linclust {input.db} {output.clstr_db} {params.outdir} {params.parameters} &> {log}
            """

    rule linclust_filter_contigs:
        name: "clustering.smk MMSeqs2-Linclust filter dereplicated viral contigs"
        input:
            db = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "linclust.db"),
            clstr_db = directory(os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "cluster_db"))
        output:
            reps = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}", "cluster_representatives.txt"),
            fa = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}.fa"),
            clstr = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}.fa.clstr")
        params:
            tmpdir = lambda wildcards, input: os.path.dirname(input.db),
            outdir = lambda wildcards, output: os.path.dirname(output.fa)
        log: os.path.join(logdir, "linclust_filter_contigs_layer_{layer}_chunk_{chunk}.log")
        conda: "../envs/mmseqs2.yml"
        threads: 1
        resources:
            mem_mb = lambda wildcards, attempt, input: max(2 * input.size_mb, 1000)
        shell:
            """
            mkdir -p {params.outdir} {params.tmpdir}

            mmseqs result2repseq {input.db} {input.clstr_db} {params.tmpdir}/rep_db &> {log}
            mmseqs convert2fasta {params.tmpdir}/rep_db {output.fa} &>> {log}

            mmseqs result2tsv {input.db} {input.clstr_db} {params.tmpdir}/cluster.tsv &>> {log}
            cp {params.tmpdir}/cluster.tsv {output.clstr}

            cut -f1 {params.tmpdir}/cluster.tsv | sort -u > {output.reps}

            rm -rf {params.tmpdir}/rep_db*
            rm -f {params.tmpdir}/cluster.tsv
            """

    rule linclust_derep_finalize:
        name: "clustering.smk MMSeqs2-Linclust finalize dereplicated vOTUs"
        localrule: True
        input:
            fa = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), f"layer_{cluster_iter}", "chunk_0.fa"),
            clstr = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), f"layer_{cluster_iter}", "chunk_0.fa.clstr")
        output:
            fa = relpath(f"identify/viral/output/derep/{cluster_method}/combined.viralcontigs.derep.fa") if cluster_method_input == "all" else relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa"),
            clstr = relpath(f"identify/viral/output/derep/{cluster_method}/combined.viralcontigs.derep.fa.clstr") if cluster_method_input == "all" else relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa.clstr"),
        params:
            outdir = lambda wildcards, output: os.path.dirname(output.fa)
        log: os.path.join(logdir, "linclust_finalize.log")
        benchmark: os.path.join(benchmarks, "linclust_finalize.log")
        shell:
            """
            mkdir -p {params.outdir}
            cp {input.fa} {output.fa}
            cp {input.clstr} {output.clstr}
            """


# ----- DNACLUST -----
if cluster_method_input in ["dnaclust", "all"]:
    cluster_method = "dnaclust"
    _method = cluster_method

    def get_iter_inputs(wildcards, method=_method):
        layer = int(wildcards.layer)
        chunk = int(wildcards.chunk)
        if layer == 1:
            return os.path.join(relpath("identify/viral/tmp/derep/cluster-splits"), f"chunk_{chunk}.fa")
        else:
            prev_layer = layer - 1
            child1 = chunk * 2
            child2 = chunk * 2 + 1
            return [
                os.path.join(relpath(f"identify/viral/tmp/derep/{method}/cluster-layers"), f"layer_{prev_layer}", f"chunk_{child1}.fa"),
                os.path.join(relpath(f"identify/viral/tmp/derep/{method}/cluster-layers"), f"layer_{prev_layer}", f"chunk_{child2}.fa")
            ]

    rule dnaclust_recursive_cluster:
        name: "clustering.smk DNACLUST recursive clustering"
        input:
            get_iter_inputs
        output:
            fa = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}.fa"),
            clstr = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}.fa.clstr")
        params:
            similarity = lambda wildcards: config.get("dnaclust-similarity", 0.95),
            extra = config.get("dnaclust-params", ""),
            outdir = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}"),
            tmpdir = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}")
        threads: 8
        resources:
            mem_mb = lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
        conda: "../envs/dnaclust.yml"
        log: os.path.join(logdir, "dnaclust_layer_{layer}_chunk_{chunk}.log")
        benchmark: os.path.join(benchmarks, "dnaclust_layer_{layer}_chunk_{chunk}.log")
        shell:
            """
            rm -rf {params.tmpdir}
            mkdir -p {params.tmpdir} {params.outdir}

            if [ "{wildcards.layer}" -eq "1" ]; then
                dnaclust {input} \
                    -s {params.similarity} \
                    {params.extra} \
                    --representatives {params.tmpdir}/reps.fa \
                    -o {params.tmpdir}/clust &>> {log}
            else
                cat {input} > {params.tmpdir}/pooled.fa
                dnaclust {params.tmpdir}/pooled.fa \
                    -s {params.similarity} \
                    {params.extra} \
                    --representatives {params.tmpdir}/reps.fa \
                    -o {params.tmpdir}/clust &>> {log}
                rm -f {params.tmpdir}/pooled.fa
            fi

            mv {params.tmpdir}/reps.fa {output.fa}
            mv {params.tmpdir}/clust_clusters.txt {output.clstr}
            rm -rf {params.tmpdir}
            """

    rule dnaclust_derep_finalize:
        name: "clustering.smk DNACLUST finalize dereplicated vOTUs"
        localrule: True
        input:
            fa = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), f"layer_{cluster_iter}", "chunk_0.fa"),
            clstr = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), f"layer_{cluster_iter}", "chunk_0.fa.clstr")
        output:
            fa = relpath(f"identify/viral/output/derep/{cluster_method}/combined.viralcontigs.derep.fa") if cluster_method_input == "all" else relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa"),
            clstr = relpath(f"identify/viral/output/derep/{cluster_method}/combined.viralcontigs.derep.fa.clstr") if cluster_method_input == "all" else relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa.clstr"),
        params:
            outdir = lambda wildcards, output: os.path.dirname(output.fa)
        log: os.path.join(logdir, "dnaclust_finalize.log")
        benchmark: os.path.join(benchmarks, "dnaclust_finalize.log")
        shell:
            """
            mkdir -p {params.outdir}
            cp {input.fa} {output.fa}
            cp {input.clstr} {output.clstr}
            """

# ----- VSEARCH -----
if cluster_method_input in ["vsearch", "all"]:
    cluster_method = "vsearch"
    _method = cluster_method

    def get_iter_inputs(wildcards, method=_method):
        layer = int(wildcards.layer)
        chunk = int(wildcards.chunk)
        if layer == 1:
            return os.path.join(relpath("identify/viral/tmp/derep/cluster-splits"), f"chunk_{chunk}.fa")
        else:
            prev_layer = layer - 1
            child1 = chunk * 2
            child2 = chunk * 2 + 1
            return [
                os.path.join(relpath(f"identify/viral/tmp/derep/{method}/cluster-layers"), f"layer_{prev_layer}", f"chunk_{child1}.fa"),
                os.path.join(relpath(f"identify/viral/tmp/derep/{method}/cluster-layers"), f"layer_{prev_layer}", f"chunk_{child2}.fa")
            ]

    rule vsearch_recursive_cluster:
        name: "clustering.smk VSEARCH recursive clustering"
        input:
            get_iter_inputs
        output:
            fa = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}.fa"),
            clstr = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}.fa.clstr")
        params:
            identity = lambda wildcards: config.get("vsearch-identity", 0.95),
            extra = config.get("vsearch-params", ""),
            outdir = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}"),
            tmpdir = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}")
        threads: 8
        resources:
            mem_mb = lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
        conda: "../envs/vsearch.yml"
        log: os.path.join(logdir, "vsearch_layer_{layer}_chunk_{chunk}.log")
        benchmark: os.path.join(benchmarks, "vsearch_layer_{layer}_chunk_{chunk}.log")
        shell:
            """
            rm -rf {params.tmpdir}
            mkdir -p {params.tmpdir} {params.outdir}

            if [ "{wildcards.layer}" -eq "1" ]; then
                vsearch --cluster_fast {input} \
                    --id {params.identity} \
                    {params.extra} \
                    --centroids {params.tmpdir}/centroids.fa \
                    --uc {params.tmpdir}/clusters.uc \
                    --threads {threads} &>> {log}
            else
                cat {input} > {params.tmpdir}/pooled.fa
                vsearch --cluster_fast {params.tmpdir}/pooled.fa \
                    --id {params.identity} \
                    {params.extra} \
                    --centroids {params.tmpdir}/centroids.fa \
                    --uc {params.tmpdir}/clusters.uc \
                    --threads {threads} &>> {log}
                rm -f {params.tmpdir}/pooled.fa
            fi

            mv {params.tmpdir}/centroids.fa {output.fa}
            mv {params.tmpdir}/clusters.uc {output.clstr}
            rm -rf {params.tmpdir}
            """

    rule vsearch_derep_finalize:
        name: "clustering.smk VSEARCH finalize dereplicated vOTUs"
        localrule: True
        input:
            fa = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), f"layer_{cluster_iter}", "chunk_0.fa"),
            clstr = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), f"layer_{cluster_iter}", "chunk_0.fa.clstr")
        output:
            fa = relpath(f"identify/viral/output/derep/{cluster_method}/combined.viralcontigs.derep.fa") if cluster_method_input == "all" else relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa"),
            clstr = relpath(f"identify/viral/output/derep/{cluster_method}/combined.viralcontigs.derep.fa.clstr") if cluster_method_input == "all" else relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa.clstr"),
        params:
            outdir = lambda wildcards, output: os.path.dirname(output.fa)
        log: os.path.join(logdir, "vsearch_finalize.log")
        benchmark: os.path.join(benchmarks, "vsearch_finalize.log")
        shell:
            """
            mkdir -p {params.outdir}
            cp {input.fa} {output.fa}
            cp {input.clstr} {output.clstr}
            """

# ----- VIRIDIC -----
if cluster_method_input in ["viridic", "all"]:
    cluster_method = "viridic"
    _method = cluster_method

    def get_iter_inputs(wildcards, method=_method):
        layer = int(wildcards.layer)
        chunk = int(wildcards.chunk)
        if layer == 1:
            return os.path.join(relpath("identify/viral/tmp/derep/cluster-splits"), f"chunk_{chunk}.fa")
        else:
            prev_layer = layer - 1
            child1 = chunk * 2
            child2 = chunk * 2 + 1
            return [
                os.path.join(relpath(f"identify/viral/tmp/derep/{method}/cluster-layers"), f"layer_{prev_layer}", f"chunk_{child1}.fa"),
                os.path.join(relpath(f"identify/viral/tmp/derep/{method}/cluster-layers"), f"layer_{prev_layer}", f"chunk_{child2}.fa")
            ]

    rule viridic_recursive_cluster:
        name: "clustering.smk VIRIDIC recursive clustering"
        input:
            get_iter_inputs
        output:
            fa = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}.fa"),
            clstr = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}.fa.clstr")
        params:
            threshold = config.get("viridic-threshold", 0.95),
            extra = config.get("viridic-params", ""),
            outdir = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}"),
            tmpdir = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), "layer_{layer}", "chunk_{chunk}"), 
            prefix = lambda wildcards, input: (os.path.splitext(os.path.basename(input[0]))[0] if wildcards.layer == "1" else "pooled")
        threads: 1
        resources:
            mem_mb = lambda wildcards, attempt, input: max(2*input.size_mb, 1000)
        conda: "../envs/viridic.yml"
        log: os.path.join(logdir, "viridic_layer_{layer}_chunk_{chunk}.log")
        benchmark: os.path.join(benchmarks, "viridic_layer_{layer}_chunk_{chunk}.log")
        shell:
            """
            rm -rf {params.tmpdir}
            mkdir -p {params.tmpdir} {params.outdir}

            if [ "{wildcards.layer}" -eq "1" ]; then
                viridic.py -i {input} \
                    -o {params.tmpdir} \
                    -t {params.threshold} \
                    {params.extra} &>> {log}
            else
                cat {input} > {params.tmpdir}/pooled.fa
                viridic.py -i {params.tmpdir}/pooled.fa \
                    -o {params.tmpdir} \
                    -t {params.threshold} \
                    {params.extra} &>> {log}
                rm -f {params.tmpdir}/pooled.fa
            fi

            # Now we know the exact output names
            prefix="{params.prefix}"
            mv {params.tmpdir}/${{prefix}}_representatives.fna {output.fa}
            mv {params.tmpdir}/${{prefix}}_clusters.tsv {output.clstr}

            rm -rf {params.tmpdir}
            """

    rule viridic_derep_finalize:
        name: "clustering.smk VIRIDIC finalize dereplicated vOTUs"
        localrule: True
        input:
            fa = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), f"layer_{cluster_iter}", "chunk_0.fa"),
            clstr = os.path.join(relpath(f"identify/viral/tmp/derep/{cluster_method}/cluster-layers"), f"layer_{cluster_iter}", "chunk_0.fa.clstr")
        output:
            fa = relpath(f"identify/viral/output/derep/{cluster_method}/combined.viralcontigs.derep.fa") if cluster_method_input == "all" else relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa"),
            clstr = relpath(f"identify/viral/output/derep/{cluster_method}/combined.viralcontigs.derep.fa.clstr") if cluster_method_input == "all" else relpath("identify/viral/output/derep/combined.viralcontigs.derep.fa.clstr"),
        params:
            outdir = lambda wildcards, output: os.path.dirname(output.fa)
        log: os.path.join(logdir, "viridic_finalize.log")
        benchmark: os.path.join(benchmarks, "viridic_finalize.log")
        shell:
            """
            mkdir -p {params.outdir}
            cp {input.fa} {output.fa}
            cp {input.clstr} {output.clstr}
            """
