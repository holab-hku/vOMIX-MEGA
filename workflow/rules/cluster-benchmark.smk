# ----------------------------------------------------------------------
# Configuration & setup
# ----------------------------------------------------------------------
vomix_module = vomix_utils.current_module
logdir = vomix_module.logdir
benchmarks = vomix_module.benchmarks
tmpd = vomix_module.tmpd

# ------------------------------------------------------------
# Setup mock dataset parameters
# ------------------------------------------------------------

DATASET_PARAMS = {
    "Mock-10K":          {"size": 10000,  "virus": 0.5,  "prok": 0.3,  "euk": 0.2,  "strain": 0, "species": 10, "seed": config.get("seed", 42)},
    "Mock-10K-HighVir":  {"size": 10000,  "virus": 1.0,  "prok": 0.0,  "euk": 0.0,  "strain": 0, "species": 10, "seed": config.get("seed", 42)},
    "Mock-10K-LowVir":   {"size": 10000,  "virus": 0.1,  "prok": 0.5,  "euk": 0.4,  "strain": 0, "species": 10, "seed": config.get("seed", 42)},
    "Mock-50K":          {"size": 50000,  "virus": 0.5,  "prok": 0.3,  "euk": 0.2,  "strain": 0, "species": 10, "seed": config.get("seed", 42)},
    "Mock-100K":         {"size": 100000, "virus": 0.5,  "prok": 0.3,  "euk": 0.2,  "strain": 0, "species": 10, "seed": config.get("seed", 42)},
    "Mock-300K":         {"size": 300000, "virus": 0.5,  "prok": 0.3,  "euk": 0.2,  "strain": 0, "species": 10, "seed": config.get("seed", 42)},
    "Mock-1000K":        {"size": 1000000,"virus": 0.5,  "prok": 0.3,  "euk": 0.2,  "strain": 0, "species": 10, "seed": config.get("seed", 42)},
    "Mock-Strain":       {"size": 20000,  "virus": 1.0,  "prok": 0.0,  "euk": 0.0,  "strain": 1, "species": 10, "seed": config.get("seed", 42)},
}

def get_resources(dataset):
    """Return memory, disk, and threads based on dataset size."""
    size = DATASET_PARAMS[dataset]["size"]
    if size <= 10000:
        return {"mem_mb": 8192,  "disk_mb": 16000, "threads": 8}
    elif size <= 50000:
        return {"mem_mb": 16384, "disk_mb": 32000, "threads": 8}
    elif size <= 100000:
        return {"mem_mb": 16384, "disk_mb": 32000, "threads": 8}
    elif size <= 300000:
        return {"mem_mb": 32768, "disk_mb": 64000, "threads": 12}
    else:  # 1M
        return {"mem_mb": 65536, "disk_mb": 128000, "threads": 16}


# ------------------------------------------------------------
# MASTER RULE
# ------------------------------------------------------------
rule cluster_benchmark_done:
    name: "cluster-benchmark.smk Done. removing tmp files"
    localrule: True
    input:
        os.path.join(datadir, "mock-data", "genomes", "refseq_viral.fna"),
        os.path.join(datadir, "mock-data", "genomes", "refseq_prok.fna"),
        os.path.join(datadir, "mock-data", "genomes", "refseq_euk.fna"),
        expand(os.path.join(datadir, "mock-data", "{dataset}.fna"), dataset=DATASET_PARAMS.keys()),
        expand(os.path.join(datadir, "mock-data", "{dataset}.ground_truth.tsv"), dataset=DATASET_PARAMS.keys()),
        # os.path.join(basedir, "cami_marine", "README.txt"),
    output:
        os.path.join(logdir, "done.log")
    shell:
        "touch {output}"


# ------------------------------------------------------------
# RULES
# ------------------------------------------------------------
rule download_refseq_viral:
    name: "cluster-benchmark.smk downloading viral genomes"
    output:
        fna = os.path.join(datadir, "mock-data", "genomes", "refseq_viral.fna"),
        tsv = os.path.join(datadir, "mock-data", "genomes", "refseq_viral_summary.tsv")
    params:
        email    = config.get("NCBI-email", ""),
        api_key  = config.get("NCBI-API-key", ""),
        script   = "../scripts/download_refseq_genomes.py",
        outdir   = datadir,
        tmpdir   = os.path.join(tmpd, "genomes", "viral"),
    conda: "../envs/seqkit-biopython.yml"
    resources:
        mem_mb=lambda wildcards, attempt, input: 4 * 10**3 * attempt
    shell:
        """
        rm -rf {params.tmpdir} {output.fna}
        mkdir -p {params.tmpdir} {params.outdir}

        python {params.script} \
            --outfile {params.tmpdir}/tmp.fna \
            --summary-tsv {params.tmpdir}/tmp.tsv \
            --email "{params.email}" \
            --api-key "{params.api_key}" \
            --mode viral \
            --genomes-per-category 0 \
            > {log} 2>&1

        mv {params.tmpdir}/tmp.fna {output.fna}
        mv {params.tmpdir}/tmp.tsv {output.tsv}
        """

rule download_prok_contaminants:
    name: "cluster-benchmark.smk downloading prokaryotic genomes"
    output:
        fna = os.path.join(datadir, "mock-data", "genomes", "refseq_prok.fna"),
        tsv = os.path.join(datadir, "mock-data", "genomes", "refseq_prok_summary.tsv")
    params:
        email    = config.get("NCBI-email", ""),
        api_key  = config.get("NCBI-API-key", ""),
        script   = "../scripts/download_refseq_genomes.py",
        outdir   = datadir,
        tmpdir   = os.path.join(tmpd, "genomes", "euk"),
    conda: "../envs/seqkit-biopython.yml"
    resources:
        mem_mb=lambda wildcards, attempt, input: 4 * 10**3 * attempt
    shell:
        """
        rm -rf {params.tmpdir}
        mkdir -p {params.tmpdir}

        python {params.script} \
            --outfile {params.tmpdir}/tmp.fna \
            --summary-tsv {params.tmpdir}/tmp.tsv \
            --email "{params.email}" \
            --api-key "{params.api_key}" \
            --mode prokaryotic \
            --verbose \
            > {log} 2>&1

        mv {params.tmpdir}/tmp.fna {output.fna}
        mv {params.tmpdir}/tmp.tsv {output.tsv}
        """

rule download_euk_contaminants:
    name: "cluster-benchmark.smk downloading eukaryotic genomes"
    output:
        fna = os.path.join(datadir, "mock-data", "genomes", "refseq_euk.fna"),
        tsv = os.path.join(datadir, "mock-data", "genomes", "refseq_euk_summary.tsv")
    params:
        email    = config.get("NCBI-email", ""),
        api_key  = config.get("NCBI-API-key", ""),
        script   = "../scripts/download_refseq_genomes.py",
        outdir   = datadir,
        tmpdir   = os.path.join(tmpd, "genomes", "euk"),
    conda: "../envs/seqkit-biopython.yml"
    resources:
        mem_mb=lambda wildcards, attempt, input: 4 * 10**3 * attempt
    shell:
        """
        rm -rf {params.tmpdir}
        mkdir -p {params.tmpdir}

        python {params.script} \
            --outfile {params.tmpdir}/tmp.fna \
            --summary-tsv {params.tmpdir}/tmp.tsv \
            --email "{params.email}" \
            --api-key "{params.api_key}" \
            --mode eukaryotic \
            --verbose \
            > {log} 2>&1

        mv {params.tmpdir}/tmp.fna {output.fna}
        mv {params.tmpdir}/tmp.tsv {output.tsv}
        """


rule mock_Mock_10K:
    name: "cluster-benchmark.smk Mock-10K dataset"
    output:
        fna = os.path.join(datadir, "mock-data", "Mock-10K.fna"),
        gt  = os.path.join(datadir, "mock-data", "Mock-10K.ground_truth.tsv")
    input:
        vir  = os.path.join(datadir, "mock-data", "genomes", "refseq_viral.fna"),
        prok = os.path.join(datadir, "mock-data", "genomes", "refseq_prok.fna"),
        euk  = os.path.join(datadir, "mock-data", "genomes", "refseq_euk.fna")
    params:
        tmpdir = os.path.join(tmpd, "Mock-10K"),
        name   = "Mock-10K",
        num    = DATASET_PARAMS["Mock-10K"]["size"],
        virus  = DATASET_PARAMS["Mock-10K"]["virus"],
        prok_f = DATASET_PARAMS["Mock-10K"]["prok"],
        euk_f  = DATASET_PARAMS["Mock-10K"]["euk"],
        strain = DATASET_PARAMS["Mock-10K"]["strain"],
        species = DATASET_PARAMS["Mock-10K"]["species"],
        seed   = DATASET_PARAMS["Mock-10K"]["seed"],
        script = "../scripts/generate_mock_clust_data.py",
    conda: "../envs/seqkit-biopython.yml"
    log: os.path.join(logdir, "mock_Mock-10K.log")
    benchmark: os.path.join(benchmarks, "mock_Mock-10K.benchmark")
    threads: get_resources("Mock-10K")["threads"]
    resources:
        mem_mb = get_resources("Mock-10K")["mem_mb"],
        disk_mb = get_resources("Mock-10K")["disk_mb"]
    shell:
        """
        set -euo pipefail
        rm -rf {params.tmpdir}
        mkdir -p {params.tmpdir}

        python {params.script} \
            --name {params.name} \
            --num-sequences {params.num} \
            --outdir {datadir}/mock-data \
            --viral-seq {input.vir} \
            --prokaryotic-seq {input.prok} \
            --eukaryotic-seq {input.euk} \
            --virus-frac {params.virus} \
            --prokaryote-frac {params.prok_f} \
            --eukaryote-frac {params.euk_f} \
            --strain-mode {params.strain} \
            --num-species {params.species} \
            --mut-rate-min 0.001 \
            --mut-rate-max 0.05 \
            --seed {params.seed} \
            --force \
            > {log} 2>&1

        touch {output.fna} {output.gt}
        """

rule mock_Mock_10K_HighVir:
    name: "cluster-benchmark.smk Mock-10K-HighVir dataset"
    output:
        fna = os.path.join(datadir, "mock-data", "Mock-10K-HighVir.fna"),
        gt  = os.path.join(datadir, "mock-data", "Mock-10K-HighVir.ground_truth.tsv")
    input:
        vir  = os.path.join(datadir, "mock-data", "genomes", "refseq_viral.fna"),
        prok = os.path.join(datadir, "mock-data", "genomes", "refseq_prok.fna"),
        euk  = os.path.join(datadir, "mock-data", "genomes", "refseq_euk.fna")
    params:
        tmpdir = os.path.join(tmpd, "Mock-10K-HighVir"),
        name   = "Mock-10K-HighVir",
        num    = DATASET_PARAMS["Mock-10K-HighVir"]["size"],
        virus  = DATASET_PARAMS["Mock-10K-HighVir"]["virus"],
        prok_f = DATASET_PARAMS["Mock-10K-HighVir"]["prok"],
        euk_f  = DATASET_PARAMS["Mock-10K-HighVir"]["euk"],
        strain = DATASET_PARAMS["Mock-10K-HighVir"]["strain"],
        species = DATASET_PARAMS["Mock-10K-HighVir"]["species"],
        seed   = DATASET_PARAMS["Mock-10K-HighVir"]["seed"],
        script = "../scripts/generate_mock_clust_data.py",
    conda: "../envs/seqkit-biopython.yml"
    log: os.path.join(logdir, "mock_Mock-10K-HighVir.log")
    benchmark: os.path.join(benchmarks, "mock_Mock-10K-HighVir.benchmark")
    threads: get_resources("Mock-10K-HighVir")["threads"]
    resources:
        mem_mb = get_resources("Mock-10K-HighVir")["mem_mb"],
        disk_mb = get_resources("Mock-10K-HighVir")["disk_mb"]
    shell:
        """
        set -euo pipefail
        rm -rf {params.tmpdir}
        mkdir -p {params.tmpdir}

        python {params.script} \
            --name {params.name} \
            --num-sequences {params.num} \
            --outdir {datadir}/mock-data \
            --tmpdir {params.tmpdir} \
            --viral-seq {input.vir} \
            --prokaryotic-seq {input.prok} \
            --eukaryotic-seq {input.euk} \
            --virus-frac {params.virus} \
            --prokaryote-frac {params.prok_f} \
            --eukaryote-frac {params.euk_f} \
            --strain-mode {params.strain} \
            --num-species {params.species} \
            --mut-rate-min 0.001 \
            --mut-rate-max 0.05 \
            --seed {params.seed} \
            --force \
            > {log} 2>&1

        touch {output.fna} {output.gt}
        """

rule mock_Mock_10K_LowVir:
    name: "cluster-benchmark.smk Mock-10K-LowVir dataset"
    output:
        fna = os.path.join(datadir, "mock-data", "Mock-10K-LowVir.fna"),
        gt  = os.path.join(datadir, "mock-data", "Mock-10K-LowVir.ground_truth.tsv")
    input:
        vir  = os.path.join(datadir, "mock-data", "genomes", "refseq_viral.fna"),
        prok = os.path.join(datadir, "mock-data", "genomes", "refseq_prok.fna"),
        euk  = os.path.join(datadir, "mock-data", "genomes", "refseq_euk.fna")
    params:
        tmpdir = os.path.join(tmpd, "Mock-10K-LowVir"),
        name   = "Mock-10K-LowVir",
        num    = DATASET_PARAMS["Mock-10K-LowVir"]["size"],
        virus  = DATASET_PARAMS["Mock-10K-LowVir"]["virus"],
        prok_f = DATASET_PARAMS["Mock-10K-LowVir"]["prok"],
        euk_f  = DATASET_PARAMS["Mock-10K-LowVir"]["euk"],
        strain = DATASET_PARAMS["Mock-10K-LowVir"]["strain"],
        species = DATASET_PARAMS["Mock-10K-LowVir"]["species"],
        seed   = DATASET_PARAMS["Mock-10K-LowVir"]["seed"],
        script = "../scripts/generate_mock_clust_data.py",
    conda: "../envs/seqkit-biopython.yml"
    log: os.path.join(logdir, "mock_Mock-10K-LowVir.log")
    benchmark: os.path.join(benchmarks, "mock_Mock-10K-LowVir.benchmark")
    threads: get_resources("Mock-10K-LowVir")["threads"]
    resources:
        mem_mb = get_resources("Mock-10K-LowVir")["mem_mb"],
        disk_mb = get_resources("Mock-10K-LowVir")["disk_mb"]
    shell:
        """
        set -euo pipefail
        rm -rf {params.tmpdir}
        mkdir -p {params.tmpdir}

        python {params.script} \
            --name {params.name} \
            --num-sequences {params.num} \
            --outdir {datadir}/mock-data \
            --tmpdir {params.tmpdir} \
            --viral-seq {input.vir} \
            --prokaryotic-seq {input.prok} \
            --eukaryotic-seq {input.euk} \
            --virus-frac {params.virus} \
            --prokaryote-frac {params.prok_f} \
            --eukaryote-frac {params.euk_f} \
            --strain-mode {params.strain} \
            --num-species {params.species} \
            --mut-rate-min 0.001 \
            --mut-rate-max 0.05 \
            --seed {params.seed} \
            --force \
            > {log} 2>&1

        touch {output.fna} {output.gt}
        """

rule mock_Mock_50K:
    name: "cluster-benchmark.smk Mock-50K dataset"
    output:
        fna = os.path.join(datadir, "mock-data", "Mock-50K.fna"),
        gt  = os.path.join(datadir, "mock-data", "Mock-50K.ground_truth.tsv")
    input:
        vir  = os.path.join(datadir, "mock-data", "genomes", "refseq_viral.fna"),
        prok = os.path.join(datadir, "mock-data", "genomes", "refseq_prok.fna"),
        euk  = os.path.join(datadir, "mock-data", "genomes", "refseq_euk.fna")
    params:
        tmpdir = os.path.join(tmpd, "Mock-50K"),
        name   = "Mock-50K",
        num    = DATASET_PARAMS["Mock-50K"]["size"],
        virus  = DATASET_PARAMS["Mock-50K"]["virus"],
        prok_f = DATASET_PARAMS["Mock-50K"]["prok"],
        euk_f  = DATASET_PARAMS["Mock-50K"]["euk"],
        strain = DATASET_PARAMS["Mock-50K"]["strain"],
        species = DATASET_PARAMS["Mock-50K"]["species"],
        seed   = DATASET_PARAMS["Mock-50K"]["seed"],
        script = "../scripts/generate_mock_clust_data.py",
    conda: "../envs/seqkit-biopython.yml"
    log: os.path.join(logdir, "mock_Mock-50K.log")
    benchmark: os.path.join(benchmarks, "mock_Mock-50K.benchmark")
    threads: get_resources("Mock-50K")["threads"]
    resources:
        mem_mb = get_resources("Mock-50K")["mem_mb"],
        disk_mb = get_resources("Mock-50K")["disk_mb"]
    shell:
        """
        set -euo pipefail
        rm -rf {params.tmpdir}
        mkdir -p {params.tmpdir}

        python {params.script} \
            --name {params.name} \
            --num-sequences {params.num} \
            --outdir {datadir}/mock-data \
            --tmpdir {params.tmpdir} \
            --viral-seq {input.vir} \
            --prokaryotic-seq {input.prok} \
            --eukaryotic-seq {input.euk} \
            --virus-frac {params.virus} \
            --prokaryote-frac {params.prok_f} \
            --eukaryote-frac {params.euk_f} \
            --strain-mode {params.strain} \
            --num-species {params.species} \
            --mut-rate-min 0.001 \
            --mut-rate-max 0.05 \
            --seed {params.seed} \
            --force \
            > {log} 2>&1

        touch {output.fna} {output.gt}
        """

rule mock_Mock_100K:
    name: "cluster-benchmark.smk Mock-100K dataset"
    output:
        fna = os.path.join(datadir, "mock-data", "Mock-100K.fna"),
        gt  = os.path.join(datadir, "mock-data", "Mock-100K.ground_truth.tsv")
    input:
        vir  = os.path.join(datadir, "mock-data", "genomes", "refseq_viral.fna"),
        prok = os.path.join(datadir, "mock-data", "genomes", "refseq_prok.fna"),
        euk  = os.path.join(datadir, "mock-data", "genomes", "refseq_euk.fna")
    params:
        tmpdir = os.path.join(tmpd, "Mock-100K"),
        name   = "Mock-100K",
        num    = DATASET_PARAMS["Mock-100K"]["size"],
        virus  = DATASET_PARAMS["Mock-100K"]["virus"],
        prok_f = DATASET_PARAMS["Mock-100K"]["prok"],
        euk_f  = DATASET_PARAMS["Mock-100K"]["euk"],
        strain = DATASET_PARAMS["Mock-100K"]["strain"],
        species = DATASET_PARAMS["Mock-100K"]["species"],
        seed   = DATASET_PARAMS["Mock-100K"]["seed"],
        script = "../scripts/generate_mock_clust_data.py",
    conda: "../envs/seqkit-biopython.yml"
    log: os.path.join(logdir, "mock_Mock-100K.log")
    benchmark: os.path.join(benchmarks, "mock_Mock-100K.benchmark")
    threads: get_resources("Mock-100K")["threads"]
    resources:
        mem_mb = get_resources("Mock-100K")["mem_mb"],
        disk_mb = get_resources("Mock-100K")["disk_mb"]
    shell:
        """
        set -euo pipefail
        rm -rf {params.tmpdir}
        mkdir -p {params.tmpdir}

        python {params.script} \
            --name {params.name} \
            --num-sequences {params.num} \
            --outdir {datadir}/mock-data \
            --tmpdir {params.tmpdir} \
            --viral-seq {input.vir} \
            --prokaryotic-seq {input.prok} \
            --eukaryotic-seq {input.euk} \
            --virus-frac {params.virus} \
            --prokaryote-frac {params.prok_f} \
            --eukaryote-frac {params.euk_f} \
            --strain-mode {params.strain} \
            --num-species {params.species} \
            --mut-rate-min 0.001 \
            --mut-rate-max 0.05 \
            --seed {params.seed} \
            --force \
            > {log} 2>&1

        touch {output.fna} {output.gt}
        """

rule mock_Mock_300K:
    name: "cluster-benchmark.smk Mock-300K dataset"
    output:
        fna = os.path.join(datadir, "mock-data", "Mock-300K.fna"),
        gt  = os.path.join(datadir, "mock-data", "Mock-300K.ground_truth.tsv")
    input:
        vir  = os.path.join(datadir, "mock-data", "genomes", "refseq_viral.fna"),
        prok = os.path.join(datadir, "mock-data", "genomes", "refseq_prok.fna"),
        euk  = os.path.join(datadir, "mock-data", "genomes", "refseq_euk.fna")
    params:
        tmpdir = os.path.join(tmpd, "Mock-300K"),
        name   = "Mock-300K",
        num    = DATASET_PARAMS["Mock-300K"]["size"],
        virus  = DATASET_PARAMS["Mock-300K"]["virus"],
        prok_f = DATASET_PARAMS["Mock-300K"]["prok"],
        euk_f  = DATASET_PARAMS["Mock-300K"]["euk"],
        strain = DATASET_PARAMS["Mock-300K"]["strain"],
        species = DATASET_PARAMS["Mock-300K"]["species"],
        seed   = DATASET_PARAMS["Mock-300K"]["seed"],
        script = "../scripts/generate_mock_clust_data.py",
    conda: "../envs/seqkit-biopython.yml"
    log: os.path.join(logdir, "mock_Mock-300K.log")
    benchmark: os.path.join(benchmarks, "mock_Mock-300K.benchmark")
    threads: get_resources("Mock-300K")["threads"]
    resources:
        mem_mb = get_resources("Mock-300K")["mem_mb"],
        disk_mb = get_resources("Mock-300K")["disk_mb"]
    shell:
        """
        set -euo pipefail
        rm -rf {params.tmpdir}
        mkdir -p {params.tmpdir}

        python {params.script} \
            --name {params.name} \
            --num-sequences {params.num} \
            --outdir {datadir}/mock-data \
            --tmpdir {params.tmpdir} \
            --viral-seq {input.vir} \
            --prokaryotic-seq {input.prok} \
            --eukaryotic-seq {input.euk} \
            --virus-frac {params.virus} \
            --prokaryote-frac {params.prok_f} \
            --eukaryote-frac {params.euk_f} \
            --strain-mode {params.strain} \
            --num-species {params.species} \
            --mut-rate-min 0.001 \
            --mut-rate-max 0.05 \
            --seed {params.seed} \
            --force \
            > {log} 2>&1

        touch {output.fna} {output.gt}
        """

rule mock_Mock_1000K:
    name: "cluster-benchmark.smk Mock-1000K dataset"
    output:
        fna = os.path.join(datadir, "mock-data", "Mock-1000K.fna"),
        gt  = os.path.join(datadir, "mock-data", "Mock-1000K.ground_truth.tsv")
    input:
        vir  = os.path.join(datadir, "mock-data", "genomes", "refseq_viral.fna"),
        prok = os.path.join(datadir, "mock-data", "genomes", "refseq_prok.fna"),
        euk  = os.path.join(datadir, "mock-data", "genomes", "refseq_euk.fna")
    params:
        tmpdir = os.path.join(tmpd, "Mock-1000K"),
        name   = "Mock-1000K",
        num    = DATASET_PARAMS["Mock-1000K"]["size"],
        virus  = DATASET_PARAMS["Mock-1000K"]["virus"],
        prok_f = DATASET_PARAMS["Mock-1000K"]["prok"],
        euk_f  = DATASET_PARAMS["Mock-1000K"]["euk"],
        strain = DATASET_PARAMS["Mock-1000K"]["strain"],
        species = DATASET_PARAMS["Mock-1000K"]["species"],
        seed   = DATASET_PARAMS["Mock-1000K"]["seed"],
        script = "../scripts/generate_mock_clust_data.py",
    conda: "../envs/seqkit-biopython.yml"
    log: os.path.join(logdir, "mock_Mock-1000K.log")
    benchmark: os.path.join(benchmarks, "mock_Mock-1000K.benchmark")
    threads: get_resources("Mock-1000K")["threads"]
    resources:
        mem_mb = get_resources("Mock-1000K")["mem_mb"],
        disk_mb = get_resources("Mock-1000K")["disk_mb"]
    shell:
        """
        set -euo pipefail
        rm -rf {params.tmpdir}
        mkdir -p {params.tmpdir}

        python {params.script} \
            --name {params.name} \
            --num-sequences {params.num} \
            --outdir {datadir}/mock-data \
            --tmpdir {params.tmpdir} \
            --viral-seq {input.vir} \
            --prokaryotic-seq {input.prok} \
            --eukaryotic-seq {input.euk} \
            --virus-frac {params.virus} \
            --prokaryote-frac {params.prok_f} \
            --eukaryote-frac {params.euk_f} \
            --strain-mode {params.strain} \
            --num-species {params.species} \
            --mut-rate-min 0.001 \
            --mut-rate-max 0.05 \
            --seed {params.seed} \
            --force \
            > {log} 2>&1

        touch {output.fna} {output.gt}
        """

rule mock_Mock_Strain:
    name: "cluster-benchmark.smk Mock-Strain dataset"
    output:
        fna = os.path.join(datadir, "mock-data", "Mock-Strain.fna"),
        gt  = os.path.join(datadir, "mock-data", "Mock-Strain.ground_truth.tsv")
    input:
        vir  = os.path.join(datadir, "mock-data", "genomes", "refseq_viral.fna"),
        prok = os.path.join(datadir, "mock-data", "genomes", "refseq_prok.fna"),
        euk  = os.path.join(datadir, "mock-data", "genomes", "refseq_euk.fna")
    params:
        tmpdir = os.path.join(tmpd, "Mock-Strain"),
        name   = "Mock-Strain",
        num    = DATASET_PARAMS["Mock-Strain"]["size"],
        virus  = DATASET_PARAMS["Mock-Strain"]["virus"],
        prok_f = DATASET_PARAMS["Mock-Strain"]["prok"],
        euk_f  = DATASET_PARAMS["Mock-Strain"]["euk"],
        strain = DATASET_PARAMS["Mock-Strain"]["strain"],
        species = DATASET_PARAMS["Mock-Strain"]["species"],
        seed   = DATASET_PARAMS["Mock-Strain"]["seed"],
        script = "../scripts/generate_mock_clust_data.py",
    conda: "../envs/seqkit-biopython.yml"
    log: os.path.join(logdir, "mock_Mock-Strain.log")
    benchmark: os.path.join(benchmarks, "mock_Mock-Strain.benchmark")
    threads: get_resources("Mock-Strain")["threads"]
    resources:
        mem_mb = get_resources("Mock-Strain")["mem_mb"],
        disk_mb = get_resources("Mock-Strain")["disk_mb"]
    shell:
        """
        set -euo pipefail
        rm -rf {params.tmpdir}
        mkdir -p {params.tmpdir}

        python {params.script} \
            --name {params.name} \
            --num-sequences {params.num} \
            --outdir {datadir}/mock-data \
            --tmpdir {params.tmpdir} \
            --viral-seq {input.vir} \
            --prokaryotic-seq {input.prok} \
            --eukaryotic-seq {input.euk} \
            --virus-frac {params.virus} \
            --prokaryote-frac {params.prok_f} \
            --eukaryote-frac {params.euk_f} \
            --strain-mode {params.strain} \
            --num-species {params.species} \
            --mut-rate-min 0.001 \
            --mut-rate-max 0.05 \
            --seed {params.seed} \
            --force \
            > {log} 2>&1

        touch {output.fna} {output.gt}
        """

# rule generate_mock:
#     name: "cluster-benchmark.smk Generating mock datasets"
#     output:
#         fna = os.path.join(datadir, "mock-data", "{dataset}.fna"),
#         gt  = os.path.join(datadir, "mock-data", "{dataset}.ground_truth.tsv")
#     input:
#         vir  = os.path.join(datadir, "mock-data", "genomes", "refseq_viral.fna"),
#         prok = os.path.join(datadir, "mock-data", "genomes", "refseq_prok.fna"),
#         euk  = os.path.join(datadir, "mock-data", "genomes", "refseq_euk.fna")
#     params:
#         tmpdir = os.path.join(tmpd, "{dataset}"),
#         name   = lambda wc: wc.dataset,
#         num    = lambda wc: DATASET_PARAMS[wc.dataset]["size"],
#         virus  = lambda wc: DATASET_PARAMS[wc.dataset]["virus"],
#         prok_f = lambda wc: DATASET_PARAMS[wc.dataset]["prok"],
#         euk_f  = lambda wc: DATASET_PARAMS[wc.dataset]["euk"],
#         strain = lambda wc: DATASET_PARAMS[wc.dataset]["strain"],
#         species = lambda wc: DATASET_PARAMS[wc.dataset]["species"],
#         seed   = lambda wc: DATASET_PARAMS[wc.dataset]["seed"],
#         script = "../scripts/generate_mock_clust_data.py",
#     conda: "../envs/seqkit-biopython.yml"
#     log: os.path.join(logdir, "mock_{dataset}.log")
#     benchmark: os.path.join(benchmarks, "mock_{dataset}.benchmark")
#     threads: lambda wc: get_resources(wc)["threads"]
#     resources:
#         mem_mb = lambda wc: get_resources(wc)["mem_mb"],
#         disk_mb = lambda wc: get_resources(wc)["disk_mb"]
#     shell:
#         """
#         set -euo pipefail
#         rm -rf {params.tmpdir}
#         mkdir -p {params.tmpdir}

#         python {params.script} \
#             --name {params.name} \
#             --num-sequences {params.num} \
#             --outdir {datadir}/mock-data \
#             --tmpdir {params.tmpdir} \
#             --viral-seq {input.vir} \
#             --prokaryotic-seq {input.prok} \
#             --eukaryotic-seq {input.euk} \
#             --virus-frac {params.virus} \
#             --prokaryote-frac {params.prok_f} \
#             --eukaryote-frac {params.euk_f} \
#             --strain-mode {params.strain} \
#             --num-species {params.species} \
#             --mut-rate-min 0.001 \
#             --mut-rate-max 0.05 \
#             --seed {params.seed} \
#             --force \
#             > {log} 2>&1

#         touch {output.fna} {output.gt}
#         """


# rule cami_marine:
#     name: "cluster-benchmark.smk - Download CAMI Marine"
#     output:
#         marker = os.path.join(basedir, "cami_marine", "README.txt")
#     params:
#         outdir   = os.path.join(basedir, "cami_marine"),
#         tmpdir   = os.path.join(tmpd, "cami_marine"),
#         url      = config.get("cami_marine_url", ""),
#     conda: "../envs/seqkit-biopython.yml"
#     log: os.path.join(logdir, "mock_cami_marine.log")
#     benchmark: os.path.join(benchmarks, "mock_cami_marine.benchmark")
#     threads: 4
#     resources: {mem_mb: 4096, disk_mb: 10000}
#     shell:
#         """
#         set -euo pipefail
#         rm -rf {params.tmpdir}
#         mkdir -p {params.tmpdir} {params.outdir}

#         python ../scripts/download_cami.py \
#             --outdir {params.outdir} \
#             --tmpdir {params.tmpdir} \
#             --url "{params.url}" \
#             > {log} 2>&1

#         touch {output.marker}
#         """


# # ================================================================
# # 6. Aggregate rule – signals all mock data is ready
# # ================================================================
# rule all_mock_done:
#     name: "cluster-benchmark.smk - All mock datasets ready"
#     localrule: True
#     input:
#         # Tier 1
#         rules.mock_10K.output.fna,
#         rules.mock_50K.output.fna,
#         rules.mock_100K.output.fna,
#         rules.mock_300K.output.fna,
#         rules.mock_Strain.output.fna,
#         # Tier 2
#         rules.cami_marine.output.marker,
#     output:
#         done = os.path.join(benchmarks, "mock_all.done")
#     shell:
#         "touch {output.done}"