# ----------------------------------------------------------------------
# Configuration & setup
# ----------------------------------------------------------------------
vomix_module = vomix_utils.current_module
logdir = vomix_module.logdir
benchmarks = vomix_module.benchmarks
tmpd = vomix_module.tmpd

# ------------------------------------------------------------
# Setup mock dataset parameters – all parameters explicit per dataset
# ------------------------------------------------------------

DATASET_PARAMS = {
    "Mock-10K": {
        "total_sequences": 10000,
        "virus": 0.5,
        "prok": 0.3,
        "euk": 0.2,
        "seed": config.get("seed", 42),
        "fragments_per_10kb": 1.0,
        "min_fragments_per_genome": 10,
        "max_fragments_per_genome": 500,
        "copies_mean": 3.0,
        "copies_min": 2,
        "copies_max": 20,
        "copies_distribution": "poisson",
        "mut_rate_min": 0.001,   # 99.9% ANI
        "mut_rate_max": 0.04,    # 93% ANI (includes MIUViG boundary)
        "indel_rate": 0.1,
        "lognormal_mu": 8.5,
        "lognormal_sigma": 1.2,
        "min_len": 500,
        "max_len": 50000,
    },
    "Mock-10K-HighVir": {
        "total_sequences": 10000,
        "virus": 1.0,
        "prok": 0.0,
        "euk": 0.0,
        "seed": config.get("seed", 42),
        "fragments_per_10kb": 1.0,
        "min_fragments_per_genome": 10,
        "max_fragments_per_genome": 500,
        "copies_mean": 3.0,
        "copies_min": 2,
        "copies_max": 20,
        "copies_distribution": "poisson",
        "mut_rate_min": 0.001,   # 99.9% ANI
        "mut_rate_max": 0.04,    # 93% ANI
        "indel_rate": 0.1,
        "lognormal_mu": 8.5,
        "lognormal_sigma": 1.2,
        "min_len": 500,
        "max_len": 50000,
    },
    "Mock-10K-LowVir": {
        "total_sequences": 10000,
        "virus": 0.1,
        "prok": 0.5,
        "euk": 0.4,
        "seed": config.get("seed", 42),
        "fragments_per_10kb": 1.0,
        "min_fragments_per_genome": 10,
        "max_fragments_per_genome": 500,
        "copies_mean": 3.0,
        "copies_min": 2,
        "copies_max": 20,
        "copies_distribution": "poisson",
        "mut_rate_min": 0.001,   # 99.9% ANI
        "mut_rate_max": 0.04,    # 93% ANI
        "indel_rate": 0.1,
        "lognormal_mu": 8.5,
        "lognormal_sigma": 1.2,
        "min_len": 500,
        "max_len": 50000,
    },
    "Mock-50K": {
        "total_sequences": 50000,
        "virus": 0.5,
        "prok": 0.3,
        "euk": 0.2,
        "seed": config.get("seed", 42),
        "fragments_per_10kb": 1.0,
        "min_fragments_per_genome": 10,
        "max_fragments_per_genome": 500,
        "copies_mean": 3.0,
        "copies_min": 2,
        "copies_max": 20,
        "copies_distribution": "poisson",
        "mut_rate_min": 0.001,   # 99.9% ANI
        "mut_rate_max": 0.04,    # 93% ANI
        "indel_rate": 0.1,
        "lognormal_mu": 8.5,
        "lognormal_sigma": 1.2,
        "min_len": 500,
        "max_len": 50000,
    },
    "Mock-100K": {
        "total_sequences": 100000,
        "virus": 0.5,
        "prok": 0.3,
        "euk": 0.2,
        "seed": config.get("seed", 42),
        "fragments_per_10kb": 1.0,
        "min_fragments_per_genome": 10,
        "max_fragments_per_genome": 500,
        "copies_mean": 3.0,
        "copies_min": 2,
        "copies_max": 20,
        "copies_distribution": "poisson",
        "mut_rate_min": 0.001,   # 99.9% ANI
        "mut_rate_max": 0.04,    # 93% ANI
        "indel_rate": 0.1,
        "lognormal_mu": 8.5,
        "lognormal_sigma": 1.2,
        "min_len": 500,
        "max_len": 50000,
    },
    "Mock-300K": {
        "total_sequences": 300000,
        "virus": 0.5,
        "prok": 0.3,
        "euk": 0.2,
        "seed": config.get("seed", 42),
        "fragments_per_10kb": 1.0,
        "min_fragments_per_genome": 10,
        "max_fragments_per_genome": 500,
        "copies_mean": 3.0,
        "copies_min": 2,
        "copies_max": 20,
        "copies_distribution": "poisson",
        "mut_rate_min": 0.001,   # 99.9% ANI
        "mut_rate_max": 0.04,    # 93% ANI
        "indel_rate": 0.1,
        "lognormal_mu": 8.5,
        "lognormal_sigma": 1.2,
        "min_len": 500,
        "max_len": 50000,
    },
    "Mock-1000K": {
        "total_sequences": 1000000,
        "virus": 0.5,
        "prok": 0.3,
        "euk": 0.2,
        "seed": config.get("seed", 42),
        "fragments_per_10kb": 1.0,
        "min_fragments_per_genome": 10,
        "max_fragments_per_genome": 500,
        "copies_mean": 3.0,
        "copies_min": 2,
        "copies_max": 20,
        "copies_distribution": "poisson",
        "mut_rate_min": 0.001,   # 99.9% ANI
        "mut_rate_max": 0.04,    # 93% ANI
        "indel_rate": 0.1,
        "lognormal_mu": 8.5,
        "lognormal_sigma": 1.2,
        "min_len": 500,
        "max_len": 50000,
    },
    "Mock-Strain": {
        "total_sequences": 20000,
        "virus": 1.0,
        "prok": 0.0,
        "euk": 0.0,
        "seed": config.get("seed", 42),
        "fragments_per_10kb": 1.0,
        "min_fragments_per_genome": 10,
        "max_fragments_per_genome": 500,
        "copies_mean": 5.0,
        "copies_min": 4,
        "copies_max": 8,
        "copies_distribution": "poisson",
        "mut_rate_min": 0.001,   # 99.9% ANI
        "mut_rate_max": 0.04,    # 93% ANI
        "indel_rate": 0.1,
        "lognormal_mu": 8.5,
        "lognormal_sigma": 1.2,
        "min_len": 500,
        "max_len": 50000,
    },
    # NEW: Genus‑level dataset (70% ANI threshold)
    "Mock-Genus": {
        "total_sequences": 20000,
        "virus": 1.0,
        "prok": 0.0,
        "euk": 0.0,
        "seed": config.get("seed", 42),
        "fragments_per_10kb": 1.0,
        "min_fragments_per_genome": 10,
        "max_fragments_per_genome": 500,
        "copies_mean": 3.0,
        "copies_min": 2,
        "copies_max": 20,
        "copies_distribution": "poisson",
        "mut_rate_min": 0.001,   # 99.9% ANI
        "mut_rate_max": 0.15,    # 70% ANI (ICTV genus threshold)
        "indel_rate": 0.1,
        "lognormal_mu": 8.5,
        "lognormal_sigma": 1.2,
        "min_len": 500,
        "max_len": 50000,
    },
}

def get_resources(dataset):
    """Return memory, disk, and threads based on dataset total_sequences."""
    size = DATASET_PARAMS[dataset]["total_sequences"]
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
        script   = "workflow/scripts/download_refseq_genomes.py",
        tmpdir   = os.path.join(tmpd, "genomes", "viral"),
        outdir   = os.path.join(datadir, "mock-data", "genomes"),
    log: os.path.join(logdir, "download_refseq_viral.log")
    conda: "../envs/seqkit-biopython.yml"
    resources:
        mem_mb=lambda wildcards, attempt, input: 4 * 10**3 * attempt
    shell:
        """
        rm -rf {params.tmpdir}
        mkdir -p {params.tmpdir} {params.outdir}

        python {params.script} \
            --outfile {params.tmpdir}/tmp.fna \
            --summary-tsv {params.tmpdir}/tmp.tsv \
            --email "{params.email}" \
            --api-key "{params.api_key}" \
            --mode viral \
            --genomes-per-category 0 \
            --verbose \
            &> {log}

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
        script   = "workflow/scripts/download_refseq_genomes.py",
        tmpdir   = os.path.join(tmpd, "genomes", "prok"),
        outdir   = os.path.join(datadir, "mock-data", "genomes"),
    log: os.path.join(logdir, "download_refseq_prok.log")
    conda: "../envs/seqkit-biopython.yml"
    resources:
        mem_mb=lambda wildcards, attempt, input: 4 * 10**3 * attempt
    shell:
        """
        rm -rf {params.tmpdir}
        mkdir -p {params.tmpdir} {params.outdir}

        python {params.script} \
            --outfile {params.tmpdir}/tmp.fna \
            --summary-tsv {params.tmpdir}/tmp.tsv \
            --email "{params.email}" \
            --api-key "{params.api_key}" \
            --mode prokaryotic \
            --verbose \
            &> {log}

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
        script   = "workflow/scripts/download_refseq_genomes.py",
        tmpdir   = os.path.join(tmpd, "genomes", "euk"),
        outdir   = os.path.join(datadir, "mock-data", "genomes"),
    log: os.path.join(logdir, "download_refseq_euk.log")
    conda: "../envs/seqkit-biopython.yml"
    resources:
        mem_mb=lambda wildcards, attempt, input: 4 * 10**3 * attempt
    shell:
        """
        rm -rf {params.tmpdir}
        mkdir -p {params.tmpdir} {params.outdir}

        python {params.script} \
            --outfile {params.tmpdir}/tmp.fna \
            --summary-tsv {params.tmpdir}/tmp.tsv \
            --email "{params.email}" \
            --api-key "{params.api_key}" \
            --mode eukaryotic \
            --verbose \
            &> {log}

        mv {params.tmpdir}/tmp.fna {output.fna}
        mv {params.tmpdir}/tmp.tsv {output.tsv}
        """


# ------------------------------------------------------------
# Mock dataset generation – explicit rules for each dataset
# ------------------------------------------------------------

# ----- Mock-10K -----
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
        total_sequences = DATASET_PARAMS["Mock-10K"]["total_sequences"],
        virus_frac = DATASET_PARAMS["Mock-10K"]["virus"],
        prok_frac = DATASET_PARAMS["Mock-10K"]["prok"],
        euk_frac = DATASET_PARAMS["Mock-10K"]["euk"],
        seed = DATASET_PARAMS["Mock-10K"]["seed"],
        fragments_per_10kb = DATASET_PARAMS["Mock-10K"]["fragments_per_10kb"],
        min_fragments_per_genome = DATASET_PARAMS["Mock-10K"]["min_fragments_per_genome"],
        max_fragments_per_genome = DATASET_PARAMS["Mock-10K"]["max_fragments_per_genome"],
        copies_mean = DATASET_PARAMS["Mock-10K"]["copies_mean"],
        copies_min = DATASET_PARAMS["Mock-10K"]["copies_min"],
        copies_max = DATASET_PARAMS["Mock-10K"]["copies_max"],
        copies_distribution = DATASET_PARAMS["Mock-10K"]["copies_distribution"],
        mut_rate_min = DATASET_PARAMS["Mock-10K"]["mut_rate_min"],
        mut_rate_max = DATASET_PARAMS["Mock-10K"]["mut_rate_max"],
        indel_rate = DATASET_PARAMS["Mock-10K"]["indel_rate"],
        lognormal_mu = DATASET_PARAMS["Mock-10K"]["lognormal_mu"],
        lognormal_sigma = DATASET_PARAMS["Mock-10K"]["lognormal_sigma"],
        min_len = DATASET_PARAMS["Mock-10K"]["min_len"],
        max_len = DATASET_PARAMS["Mock-10K"]["max_len"],
        script = "workflow/scripts/generate_mock_clust_data.py",
        outdir = os.path.join(datadir, "mock-data"),
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
        mkdir -p {params.tmpdir} {params.outdir}

        python {params.script} \
            --name {params.name} \
            --total-sequences {params.total_sequences} \
            --outdir {params.outdir} \
            --viral-seq {input.vir} \
            --prokaryotic-seq {input.prok} \
            --eukaryotic-seq {input.euk} \
            --virus-frac {params.virus_frac} \
            --prokaryote-frac {params.prok_frac} \
            --eukaryote-frac {params.euk_frac} \
            --fragments-per-10kb {params.fragments_per_10kb} \
            --min-fragments-per-genome {params.min_fragments_per_genome} \
            --max-fragments-per-genome {params.max_fragments_per_genome} \
            --copies-mean {params.copies_mean} \
            --copies-min {params.copies_min} \
            --copies-max {params.copies_max} \
            --copies-distribution {params.copies_distribution} \
            --mut-rate-min {params.mut_rate_min} \
            --mut-rate-max {params.mut_rate_max} \
            --indel-rate {params.indel_rate} \
            --lognormal-mu {params.lognormal_mu} \
            --lognormal-sigma {params.lognormal_sigma} \
            --min-len {params.min_len} \
            --max-len {params.max_len} \
            --seed {params.seed} \
            --verbose \
            --force \
            &> {log}
        """

# ----- Mock-10K-HighVir -----
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
        total_sequences = DATASET_PARAMS["Mock-10K-HighVir"]["total_sequences"],
        virus_frac = DATASET_PARAMS["Mock-10K-HighVir"]["virus"],
        prok_frac = DATASET_PARAMS["Mock-10K-HighVir"]["prok"],
        euk_frac = DATASET_PARAMS["Mock-10K-HighVir"]["euk"],
        seed = DATASET_PARAMS["Mock-10K-HighVir"]["seed"],
        fragments_per_10kb = DATASET_PARAMS["Mock-10K-HighVir"]["fragments_per_10kb"],
        min_fragments_per_genome = DATASET_PARAMS["Mock-10K-HighVir"]["min_fragments_per_genome"],
        max_fragments_per_genome = DATASET_PARAMS["Mock-10K-HighVir"]["max_fragments_per_genome"],
        copies_mean = DATASET_PARAMS["Mock-10K-HighVir"]["copies_mean"],
        copies_min = DATASET_PARAMS["Mock-10K-HighVir"]["copies_min"],
        copies_max = DATASET_PARAMS["Mock-10K-HighVir"]["copies_max"],
        copies_distribution = DATASET_PARAMS["Mock-10K-HighVir"]["copies_distribution"],
        mut_rate_min = DATASET_PARAMS["Mock-10K-HighVir"]["mut_rate_min"],
        mut_rate_max = DATASET_PARAMS["Mock-10K-HighVir"]["mut_rate_max"],
        indel_rate = DATASET_PARAMS["Mock-10K-HighVir"]["indel_rate"],
        lognormal_mu = DATASET_PARAMS["Mock-10K-HighVir"]["lognormal_mu"],
        lognormal_sigma = DATASET_PARAMS["Mock-10K-HighVir"]["lognormal_sigma"],
        min_len = DATASET_PARAMS["Mock-10K-HighVir"]["min_len"],
        max_len = DATASET_PARAMS["Mock-10K-HighVir"]["max_len"],
        script = "workflow/scripts/generate_mock_clust_data.py",
        outdir = os.path.join(datadir, "mock-data"),
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
        mkdir -p {params.tmpdir} {params.outdir}

        python {params.script} \
            --name {params.name} \
            --total-sequences {params.total_sequences} \
            --outdir {params.outdir} \
            --viral-seq {input.vir} \
            --prokaryotic-seq {input.prok} \
            --eukaryotic-seq {input.euk} \
            --virus-frac {params.virus_frac} \
            --prokaryote-frac {params.prok_frac} \
            --eukaryote-frac {params.euk_frac} \
            --fragments-per-10kb {params.fragments_per_10kb} \
            --min-fragments-per-genome {params.min_fragments_per_genome} \
            --max-fragments-per-genome {params.max_fragments_per_genome} \
            --copies-mean {params.copies_mean} \
            --copies-min {params.copies_min} \
            --copies-max {params.copies_max} \
            --copies-distribution {params.copies_distribution} \
            --mut-rate-min {params.mut_rate_min} \
            --mut-rate-max {params.mut_rate_max} \
            --indel-rate {params.indel_rate} \
            --lognormal-mu {params.lognormal_mu} \
            --lognormal-sigma {params.lognormal_sigma} \
            --min-len {params.min_len} \
            --max-len {params.max_len} \
            --seed {params.seed} \
            --verbose \
            --force \
            &> {log}
        """

# ----- Mock-10K-LowVir -----
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
        total_sequences = DATASET_PARAMS["Mock-10K-LowVir"]["total_sequences"],
        virus_frac = DATASET_PARAMS["Mock-10K-LowVir"]["virus"],
        prok_frac = DATASET_PARAMS["Mock-10K-LowVir"]["prok"],
        euk_frac = DATASET_PARAMS["Mock-10K-LowVir"]["euk"],
        seed = DATASET_PARAMS["Mock-10K-LowVir"]["seed"],
        fragments_per_10kb = DATASET_PARAMS["Mock-10K-LowVir"]["fragments_per_10kb"],
        min_fragments_per_genome = DATASET_PARAMS["Mock-10K-LowVir"]["min_fragments_per_genome"],
        max_fragments_per_genome = DATASET_PARAMS["Mock-10K-LowVir"]["max_fragments_per_genome"],
        copies_mean = DATASET_PARAMS["Mock-10K-LowVir"]["copies_mean"],
        copies_min = DATASET_PARAMS["Mock-10K-LowVir"]["copies_min"],
        copies_max = DATASET_PARAMS["Mock-10K-LowVir"]["copies_max"],
        copies_distribution = DATASET_PARAMS["Mock-10K-LowVir"]["copies_distribution"],
        mut_rate_min = DATASET_PARAMS["Mock-10K-LowVir"]["mut_rate_min"],
        mut_rate_max = DATASET_PARAMS["Mock-10K-LowVir"]["mut_rate_max"],
        indel_rate = DATASET_PARAMS["Mock-10K-LowVir"]["indel_rate"],
        lognormal_mu = DATASET_PARAMS["Mock-10K-LowVir"]["lognormal_mu"],
        lognormal_sigma = DATASET_PARAMS["Mock-10K-LowVir"]["lognormal_sigma"],
        min_len = DATASET_PARAMS["Mock-10K-LowVir"]["min_len"],
        max_len = DATASET_PARAMS["Mock-10K-LowVir"]["max_len"],
        script = "workflow/scripts/generate_mock_clust_data.py",
        outdir = os.path.join(datadir, "mock-data"),
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
        mkdir -p {params.tmpdir} {params.outdir}

        python {params.script} \
            --name {params.name} \
            --total-sequences {params.total_sequences} \
            --outdir {params.outdir} \
            --viral-seq {input.vir} \
            --prokaryotic-seq {input.prok} \
            --eukaryotic-seq {input.euk} \
            --virus-frac {params.virus_frac} \
            --prokaryote-frac {params.prok_frac} \
            --eukaryote-frac {params.euk_frac} \
            --fragments-per-10kb {params.fragments_per_10kb} \
            --min-fragments-per-genome {params.min_fragments_per_genome} \
            --max-fragments-per-genome {params.max_fragments_per_genome} \
            --copies-mean {params.copies_mean} \
            --copies-min {params.copies_min} \
            --copies-max {params.copies_max} \
            --copies-distribution {params.copies_distribution} \
            --mut-rate-min {params.mut_rate_min} \
            --mut-rate-max {params.mut_rate_max} \
            --indel-rate {params.indel_rate} \
            --lognormal-mu {params.lognormal_mu} \
            --lognormal-sigma {params.lognormal_sigma} \
            --min-len {params.min_len} \
            --max-len {params.max_len} \
            --seed {params.seed} \
            --verbose \
            --force \
            &> {log}
        """

# ----- Mock-50K -----
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
        total_sequences = DATASET_PARAMS["Mock-50K"]["total_sequences"],
        virus_frac = DATASET_PARAMS["Mock-50K"]["virus"],
        prok_frac = DATASET_PARAMS["Mock-50K"]["prok"],
        euk_frac = DATASET_PARAMS["Mock-50K"]["euk"],
        seed = DATASET_PARAMS["Mock-50K"]["seed"],
        fragments_per_10kb = DATASET_PARAMS["Mock-50K"]["fragments_per_10kb"],
        min_fragments_per_genome = DATASET_PARAMS["Mock-50K"]["min_fragments_per_genome"],
        max_fragments_per_genome = DATASET_PARAMS["Mock-50K"]["max_fragments_per_genome"],
        copies_mean = DATASET_PARAMS["Mock-50K"]["copies_mean"],
        copies_min = DATASET_PARAMS["Mock-50K"]["copies_min"],
        copies_max = DATASET_PARAMS["Mock-50K"]["copies_max"],
        copies_distribution = DATASET_PARAMS["Mock-50K"]["copies_distribution"],
        mut_rate_min = DATASET_PARAMS["Mock-50K"]["mut_rate_min"],
        mut_rate_max = DATASET_PARAMS["Mock-50K"]["mut_rate_max"],
        indel_rate = DATASET_PARAMS["Mock-50K"]["indel_rate"],
        lognormal_mu = DATASET_PARAMS["Mock-50K"]["lognormal_mu"],
        lognormal_sigma = DATASET_PARAMS["Mock-50K"]["lognormal_sigma"],
        min_len = DATASET_PARAMS["Mock-50K"]["min_len"],
        max_len = DATASET_PARAMS["Mock-50K"]["max_len"],
        script = "workflow/scripts/generate_mock_clust_data.py",
        outdir = os.path.join(datadir, "mock-data"),
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
        mkdir -p {params.tmpdir} {params.outdir}

        python {params.script} \
            --name {params.name} \
            --total-sequences {params.total_sequences} \
            --outdir {params.outdir} \
            --viral-seq {input.vir} \
            --prokaryotic-seq {input.prok} \
            --eukaryotic-seq {input.euk} \
            --virus-frac {params.virus_frac} \
            --prokaryote-frac {params.prok_frac} \
            --eukaryote-frac {params.euk_frac} \
            --fragments-per-10kb {params.fragments_per_10kb} \
            --min-fragments-per-genome {params.min_fragments_per_genome} \
            --max-fragments-per-genome {params.max_fragments_per_genome} \
            --copies-mean {params.copies_mean} \
            --copies-min {params.copies_min} \
            --copies-max {params.copies_max} \
            --copies-distribution {params.copies_distribution} \
            --mut-rate-min {params.mut_rate_min} \
            --mut-rate-max {params.mut_rate_max} \
            --indel-rate {params.indel_rate} \
            --lognormal-mu {params.lognormal_mu} \
            --lognormal-sigma {params.lognormal_sigma} \
            --min-len {params.min_len} \
            --max-len {params.max_len} \
            --seed {params.seed} \
            --verbose \
            --force \
            &> {log}
        """

# ----- Mock-100K -----
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
        total_sequences = DATASET_PARAMS["Mock-100K"]["total_sequences"],
        virus_frac = DATASET_PARAMS["Mock-100K"]["virus"],
        prok_frac = DATASET_PARAMS["Mock-100K"]["prok"],
        euk_frac = DATASET_PARAMS["Mock-100K"]["euk"],
        seed = DATASET_PARAMS["Mock-100K"]["seed"],
        fragments_per_10kb = DATASET_PARAMS["Mock-100K"]["fragments_per_10kb"],
        min_fragments_per_genome = DATASET_PARAMS["Mock-100K"]["min_fragments_per_genome"],
        max_fragments_per_genome = DATASET_PARAMS["Mock-100K"]["max_fragments_per_genome"],
        copies_mean = DATASET_PARAMS["Mock-100K"]["copies_mean"],
        copies_min = DATASET_PARAMS["Mock-100K"]["copies_min"],
        copies_max = DATASET_PARAMS["Mock-100K"]["copies_max"],
        copies_distribution = DATASET_PARAMS["Mock-100K"]["copies_distribution"],
        mut_rate_min = DATASET_PARAMS["Mock-100K"]["mut_rate_min"],
        mut_rate_max = DATASET_PARAMS["Mock-100K"]["mut_rate_max"],
        indel_rate = DATASET_PARAMS["Mock-100K"]["indel_rate"],
        lognormal_mu = DATASET_PARAMS["Mock-100K"]["lognormal_mu"],
        lognormal_sigma = DATASET_PARAMS["Mock-100K"]["lognormal_sigma"],
        min_len = DATASET_PARAMS["Mock-100K"]["min_len"],
        max_len = DATASET_PARAMS["Mock-100K"]["max_len"],
        script = "workflow/scripts/generate_mock_clust_data.py",
        outdir = os.path.join(datadir, "mock-data"),
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
        mkdir -p {params.tmpdir} {params.outdir}

        python {params.script} \
            --name {params.name} \
            --total-sequences {params.total_sequences} \
            --outdir {params.outdir} \
            --viral-seq {input.vir} \
            --prokaryotic-seq {input.prok} \
            --eukaryotic-seq {input.euk} \
            --virus-frac {params.virus_frac} \
            --prokaryote-frac {params.prok_frac} \
            --eukaryote-frac {params.euk_frac} \
            --fragments-per-10kb {params.fragments_per_10kb} \
            --min-fragments-per-genome {params.min_fragments_per_genome} \
            --max-fragments-per-genome {params.max_fragments_per_genome} \
            --copies-mean {params.copies_mean} \
            --copies-min {params.copies_min} \
            --copies-max {params.copies_max} \
            --copies-distribution {params.copies_distribution} \
            --mut-rate-min {params.mut_rate_min} \
            --mut-rate-max {params.mut_rate_max} \
            --indel-rate {params.indel_rate} \
            --lognormal-mu {params.lognormal_mu} \
            --lognormal-sigma {params.lognormal_sigma} \
            --min-len {params.min_len} \
            --max-len {params.max_len} \
            --seed {params.seed} \
            --verbose \
            --force \
            &> {log}
        """

# ----- Mock-300K -----
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
        total_sequences = DATASET_PARAMS["Mock-300K"]["total_sequences"],
        virus_frac = DATASET_PARAMS["Mock-300K"]["virus"],
        prok_frac = DATASET_PARAMS["Mock-300K"]["prok"],
        euk_frac = DATASET_PARAMS["Mock-300K"]["euk"],
        seed = DATASET_PARAMS["Mock-300K"]["seed"],
        fragments_per_10kb = DATASET_PARAMS["Mock-300K"]["fragments_per_10kb"],
        min_fragments_per_genome = DATASET_PARAMS["Mock-300K"]["min_fragments_per_genome"],
        max_fragments_per_genome = DATASET_PARAMS["Mock-300K"]["max_fragments_per_genome"],
        copies_mean = DATASET_PARAMS["Mock-300K"]["copies_mean"],
        copies_min = DATASET_PARAMS["Mock-300K"]["copies_min"],
        copies_max = DATASET_PARAMS["Mock-300K"]["copies_max"],
        copies_distribution = DATASET_PARAMS["Mock-300K"]["copies_distribution"],
        mut_rate_min = DATASET_PARAMS["Mock-300K"]["mut_rate_min"],
        mut_rate_max = DATASET_PARAMS["Mock-300K"]["mut_rate_max"],
        indel_rate = DATASET_PARAMS["Mock-300K"]["indel_rate"],
        lognormal_mu = DATASET_PARAMS["Mock-300K"]["lognormal_mu"],
        lognormal_sigma = DATASET_PARAMS["Mock-300K"]["lognormal_sigma"],
        min_len = DATASET_PARAMS["Mock-300K"]["min_len"],
        max_len = DATASET_PARAMS["Mock-300K"]["max_len"],
        script = "workflow/scripts/generate_mock_clust_data.py",
        outdir = os.path.join(datadir, "mock-data"),
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
        mkdir -p {params.tmpdir} {params.outdir}

        python {params.script} \
            --name {params.name} \
            --total-sequences {params.total_sequences} \
            --outdir {params.outdir} \
            --viral-seq {input.vir} \
            --prokaryotic-seq {input.prok} \
            --eukaryotic-seq {input.euk} \
            --virus-frac {params.virus_frac} \
            --prokaryote-frac {params.prok_frac} \
            --eukaryote-frac {params.euk_frac} \
            --fragments-per-10kb {params.fragments_per_10kb} \
            --min-fragments-per-genome {params.min_fragments_per_genome} \
            --max-fragments-per-genome {params.max_fragments_per_genome} \
            --copies-mean {params.copies_mean} \
            --copies-min {params.copies_min} \
            --copies-max {params.copies_max} \
            --copies-distribution {params.copies_distribution} \
            --mut-rate-min {params.mut_rate_min} \
            --mut-rate-max {params.mut_rate_max} \
            --indel-rate {params.indel_rate} \
            --lognormal-mu {params.lognormal_mu} \
            --lognormal-sigma {params.lognormal_sigma} \
            --min-len {params.min_len} \
            --max-len {params.max_len} \
            --seed {params.seed} \
            --verbose \
            --force \
            &> {log}
        """

# ----- Mock-1000K -----
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
        total_sequences = DATASET_PARAMS["Mock-1000K"]["total_sequences"],
        virus_frac = DATASET_PARAMS["Mock-1000K"]["virus"],
        prok_frac = DATASET_PARAMS["Mock-1000K"]["prok"],
        euk_frac = DATASET_PARAMS["Mock-1000K"]["euk"],
        seed = DATASET_PARAMS["Mock-1000K"]["seed"],
        fragments_per_10kb = DATASET_PARAMS["Mock-1000K"]["fragments_per_10kb"],
        min_fragments_per_genome = DATASET_PARAMS["Mock-1000K"]["min_fragments_per_genome"],
        max_fragments_per_genome = DATASET_PARAMS["Mock-1000K"]["max_fragments_per_genome"],
        copies_mean = DATASET_PARAMS["Mock-1000K"]["copies_mean"],
        copies_min = DATASET_PARAMS["Mock-1000K"]["copies_min"],
        copies_max = DATASET_PARAMS["Mock-1000K"]["copies_max"],
        copies_distribution = DATASET_PARAMS["Mock-1000K"]["copies_distribution"],
        mut_rate_min = DATASET_PARAMS["Mock-1000K"]["mut_rate_min"],
        mut_rate_max = DATASET_PARAMS["Mock-1000K"]["mut_rate_max"],
        indel_rate = DATASET_PARAMS["Mock-1000K"]["indel_rate"],
        lognormal_mu = DATASET_PARAMS["Mock-1000K"]["lognormal_mu"],
        lognormal_sigma = DATASET_PARAMS["Mock-1000K"]["lognormal_sigma"],
        min_len = DATASET_PARAMS["Mock-1000K"]["min_len"],
        max_len = DATASET_PARAMS["Mock-1000K"]["max_len"],
        script = "workflow/scripts/generate_mock_clust_data.py",
        outdir = os.path.join(datadir, "mock-data"),
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
        mkdir -p {params.tmpdir} {params.outdir}

        python {params.script} \
            --name {params.name} \
            --total-sequences {params.total_sequences} \
            --outdir {params.outdir} \
            --viral-seq {input.vir} \
            --prokaryotic-seq {input.prok} \
            --eukaryotic-seq {input.euk} \
            --virus-frac {params.virus_frac} \
            --prokaryote-frac {params.prok_frac} \
            --eukaryote-frac {params.euk_frac} \
            --fragments-per-10kb {params.fragments_per_10kb} \
            --min-fragments-per-genome {params.min_fragments_per_genome} \
            --max-fragments-per-genome {params.max_fragments_per_genome} \
            --copies-mean {params.copies_mean} \
            --copies-min {params.copies_min} \
            --copies-max {params.copies_max} \
            --copies-distribution {params.copies_distribution} \
            --mut-rate-min {params.mut_rate_min} \
            --mut-rate-max {params.mut_rate_max} \
            --indel-rate {params.indel_rate} \
            --lognormal-mu {params.lognormal_mu} \
            --lognormal-sigma {params.lognormal_sigma} \
            --min-len {params.min_len} \
            --max-len {params.max_len} \
            --seed {params.seed} \
            --verbose \
            --force \
            &> {log}
        """

# ----- Mock-Strain -----
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
        total_sequences = DATASET_PARAMS["Mock-Strain"]["total_sequences"],
        virus_frac = DATASET_PARAMS["Mock-Strain"]["virus"],
        prok_frac = DATASET_PARAMS["Mock-Strain"]["prok"],
        euk_frac = DATASET_PARAMS["Mock-Strain"]["euk"],
        seed = DATASET_PARAMS["Mock-Strain"]["seed"],
        fragments_per_10kb = DATASET_PARAMS["Mock-Strain"]["fragments_per_10kb"],
        min_fragments_per_genome = DATASET_PARAMS["Mock-Strain"]["min_fragments_per_genome"],
        max_fragments_per_genome = DATASET_PARAMS["Mock-Strain"]["max_fragments_per_genome"],
        copies_mean = DATASET_PARAMS["Mock-Strain"]["copies_mean"],
        copies_min = DATASET_PARAMS["Mock-Strain"]["copies_min"],
        copies_max = DATASET_PARAMS["Mock-Strain"]["copies_max"],
        copies_distribution = DATASET_PARAMS["Mock-Strain"]["copies_distribution"],
        mut_rate_min = DATASET_PARAMS["Mock-Strain"]["mut_rate_min"],
        mut_rate_max = DATASET_PARAMS["Mock-Strain"]["mut_rate_max"],
        indel_rate = DATASET_PARAMS["Mock-Strain"]["indel_rate"],
        lognormal_mu = DATASET_PARAMS["Mock-Strain"]["lognormal_mu"],
        lognormal_sigma = DATASET_PARAMS["Mock-Strain"]["lognormal_sigma"],
        min_len = DATASET_PARAMS["Mock-Strain"]["min_len"],
        max_len = DATASET_PARAMS["Mock-Strain"]["max_len"],
        script = "workflow/scripts/generate_mock_clust_data.py",
        outdir = os.path.join(datadir, "mock-data"),
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
        mkdir -p {params.tmpdir} {params.outdir}

        python {params.script} \
            --name {params.name} \
            --total-sequences {params.total_sequences} \
            --outdir {params.outdir} \
            --viral-seq {input.vir} \
            --prokaryotic-seq {input.prok} \
            --eukaryotic-seq {input.euk} \
            --virus-frac {params.virus_frac} \
            --prokaryote-frac {params.prok_frac} \
            --eukaryote-frac {params.euk_frac} \
            --fragments-per-10kb {params.fragments_per_10kb} \
            --min-fragments-per-genome {params.min_fragments_per_genome} \
            --max-fragments-per-genome {params.max_fragments_per_genome} \
            --copies-mean {params.copies_mean} \
            --copies-min {params.copies_min} \
            --copies-max {params.copies_max} \
            --copies-distribution {params.copies_distribution} \
            --mut-rate-min {params.mut_rate_min} \
            --mut-rate-max {params.mut_rate_max} \
            --indel-rate {params.indel_rate} \
            --lognormal-mu {params.lognormal_mu} \
            --lognormal-sigma {params.lognormal_sigma} \
            --min-len {params.min_len} \
            --max-len {params.max_len} \
            --seed {params.seed} \
            --verbose \
            --force \
            &> {log}
        """

# ----- Mock-Genus (70% ANI genus threshold) -----
rule mock_Mock_Genus:
    name: "cluster-benchmark.smk Mock-Genus dataset"
    output:
        fna = os.path.join(datadir, "mock-data", "Mock-Genus.fna"),
        gt  = os.path.join(datadir, "mock-data", "Mock-Genus.ground_truth.tsv")
    input:
        vir  = os.path.join(datadir, "mock-data", "genomes", "refseq_viral.fna"),
        prok = os.path.join(datadir, "mock-data", "genomes", "refseq_prok.fna"),
        euk  = os.path.join(datadir, "mock-data", "genomes", "refseq_euk.fna")
    params:
        tmpdir = os.path.join(tmpd, "Mock-Genus"),
        name   = "Mock-Genus",
        total_sequences = DATASET_PARAMS["Mock-Genus"]["total_sequences"],
        virus_frac = DATASET_PARAMS["Mock-Genus"]["virus"],
        prok_frac = DATASET_PARAMS["Mock-Genus"]["prok"],
        euk_frac = DATASET_PARAMS["Mock-Genus"]["euk"],
        seed = DATASET_PARAMS["Mock-Genus"]["seed"],
        fragments_per_10kb = DATASET_PARAMS["Mock-Genus"]["fragments_per_10kb"],
        min_fragments_per_genome = DATASET_PARAMS["Mock-Genus"]["min_fragments_per_genome"],
        max_fragments_per_genome = DATASET_PARAMS["Mock-Genus"]["max_fragments_per_genome"],
        copies_mean = DATASET_PARAMS["Mock-Genus"]["copies_mean"],
        copies_min = DATASET_PARAMS["Mock-Genus"]["copies_min"],
        copies_max = DATASET_PARAMS["Mock-Genus"]["copies_max"],
        copies_distribution = DATASET_PARAMS["Mock-Genus"]["copies_distribution"],
        mut_rate_min = DATASET_PARAMS["Mock-Genus"]["mut_rate_min"],
        mut_rate_max = DATASET_PARAMS["Mock-Genus"]["mut_rate_max"],
        indel_rate = DATASET_PARAMS["Mock-Genus"]["indel_rate"],
        lognormal_mu = DATASET_PARAMS["Mock-Genus"]["lognormal_mu"],
        lognormal_sigma = DATASET_PARAMS["Mock-Genus"]["lognormal_sigma"],
        min_len = DATASET_PARAMS["Mock-Genus"]["min_len"],
        max_len = DATASET_PARAMS["Mock-Genus"]["max_len"],
        script = "workflow/scripts/generate_mock_clust_data.py",
        outdir = os.path.join(datadir, "mock-data"),
    conda: "../envs/seqkit-biopython.yml"
    log: os.path.join(logdir, "mock_Mock-Genus.log")
    benchmark: os.path.join(benchmarks, "mock_Mock-Genus.benchmark")
    threads: get_resources("Mock-Genus")["threads"]
    resources:
        mem_mb = get_resources("Mock-Genus")["mem_mb"],
        disk_mb = get_resources("Mock-Genus")["disk_mb"]
    shell:
        """
        set -euo pipefail
        rm -rf {params.tmpdir}
        mkdir -p {params.tmpdir} {params.outdir}

        python {params.script} \
            --name {params.name} \
            --total-sequences {params.total_sequences} \
            --outdir {params.outdir} \
            --viral-seq {input.vir} \
            --prokaryotic-seq {input.prok} \
            --eukaryotic-seq {input.euk} \
            --virus-frac {params.virus_frac} \
            --prokaryote-frac {params.prok_frac} \
            --eukaryote-frac {params.euk_frac} \
            --fragments-per-10kb {params.fragments_per_10kb} \
            --min-fragments-per-genome {params.min_fragments_per_genome} \
            --max-fragments-per-genome {params.max_fragments_per_genome} \
            --copies-mean {params.copies_mean} \
            --copies-min {params.copies_min} \
            --copies-max {params.copies_max} \
            --copies-distribution {params.copies_distribution} \
            --mut-rate-min {params.mut_rate_min} \
            --mut-rate-max {params.mut_rate_max} \
            --indel-rate {params.indel_rate} \
            --lognormal-mu {params.lognormal_mu} \
            --lognormal-sigma {params.lognormal_sigma} \
            --min-len {params.min_len} \
            --max-len {params.max_len} \
            --seed {params.seed} \
            --verbose \
            --force \
            &> {log}
        """


# ------------------------------------------------------------
# CAMI Marine (commented out, kept for reference)
# ------------------------------------------------------------
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
#
#         python workflow/scripts/download_cami.py \
#             --outdir {params.outdir} \
#             --tmpdir {params.tmpdir} \
#             --url "{params.url}" \
#             &> {log}
#
#         touch {output.marker}
#         """