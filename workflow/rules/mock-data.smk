import os

logdir     = os.path.join(config['datadir'], "mock-data/logs")
benchmarks = os.path.join(config['datadir'], "mock-data/benchmarks")
tmpd       = os.path.join(config['datadir'], "mock-data/tmp")
mock_out   = os.path.join(config['datadir'], "mock-data")

os.makedirs(logdir,     exist_ok=True)
os.makedirs(benchmarks, exist_ok=True)
os.makedirs(tmpd,       exist_ok=True)
os.makedirs(mock_out,   exist_ok=True)


rule done:
    name: "mock-data.smk Done. removing tmp files"
    localrule: True
    input:
        os.path.join(mock_out, "refseq_viral.fna"),
        os.path.join(mock_out, "refseq_contaminants.fna"),
        # os.path.join(mock_out, "Mock-10K.fna"),
        # os.path.join(mock_out, "Mock-50K.fna"),
        # os.path.join(mock_out, "Mock-100K.fna"),
        # os.path.join(mock_out, "Mock-300K.fna"),
        # os.path.join(mock_out, "Mock-Strain.fna"),
        # os.path.join(mock_out, "cami_marine", "README.txt"),
    output:
        os.path.join(logdir, "done.log")
    shell:
        "touch {output}"


rule download_refseq_viral:
    output:
        fna = os.path.join(mock_out, "refseq_viral.fna")
    params:
        email    = config.get("NCBI-email", ""),
        api_key  = config.get("NCBI-API-key", ""),
        script   = "../scripts/download_refseq_genomes.py",
        outdir   = mock_out,
        tmpdir   = os.path.join(tmpd, "refseq_viral"),
    conda: "../envs/seqkit-biopython.yml"
    shell:
        """
        rm -rf {params.tmpdir} {output.fna}
        mkdir -p {params.tmpdir} {params.outdir}

        python {params.script} \
            --outfile {params.tmpdir}/tmp.fna \
            --email "{params.email}" \
            --api-key "{params.api_key}" \
            --mode viral \
            --genomes-per-category 0 \
            > {log} 2>&1

        mv {params.tmpdir}/tmp.fna {output.fna}
        """

rule download_refseq_contaminants:
    output:
        fna = os.path.join(mock_out, "refseq_contaminants.fna")
    params:
        email    = config.get("NCBI-email", ""),
        api_key  = config.get("NCBI-API-key", ""),
        script   = "../scripts/download_refseq_genomes.py",
        outdir   = mock_out,
        tmpdir   = os.path.join(tmpd, "refseq_contaminants"),
    conda: "../envs/seqkit-biopython.yml"
    shell:
        """
        rm -rf {params.tmpdir}
        mkdir -p {params.tmpdir}

        python {params.script} \
            --outfile {params.tmpdir}/tmp.fna \
            --email "{params.email}" \
            --api-key "{params.api_key}" \
            --mode contaminants \
            > {log} 2>&1

        mv {params.tmpdir}/tmp.fna {output.fna}
        """

# rule mock_10K:
#     name: "mock-data.smk - Generate Mock-10K"
#     output:
#         fna = os.path.join(mock_out, "Mock-10K.fna")
#     input:
#         viral   = rules.download_refseq_viral.output.cache,
#         contam  = rules.download_refseq_contaminants.output.cache 
#     params:
#         outdir   = mock_out,
#         tmpdir   = os.path.join(tmpd, "Mock-10K"),
#         name     = "Mock-10K",
#         num      = 10000,
#         strain   = 0,
#         seed     = 42,
#     conda: "../envs/seqkit-biopython.yml"
#     log: os.path.join(logdir, "mock_Mock-10K.log")
#     benchmark: os.path.join(benchmarks, "mock_Mock-10K.benchmark")
#     threads: 4
#     resources:
#         mem_mb = 8192,
#     shell:
#         """
#         set -euo pipefail
#         rm -rf {params.tmpdir}
#         mkdir -p {params.tmpdir}

#         python ../scripts/generate_mock.py \
#             --name {params.name} \
#             --num-sequences {params.num} \
#             --outdir {params.outdir} \
#             --tmpdir {params.tmpdir} \
#             --viral-cache {input.viral} \
#             --contam-cache {input.contam} \
#             --strain-mode {params.strain} \
#             --seed {params.seed} \
#             > {log} 2>&1

#         touch {output.fna}
#         """

# # ---- Mock-50K ----
# rule mock_50K:
#     name: "mock-data.smk - Generate Mock-50K"
#     output: fna = os.path.join(mock_out, "Mock-50K.fna")
#     input:
#         viral   = rules.download_refseq_viral.output.cache,
#         contam  = rules.download_refseq_contaminants.output.cache
#     params:
#         outdir   = mock_out,
#         tmpdir   = os.path.join(tmpd, "Mock-50K"),
#         name     = "Mock-50K",
#         num      = 50000,
#         strain   = 0,
#         seed     = 123,
#         viral_cache = input.viral,
#         contam_cache = input.contam,
#     conda: "../envs/seqkit-biopython.yml"
#     log: os.path.join(logdir, "mock_Mock-50K.log")
#     benchmark: os.path.join(benchmarks, "mock_Mock-50K.benchmark")
#     threads: 4
#     resources: {mem_mb: 8192, disk_mb: 5000}
#     shell: """ ... same as above with appropriate params ... """

# # (Repeat for Mock-100K, Mock-300K with different num and seed)

# rule mock_100K:
#     output: fna = os.path.join(mock_out, "Mock-100K.fna")
#     params:
#         outdir   = mock_out,
#         tmpdir   = os.path.join(tmpd, "Mock-100K"),
#         name     = "Mock-100K",
#         num      = 100000,
#         strain   = 0,
#         seed     = 456,
#         viral_cache = input.viral,
#         contam_cache = input.contam,
#     # ... rest same

# rule mock_300K:
#     output: fna = os.path.join(mock_out, "Mock-300K.fna")
#     params:
#         outdir   = mock_out,
#         tmpdir   = os.path.join(tmpd, "Mock-300K"),
#         name     = "Mock-300K",
#         num      = 300000,
#         strain   = 0,
#         seed     = 789,
#         viral_cache = input.viral,
#         contam_cache = input.contam,
#     resources: {mem_mb: 16000, disk_mb: 20000}   # larger dataset


# rule mock_Strain:
#     output: fna = os.path.join(mock_out, "Mock-Strain.fna")
#     params:
#         outdir   = mock_out,
#         tmpdir   = os.path.join(tmpd, "Mock-Strain"),
#         name     = "Mock-Strain",
#         num      = 20000,
#         strain   = 1,      # enable strain‑level variation
#         seed     = 999,
#         viral_cache = input.viral,
#         contam_cache = input.contam,   # not really used for this dataset
#     resources: {mem_mb: 8192, disk_mb: 10000}


# rule cami_marine:
#     name: "mock-data.smk - Download CAMI Marine"
#     output:
#         marker = os.path.join(mock_out, "cami_marine", "README.txt")
#     params:
#         outdir   = os.path.join(mock_out, "cami_marine"),
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
#     name: "mock-data.smk - All mock datasets ready"
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