# vOMIX-snakemake --use-conda

vOMIX-snakemake --use-conda is a reproducible Snakemake workflow for viral metagenomics, including preprocessing, assembly, viral identification, taxonomy, host assignment, annotation, benchmarking, and end-to-end analysis.

## Quick install

### Conda

```bash
conda create -n vomix -c conda-forge snakemake conda=8.25.5 biopython=1.84 snakemake-executor-plugin-cluster-generic=1.0.9 -y
conda activate vomix
snakemake --use-conda -v
```

### Mamba

```bash
mamba create -n vomix -c conda-forge snakemake conda=8.25.5 biopython=1.84 snakemake-executor-plugin-cluster-generic=1.0.9 -y
mamba activate vomix
snakemake -v
```

### Clone the repository

```bash
git clone https://github.com/holab-hku/vomix-snakemake
cd vomix-snakemake
```

## Quick start

```bash
snakemake --use-conda --sdm conda --config module="viral-identify" outdir="test_res" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" -j 4 --latency-wait 20
```

## Documentation

The full documentation is available at:

<<https://vomix-snakemake> --use-conda.readthedocs.io/>

It includes installation guides, run configuration, module-specific workflows, and troubleshooting.
