# `viral-taxonomy`

```{image} ../_static/pipeline_viral-taxonomy.svg
:width: 550
:class: no-scaled-link
:align: center
```

The `viral-taxonomy` module focuses on determining the taxonomic identity of each sequence in your vOTU database or user-provided viral FASTA sequences, leveraging a dual-engine framework: `geNomad` and `PhaGCN`. By default, the pipeline runs both tools to maximize taxonomic resolution. `geNomad` serves as the baseline layer, providing highly accurate assignments up to the Family level. To push past this boundary, the pipeline incorporates `PhaGCN`, a graph convolutional network model capable of annotating sequences down to the Species level, though its accuracy at this depth is still undergoing rigorous benchmarking with high-quality data.

To maintain an optimized balance between computational efficiency and taxonomic yield, alternative methodologies were thoroughly tested and intentionally excluded from the final architecture:

* **NCBI-database alignment approach:** While this method yielded highly similar taxonomic profiles to PhaGCN, it relies on a massive, resource-heavy reference database. Because PhaGCN utilizes a pre-trained model trained on the same data, it delivers equivalent results in a significantly smaller footprint and with far lower computational overhead.
* **VIRify ViPhOG HMM-approach:** Although its classifications aligned closely with geNomad, this method managed to annotate only 10% of the sequences that geNomad successfully resolved. Combined with its high memory demands and intense computational footprint, its marginal utility did not justify the resource expense.

## Quick Run

::::{tab-set}
:::{tab-item} Sample List

**`Conda`**

```bash
snakemake --sdm conda --use-conda --config module="viral-taxonomy" samplelist="sample/sample_list.csv" datadir="sample/fastq" outdir="sample/results" -j 16 --latency-wait 20
```

**`Apptainer`**

```bash
snakemake --sdm conda apptainer --config module="viral-taxonomy" samplelist="sample/sample_list.csv" datadir="sample/fastq" outdir="sample/results" -j 16 --latency-wait 20
```

**`HPC (PBS) (Conda)`**

```bash
EMAIL="youremail@protonmail.com"
snakemake --sdm conda --use-conda --config module="viral-taxonomy" samplelist="sample/sample_list.csv" datadir="sample/fastq" outdir="sample/results" -j 88 --latency-wait 20 --executor cluster-generic --cluster-generic-submit-cmd "qsub -N {log} -l nodes=1:ppn={threads} -l mem={resources.mem_mb}m -l walltime=120:00:00 -M $EMAIL -q cgsd -o qsub.log -e qsub.log -m a"
```

**`HPC (PBS) (Apptainer)`**

```bash
EMAIL="youremail@protonmail.com"
snakemake --sdm conda apptainer --config module="viral-taxonomy" samplelist="sample/sample_list.csv" datadir="sample/fastq" outdir="sample/results" -j 88 --latency-wait 20 --executor cluster-generic --cluster-generic-submit-cmd "qsub -N {log} -l nodes=1:ppn={threads} -l mem={resources.mem_mb}m -l walltime=120:00:00 -M $EMAIL -q cgsd -o qsub.log -e qsub.log -m a"
```

:::
:::{tab-item} Single Fasta

**`Conda`**

```bash
snakemake --sdm conda --use-conda --config module="viral-taxonomy" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results" -j 16 --latency-wait 20
```

**`Apptainer`**

```bash
snakemake --sdm conda apptainer --config module="viral-taxonomy" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results" -j 16 --latency-wait 20
```

**`HPC (PBS) (Conda)`**

```bash
EMAIL="youremail@protonmail.com"
snakemake --sdm conda --use-conda --config module="viral-taxonomy" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results" -j 88 --latency-wait 20 --executor cluster-generic --cluster-generic-submit-cmd "qsub -N {log} -l nodes=1:ppn={threads} -l mem={resources.mem_mb}m -l walltime=120:00:00 -M $EMAIL -q cgsd -o qsub.log -e qsub.log -m a"
```

**`HPC (PBS) (Apptainer)`**

```bash
EMAIL="youremail@protonmail.com"
snakemake --sdm conda apptainer --config module="viral-taxonomy" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results" -j 88 --latency-wait 20 --executor cluster-generic --cluster-generic-submit-cmd "qsub -N {log} -l nodes=1:ppn={threads} -l mem={resources.mem_mb}m -l walltime=120:00:00 -M $EMAIL -q cgsd -o qsub.log -e qsub.log -m a"
```

:::
::::

``` {admonition} Dry Run
:class: Tip
Use the `-n` flag with your command to perform a dry-run first, listing all analyses that will be performed 
```

## Configuration

* **`PhaBox2-db`:** The directory path where the PhaBox2 database is installed or will be downloaded. Defaults to the Snakemake base directory under workflow/databases. (default: "database/phabox_db_v2")
* **`phabox2-db-name`:** The designated database name or identifier file package required for the PhaBox2 classification tool execution. (default: "phabox_db_v2")
* **`phabox2-db-baselink`:** The primary remote server base link URL used to fetch and download resource updates for the PhaBox2 database. (default: "[https://github.com/KennthShang/PhaBOX/releases/download/v2](https://github.com/KennthShang/PhaBOX/releases/download/v2)")
* **`phagcn-min-len`:** The minimum allowed contig length for evaluation using the PhaGCN taxonomic assignment engine. (default: 1000)
* **`phagcn-params`:** Additional operational arguments passed to the PhaGCN classification instance. (default: "")
* **`genomad-db`:** Path to geNomad database directory. (default: "workflow/database/genomad")
* **`genomad-params-tax`:** Additional operational configurations passed to geNomad during viral taxonomic assignment. (default: "--enable-score-calibration --relaxed")

## Outputs

```bash
taxonomy/viral
├── output
    ├── viral_taxonomy_merged.tsv
    ├── geNomad
    └── PhaGCN
├── logs
├── benchmarks
└── tmp
```

The outputs from both `geNomad` and `PhaGCN` are automatically cross-referenced and consolidated into a unified taxonomic profile for each vOTU, exported as a clean, ready-to-analyze `.csv` matrix.

## {octicon}`book;0.85em` Troubleshooting Guide

We have specific guidelines for troubleshooting vomix-snakemake so we can help you out in your analysis journey as efficiently as possible! If you run into any unexpected errors, warnings, etc. please visit our [Troubleshooting Guide](/troubleshoot.md).

## {octicon}`bug;0.85em` Report a bug to us

Have any questions or you've found a bug during your analysis? Please don't hesitate to report it to us by making an issue on our [{octicon}`mark-github;0.95em` GitHub repository](https://github.com/holab-hku/vomix-snakemake/issues/new).
