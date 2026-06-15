# Quick Start

vOMIX-snakemake is a comprehensive pipeline, and we know that big-data analysis can feel intimidating. We'd like to make that as easy as possible for you. vOMIX-snakemake simplifies analysis by removing tedious ad-hoc coding, and speeds up viral metagenomic analysis through job optimization and mecnhmarking, all while maintaining full freedom to tune parameters to your liking. Here are three different quick analysese you can do. 

## {octicon}`checkbox;0.85em` Quick Tutorial #1: Identify Viruses in Mock Contig File

You can quickly analyse a mock dataset of 988 viral and non-viral mixed contigs using the `viral-identify` module.The sample data should already be included in your directory, or can be downloaded via `wget https://github.com/holab-hku/vomix-snakemake/tree/main/sample`. 

::::{tab-set}
:::{tab-item} Conda
```bash
# Activate conda environment
conda activate vomix

# Use more memory (22 GB) but run faster (~10 mins)
snakemake --use-conda --config module="viral-identify" outdir="test_res" splits=0 fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" -j 64 --latency-wait 20

# Use less memory (8 GB) but run slower (~40 mins)
snakemake --use-conda --config module="viral-identify" outdir="test_res" splits=8 fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" -j 64 --latency-wait 20

# Quick look at results
head -n 20 test_res/identify/viral/output/classification_summary_vOTUs.csv
```
:::
:::{tab-item} Docker
```bash
# Activate conda environment
conda activate vomix

# Use more memory (22 GB) but run faster (~10 mins)
snakemake --use-container --sdm conda --config module="viral-identify" outdir="test_res" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" splits=0 -j 64 --latency-wait 20

# Use less memory (8 GB) but run slower (~40 mins)
snakemake --use-container --sdm conda --config module="viral-identify" outdir="test_res" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" splits=8 -j 64 --latency-wait 20

# Quick look at results
head -n 20 test_res/identify/viral/output/classification_summary_vOTUs.csv 
```
:::
:::{tab-item} Apptainer
```bash
# Activate conda environment
conda activate vomix

# Use more memory (22 GB) but run faster (~10 mins)
snakemake --sdm conda apptainer --config module="viral-identify" outdir="test_res" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" splits=0 -j 64 --latency-wait 20

# Use less memory (8 GB) but run slower (~40 mins)
snakemake --sdm conda apptainer --config module="viral-identify" outdir="test_res" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" splits=8 -j 64 --latency-wait 20

# Quick look at results
head -n 20 test_res/identify/viral/output/classification_summary_vOTUs.csv 
```
:::
::::

``` {admonition} Databases & Environments
:class: Note
Conda environments for each module will be installed automatically at the beginning of the command. Databases will be downloaded for each software automatically before the command is run.
```

```{admonition} Memory Footprint
:class: note
The standard viral identification analysis of vOMIX-snakemake is desinged to be at a maximum of 24Gb, used by geNomad during viral contig identification. If you are on a small computer, you can further reduce this number by introducing the `--config splits=8` which will reduce memory use at the expense of computation time.
```

## {octicon}`checkbox;0.85em` Quick Tutorial #2: Benchmark 6 Different Viral Contig Identification Tools

Although not part of the standard vOMIX-snakemake analysis, the `viral-benchmark` module for you to quickly benchmark 6 different viral contig identification tools including geNomad, DeepVirFinder, Phamer, VirSorter2, VirFinder, and VIBRANT. You can view all softwares cited via our [references](/reference.md) page. Please note that we have extensively [benchmarked](/benchmark.md) geNomad to be the best performing general tool for viral contig identification.

::::{tab-set}
:::{tab-item} Conda
```bash
# Activate conda environment
conda activate vomix

# Use more memory (22 GB) but run faster (~10 mins)
snakemake --use-conda --config module="viral-benchmark" outdir="test_res" splits=0 fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" -j 4 --latency-wait 20

# Use less memory (8 GB) but run slower (~40 mins)
snakemake --use-conda --config module="viral-benchmark" outdir="test_res" splits=8 fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" -j 4 --latency-wait 20

# Quick look at results
head -n 20 test_res/identify/viral/output/benchmark_results_merged.csv
```
:::
:::{tab-item} Docker
```bash
# Activate conda environment
conda activate vomix

# Use more memory (22 GB) but run faster (~10 mins)
snakemake --use-container --sdm conda --config module="viral-benchmark" outdir="test_res" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" splits=0 -j 4 --latency-wait 20

# Use less memory (8 GB) but run slower (~40 mins)
snakemake --use-container --sdm conda --config module="viral-benchmark" outdir="test_res" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" splits=8 -j 4 --latency-wait 20

# Quick look at results
head -n 20 test_res/identify/viral/output/benchmark_results_merged.csv
```
:::
:::{tab-item} Apptainer
```bash
# Activate conda environment
conda activate vomix

# Use more memory (22 GB) but run faster (~10 mins)
snakemake --sdm conda apptainer --config module="viral-benchmark" outdir="test_res" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" splits=0 -j 4 --latency-wait 20

# Use less memory (8 GB) but run slower (~40 mins)
snakemake --sdm conda apptainer --config module="viral-benchmark" outdir="test_res" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" splits=8 -j 4 --latency-wait 20

# Quick look at results
head -n 20 test_res/identify/viral/output/benchmark_results_merged.csv
```
:::
::::

## {octicon}`checkbox;.85em` Quick Tutorial #003: Perform End-to-End Analysis on SRA Samples

The main purpose of vOMIX-snakemake is to make end-to-end viral metagenomic analysis straightforward. You can provide a list of SRA accessions or combine them with local fastq.gz  files via the `sample_list.csv`. The `viral-end-to-end` module will automatically run `preprocess`, `assembly`, `viral-identify`, `viral-taxonomy`, `viral-host`, and `viral-annotate`.

::::{tab-set}
:::{tab-item} Conda
```bash
# Activate conda environment
conda activate vomix

# Dry run (check jobs)
snakemake --use-conda --configfile config/config.yml --config module="end-to-end" decontam-host=False outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" -j 64 --dry-run

# Run your pipeline
snakemake --use-conda --configfile config/config.yml --config module="end-to-end" decontam-host=False outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" -j 64 --latency-wait 20

# Quick look at results
head -n 20 test_res/identify/viral/output/classification_summary_vOTUs.csv 
```
:::
:::{tab-item} Docker
```bash
# Activate conda environment
conda activate vomix

snakemake --use-container --sdm condaa --configfile config/config.yml --config module="end-to-end" decontam-host=False outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" -j 64 --latency-wait 20

# Quick look at results
head -n 20 test_res/identify/viral/output/classification_summary_vOTUs.csv 
```
:::
:::{tab-item} Apptainer
```bash
# Activate conda environment
conda activate vomix

snakemake --sdm conda apptainer --configfile config/config.yml --config module="end-to-end" decontam-host=False outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" -j 64 --latency-wait 20

# Quick look at results
head -n 20 test_res/identify/viral/output/classification_summary_vOTUs.csv 
```
:::
::::

```{admonition} Sample List
:class: Tip 
The `sample_list.csv` takes SRA accessions and automatically downloads and processes the data for you. Please see the [sample list guide](/advanced.md#sample-list-csv) to setting up your own custom data for analysis. 
```
```{admonition} Custom Configuration
:class: Tip
You can configure all underlying commands and configuration of vOMIX-snakemake either via the `--configfile config.yml` file, or by providing specific parameters via the `--config {name}={value} snakemake directive. Please see the full [configuration guide](/advanced.md#configuration) to customize your run.`
```
```{admonition} Further Analysis
:class: Tip
vOMIX-snakemake can perform more than just what the `viral-end-to-end` module can do. To view the full module list, visit the [Run & Configuration Guide](/run.md).
```

## {octicon}`book;0.85em` Troubleshooting Guide

We have specific guidelines for troubleshooting vOMIX-snakemake so we can help you out in your analysis journey as efficiently as possible! If you run into any unexpected errors, warnings, etc. please visit our [Troubleshooting Guide](/troubleshoot.md).

### {octicon}`bug;0.85em` Report a bug to us !

Have any questions or you've found a bug during your analysis? Please don't hesitate to report it to us by making an issue on our [{octicon}`mark-github;0.95em` GitHub repository](https://github.com/holab-hku/vOMIX-MEGA/issues/new).
