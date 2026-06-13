# Quick Start

vOMIX-snakemake is a comprehensive pipeline, and we know that big-data analysis can feel intimidating. Still, we'd like to show you how easy it can be! vOMIX-MEGA is here to simplify everything, remove tedious ad-hoc coding, and speed up viral metagenomic analysis by miles, all while maintaining full freedom to tune parameters to your liking. Here are three different quick analysese you can do. 

## Tutorial 1: Identify Viruses in Contig File

You can quickly analyse and analyse a mock dataset of viral and non-viral mixed contigs using the `viral-identify` module.The sample data should already be included in your directory, or can be downloaded via `wget https://github.com/holab-hku/vomix-snakemake/tree/main/sample`.

``` {admonition} Databases & Environments
:class: Note
Conda environments for each module will be installed automatically at the beginning of the command. Databases will be downloaded for each software automatically before the command is run.
```

``` {admonition} Skip Database Setup
:class: Note
If you would like to skip database setup, you can add `--config setup-database=False` to your command.
```

::::{tab-set}
:::{tab-item} Conda
```bash
# Use more memory (22 GB) but run faster (~10 mins)
snakemake --use-conda --config module="viral-identify" outdir="test_res" splits=0 fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" -j 4 --latency-wait 20

# Use less memory (8 GB) but run slower (~40 mins)
snakemake --use-conda --config module="viral-identify" outdir="test_res" splits=8 fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" -j 4 --latency-wait 20
```
:::
:::{tab-item} Docker
```bash
snakemake --use-container --sdm conda --config module="viral-identify" outdir="test_res" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" splits=0 -j 4 --latency-wait 20
```
:::
:::{tab-item} Apptainer
```bash
snakemake --sdm conda apptainer  --config module="viral-identify" outdir="test_res" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" splits=0 -j 4 --latency-wait 20
```
:::
::::


```{admonition} Inputs
:class: note 
the `viral-identify` module can either take a single fasta file as input via the `--config fasta="sample.fasta"` input, or a directory of different fasta files corresponding to contig files generated from multiple samples via the `--config fastadir="sample/fastadir/" input.
``` 

::::{tab-set}
:::{tab-item} Single Fasta File Input
```bash
snakemake --config module="viral-identify" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="test_res"  splits=0 -j 4 --latency-wait 20
```
:::
:::{tab-item} Fasta Directory Input
```bash 
snakemake --config module="viral-identify" fastadir="sample/contigs" outdir="test_res" splits=0 -j 4 --latency-wait 20
```
:::
::::

```{admonition} Memory Footprint
:class: note
The standard viral identification analysis of vOMIX-snakemake is desinged to be at a maximum of 24Gb, used by geNomad during viral contig identification. If you are on a small computer, you can further reduce this number by introducing the `--config splits=8` which will reduce memory use at the expense of computation time.
```

## Tutorial 2: Benchmark Viral Contig Tools


## Tutorial 3: Perform End-to-End Analysis on SRA Samples

You can to full end-to-end viral analysis on a set of samples using the `

You can quickly analyse and analyse a mock dataset of viral and non-viral mixed contigs using the `viral-identify` module.The sample data should already be included in your directory, or can be downloaded via `wget https://github.com/holab-hku/vomix-snakemake/tree/main/sample`.

::::{tab-set}
:::{tab-item} Conda
```bash
# Dry run (check jobs)
snakemake --use-conda --configfile config/config.yml --config module="end-to-end" decontam-host=False outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" -j 4 --dry-run

# Run your pipeline
snakemake --use-conda --configfile config/config.yml --config module="end-to-end" decontam-host=False outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" -j 4 --latency-wait 20
```
:::
:::{tab-item} Docker
```bash
snakemake --use-container --sdm condaa --configfile config/config.yml --config module="end-to-end" decontam-host=False outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" -j 4 --latency-wait 20
```
:::
:::{tab-item} Apptainer
```bash
snakemake --sdm conda apptainer --configfile config/config.yml --config module="end-to-end" decontam-host=False outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" -j 4 --latency-wait 20
```
:::
::::

```{admonition} Sample List Format
:class: Tip 
The `sample_list.csv` takes SRA accessions and automatically downloads and processes the data for you. Please see the [guide](/advanced.md#sample-list-csv) to setting up your sample list. 
```



> **NOTE:** To see all snakemake parameters run `snakemake -h`. Configurations nested after `--config` are all specific to vOMIX-MEGA and can only be accessed via the `config/config.yml` file. Please read more in our Run & Configuration documentation page.

The pipeline handles all processes. To run on a cluster, you can visit Advanced Usage for a detailed setup guide. To change the pipeline configuration, you may change the `config.yml` file or visit Run & Configuration for a detailed run and configuration explanation.

> **Note:** Conda environments for each module will be installed automatically at the beginning of the command. Databases will be downloaded for each software automatically before the command is run.

## ⚙️ Troubleshooting Guide

We have specific guidelines for troubleshooting vOMIX-MEGA so we can help you out in your analysis journey as efficiently as possible! If you run into any unexpected errors, warnings, etc. please visit our Troubleshooting Guide.

