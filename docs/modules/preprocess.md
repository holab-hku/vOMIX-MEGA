# `preprocess`

```{image} ../_static/figures/preprocess.svg
:width: 550
:class: no-scaled-link
:align: center
```

insert quick description of module 


## Quick Run

::::{tab-set}
:::{tab-item} Sample List
```bash
# Conda Run
snakemake --config module="preprocess" decontam-host=False outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20

# Docker Run

# Apptainer
```
::: 
::::

``` {admonition} Dry Run
:class: Tip
Use the `-n` flag with your command to perform a dry-run first, listing all analyses that will be performed 
```
``` {admonition} Dry Run
:class: Tip
Use the `-n` flag with your command to perform a dry-run first, listing all analyses that will be performed
```

## Configuration

- **`dwnldparams`:** Parameters for fasterq-dump for downloading from sra tools https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump || default: ""
- **`pigzparams`:** Parameters of pigz for compressing downloaded fastq files https://github.com/madler/pigz || default: ""
- **`fastpparams`:** Parameters to pass on fastp software https://github.com/OpenGene/fastp || default: ""
- **`hostileparams`:** Parameters for hostile decontamination https://github.com/bede/hostile|| default: ""
- **`hostilealigner`:** Which mapper to use for host decontamination- bowtie2 or minimap2 (recommended) || default: "minimap2"
- **`alignerparams`:** PLEASE DO NOT change the -x sr for minimap2 to make sure it can accurately map short reads || default: "-x sr"
- **`indexpath`:** Path to host de-contamination index for hostile || default: "./workflow/database/hostile/human-t2t-hla.fa.gz"

## Outputs

```bash
preprocess
├── samples
    ├── SampleA
        ├── output
            ├── SampleA_R1_cut.trim.filt.fastq.gz
            ├── SampleA_R2_cut.trim.filt.fastq.gz
            ├── SampleA_R1_cut.trim.filt.nodecontam.fastq.gz # If host-decontamination is performed
            └── SampleA_R2_cut.trim.filt.nodecontam.fastq.gz # If host-decontamination is performed
        ├── report.fastp.html
        ├── report.fastp.json
        ├── SampleA_R1.fastq.gz # Symbolic link
        └── SampleA_R2.fastq.gz # Symbolic link
    ├── SampleB
    ├── SampleC
├── logs
├── benchmarks
├── reports
    ├── preprocess_report.html
    └── library_size_stats.csv
└── tmp
```

The main output of this module is cleaned paired-end fastq files `SampleA_R1.fastq.gz` and `SampleA_R2.fastq.gz`. You can also find an HTML report on the quality filtering at `reports/preprocess_report.html`. Before continuing the analysis, you may filter out certain samples from your `sample_list.csv` file.

## {octicon}`book;0.85em` Troubleshooting Guide

We have specific guidelines for troubleshooting vOMIX-snakemake so we can help you out in your analysis journey as efficiently as possible! If you run into any unexpected errors, warnings, etc. please visit our [Troubleshooting Guide](/troubleshoot.md).

### {octicon}`bug;0.85em` Report a bug to us !

Have any questions or you've found a bug during your analysis? Please don't hesitate to report it to us by making an issue on our [{octicon}`mark-github;0.95em` GitHub repository](https://github.com/holab-hku/vOMIX-MEGA/issues/new).
