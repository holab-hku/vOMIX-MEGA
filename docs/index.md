# vOMIX-snakemake

vOMIX-snakemake is the back-end pipline of the command-line tool vOMIX-MEGA: a reproducible, scalable, and fast viral metagenomic pipeline for analyzing large-scale bulk-metagenomic and viromic data. We have engineered multiple bottlenecks in current state-ofthe-art software to allow rapid and well-benchmarked viral metagenomic analysis. vOMIX-snakemake and vOMIX-MEGA (referred collectively as vOMIX) provides numerous modules for end-to-end. Here are some of it's handy features !

::::{grid} :gutter: 2

:::{grid-item-card} :columns: 12 12 4 4 High Speed ^^^ vOMIX operates 10-1000 times faster than current pipelines due to engineering of multiple slow underlying software:::

:::{grid-item-card} :columns: 12 12 4 4 Stable Memory Footprint ^^^ The standard viral end-to-end analysis takes a maximum memory of 24 Gb due to selection and fine-tuning of tools:::

:::{grid-item-card} :columns: 12 12 4 4 Modular Personalizable Analysis  ^^^  Each module can be separately used with a variety of different inputs, and each software can be fine-tuned for your needs :::


::::

::::{grid} :gutter: 2

:::{grid-item-card} :columns: 12 12 4 4 Benchmarked  ^^^ We've benchmarked our viral identification on experimental as well as mock-data and use the best performing tools  :::

:::{grid-item-card} :columns: 12 12 4 4 Reproducible & Dockerized ^^^ All analysis is done via a snakemake back-end which is configured with Docker to allow full reproducibility:::

:::{grid-item-card} :columns: 12 12 4 4 Easy SRA Input  ^^^ Simply feed vOMIX your SRA accession codes and it will download, process, and analyse the viral community of your samples automatically. :::

::::


# {rocket}`rocket;0.85em` Getting Started

:::{card} Installation
:link: installation
:link-type: doc
Instructions on how to install vOMIX-snakemake on your computer or server.
:::

:::{card} Quickstart
:link: quickstart
:link-type: doc
Learn how to run vOMIX-snakemake on a sample dataset.
:::


# {octicon}`bookmark;0.85em` Citing vOMIX-MEGA

If you use vOMIX-MEGA in your work, please consider citing its pre-print manuscript:

:::{card}
:link: insert link heret

**vOMIX-MEGA: A critical speed enhancement for end-to-end viral metagenomics**
Erfan Shekarriz, Elsa VIJENDRAN, Joshua WK Ho  — *bioRxiv* (2026), DOI: XXXXXXXXXXXXXXXXXXX.
:::

## {octicon}`bug;0.85em` Report a bug to us ! 

Have any questions or you've found a bug during your analysis? Please don't hesitate to report it to us by making an issue on our [{octicon}`mark-github;0.95em` GitHub repository](https://github.com/holab-hku/vOMIX-MEGA).




```{toctree}
:hidden:
:caption: Introduction

self
```


```{toctree}
:maxdepth: 3
:caption: vomix-snakemake

install
quickstart
run
outputs
troubleshoot
```

```{toctree}
:maxdepth: 3 
:caption: Advanced Usage

advanced
```
