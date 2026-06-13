# Installation

## Conda & Mamba

You can install vOMIX-snakemake in you computer using general-purpose package managers ( `mamba`, `conda`). .

::::{tab-set}

:::{tab-item} Mamba
[Mamba](https://mamba.readthedocs.io/en/latest/) is a package manager that handles all your dependencies for you. To install vOMIX-snakemake mad using Mamba, you need to create a new environment and activate it before running the `snakemake` command.

```bash
# Install base environment
mamba create -n vomix -c conda-forge snakemake=8.25.5 biopython=1.84 -y

# Activate environment
mamba activate vomix

# Verify Installation
snakemake -v
```
:::

:::{tab-item} Conda
[Conda](https://docs.conda.io/projects/conda/en/stable/) is a package manager that handles all your dependencies for you. To install vOMIX-snakemake using Conda, you need to create a new environment and activate it before running the `snakemake` command.

```bash
# Install base environment
conda create -n vomix -c conda-forge snakemake=8.25.5 biopython=1.84 -y

# Activate environment
conda activate vomix

# Verify Installation
snakemake -v
```
:::

::::

```{admonition} Conda Channel Priorities
:class: attention
If you are using conda or mamba, make sure to set channel orders correctly and set channel priority to strict. Via the `conda config --add channels bioconda`, `conda config --add channels conda-forge`, and `conda config --set channel_priority strict` respectively. For mamba replace `conda` with `mamba` respectively.  
```

## Docker & Apptainer

vOMIX-snakemake is built on a snakemake back-end, which facilitates native `Docker` and `Apptainer` employment. The container image generated a container image contains all required environments. This makes our pipeline truly and fully-reproducible, whereas `conda` and `mamba` installations might not work on Windows or Mac-ARM systems. 

::::{tab-set}

:::{tab-item} Docker
```bash
snakemake --use-container --use-conda --cores 4
```
:::

:::{tab-item} Apptainer
```bash
snakemake --software-deployment-method conda apptainer
```
:::

::::

```{admonition} Using Conda within Containers
:class: attention
All jobs within vOMIX-snakemake have specified conda environments. The containers built via snakemake is then just a light-weight image of a standard linux OS with conda environments being built on top, meaning it still relies on conda environment installation. You must set `--software-deployment-method conda apptainer` or `--sdm conda apptainer`, even if you are running on containers. 
```

## ⚙️ Troubleshooting Guide

We have specific guidelines for troubleshooting vOMIX-snakemake so we can help you out in your analysis journey as efficiently as possible! If you run into any unexpected errors, warnings, etc. please visit our [Troubleshooting Guide](/troubleshoot.md).

