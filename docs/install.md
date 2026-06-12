# Installation

## Installing vOMIX-snakemake

You can install geNomad in you computer using either a general-purpose package manager (`pixi`, `mamba`, `conda`) or a Python-specific package manager (`uv` or `pip`).

::::{tab-set}

:::{tab-item} Pixi
[Pixi](https://pixi.sh/latest/) is probably the simplest way to install geNomad. It takes care of all the dependencies for you and doesn't require any environment management, meaning the `genomad` command will be available globally.

```bash
pixi global install -c conda-forge -c bioconda genomad
genomad --help
```
:::

:::{tab-item} Mamba
[Mamba](https://mamba.readthedocs.io/en/latest/) is a package manager that handles all your dependencies for you. To install vOMIX-snakemake mad using Mamba, you need to create a new environment and activate it before running the `snakemake` command.

```bash
# Make sure conda is updated for snakemake compatibility [IMPORTANT]
# Set channel priority to strict before running vOMIX-snakemake to ensure reproducibility [IMPORTANT]
mamba config --add channels bioconda
mamba config --add channels conda-forge
mamba config --set channel_priority strict

# Install base environment
mamba create -n vomix -c conda-forge snakemake=8.25.5 biopython=1.84 -y
mamba activate vomix

# Verify the two essential base tools are running
snakemake -v
```
:::

:::{tab-item} Conda
[Conda](https://docs.conda.io/projects/conda/en/stable/) is a package manager that handles all your dependencies for you. To install vOMIX-snakemake using Conda, you need to create a new environment and activate it before running the `snakemake` command.

```bash
# Make sure conda is updated for snakemake compatibility [IMPORTANT]
# Set channel priority to strict before running vOMIX-snakemake to ensure reproducibility [IMPORTANT]
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# Install base environment
conda create -n vomix -c conda-forge snakemake=8.25.5 biopython=1.84 -y
conda activate vomix

# Verify the two essential base tools are running
snakemake -v
```
:::

::::

```{admonition} Non-python dependencies
:class: attention
Pixi, Mamba, and Conda will install both the Python dependencies and the non-Python dependencies required by geNomad. If you install geNomad using `uv` or `pip`, make sure to add [`MMseqs2`](https://github.com/soedinglab/MMseqs2/) and [`ARAGORN`](http://www.ansikte.se/ARAGORN/) to your `$PATH`.
```

## Running geNomad using containers

You can also execute geNomad using containerization tools, such as Docker and Podman. To pull the image, execute the command below.

::::{tab-set}

:::{tab-item} Docker
```bash
docker pull antoniopcamargo/genomad
```
:::

:::{tab-item} Podman
```bash
podman pull docker.io/antoniopcamargo/genomad
```
:::

::::

To start a geNomad container you have to mount a folder from the host system into the container with the `-v` argument. The following command mounts the current working directory (`$(pwd)`) under `/app` inside the container and then executes the `genomad download-database` and `genomad end-to-end` commands.

::::{tab-set}

:::{tab-item} Docker
```bash
docker run -ti --rm -v "$(pwd):/app" antoniopcamargo/genomad download-database .
docker run -ti --rm -v "$(pwd):/app" antoniopcamargo/genomad end-to-end input.fna output genomad_db
```
:::

:::{tab-item} Podman
```bash
podman run -u 0 -ti --rm -v "$(pwd):/app" antoniopcamargo/genomad download-database .
podman run -u 0 -ti --rm -v "$(pwd):/app" antoniopcamargo/genomad end-to-end input.fna output genomad_db
```
:::

::::

## ⚙️ Troubleshooting Guide

We have specific guidelines for troubleshooting vOMIX-snakemake so we can help you out in your analysis journey as efficiently as possible! If you run into any unexpected errors, warnings, etc. please visit our Troubleshooting Guide.

