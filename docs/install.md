# Installation

## Conda, Mamba, & Pip

You can install **vOMIX‑MEGA** using general‑purpose package managers like **Conda** or **Mamba**, a Python‑only dependency manager like **pip**, or a fully containerised approach via **Apptainer**. We provide a **Bioconda** package for stable releases, while the latest development version can be installed directly from GitHub.

::::{tab-set}

:::{tab-item} Conda
[Conda](https://docs.conda.io/projects/conda/en/stable/) is a package manager that handles all your dependencies for you. To install vOMIX-MEGA using Conda, you can create the environment from the repository environment file.

```bash
# Install via Bioconda
conda create -n vomix -c bioconda vomix-mega
conda activate vomix

# Verify Installation
vomix -h
```

:::

:::{tab-item} Mamba
[Mamba](https://mamba.readthedocs.io/en/latest/index.html) is a re‑implementation of Conda that uses a faster, more robust dependency solver. The commands are almost identical:

```bash
# Install via Bioconda
mamba create -n vomix -c bioconda vomix-mega
mamba activate vomix

# Verify Installation
vomix -h
```

:::

:::{tab-item} Conda Lock (Source)
[Conda‑Lock](https://github.com/conda/conda-lock) installs from a lock file with fully pinned dependencies. This is often the most reliable fallback if the standard conda/mamba installation fails.

```bash
# Download GitHub directory
git clone https://github.com/holab-hku/vOMIX-MEGA.git
cd vOMIX-MEGA

# Create a conda-lock environment 
# if you don't want to install into base environment
conda create -n conda-lock -c conda-forge conda-lock=4.0.2 -y
conda activate conda-lock

# Use conda-lock to install the environment from the repository lock file
conda-lock install --name vomix conda-lock.yml
conda deactivate # deactive conda-lock environment
conda activate vomix # activate vomix environment

# Install using pip
pip install .

# Verify Installation
vomix -h
```

:::

:::{tab-item} pip (Source)
[Pip](https://pypi.org/project/pip/) is a package manager written in Python and is used to install and manage software packages. The Python Software Foundation recommends using pip to install Python applications and its dependencies during deployment.

```{admonition} Pip Installation
:class: attention
Although  base dependcies of `vOMIX-MEGA` are all pip-installable, underlying module tools rely on non-python dependencies such as `megahit`, `geNomad`, and more. Running specific modules ultimately requires either `conda` or `apptainer` installed. 
```

```bash
# Download GitHub directory
git clone https://github.com/holab-hku/vOMIX-MEGA.git
cd vOMIX-MEGA

# Install base environment
conda create -n vomix

# Activate environment
conda activate vomix

# Install using pip
pip install .

# Verify Installation
vomix -h
```

:::

::::

```{admonition} Conda update
:class: attention
If you are using conda or mamba, make sure your conda installation is up to date before proceeding. You can update it with `conda update -n base -c defaults conda`
```

```{admonition} Conda Channel Priorities
:class: attention
If you are using conda or mamba, make sure to set channel orders correctly and set channel priority to strict. Via the `conda config --add channels defaults`, `conda config --add channels bioconda`, `conda config --add channels conda-forge`, and `conda config --set channel_priority strict` respectively. For mamba replace `conda` with `mamba` respectively.  
```

```{admonition} Conda Lock as a fallback
:class: note
If the standard conda or mamba installation methods do not work, `conda-lock` is usually the best approach because it installs the pinned dependency set from the repository lock file. Snakemake and many bioinformatics tools rely heavily on POSIX-compliant workflows, so they are currently best supported on non-Windows operating systems.
```

## Apptainer Installation

vOMIX-MEGA is built on a snakemake back-end, which facilitates native containerized deployment via an `Apptainer` (formerly `Singularity`) `.sif` image. The container image generated contains explicitly each conda environment mounted on top of a base operating system. Containers are preferred for the most robust forms of reproducibility, whereas `conda` and `mamba` installations might not work on Windows or Mac-ARM systems.

```{admonition} Install Apptainer
:class: tip
Make sure you have `Apptainer` installed before you run the following commands. The advantage of Apptainer over docker is that it allows non-privilieged installation without root permissions. View the full documentation at the [Apptainer Wiki](https://apptainer.org/docs/admin/main/installation.html). Note that while the prebuilt image does not require root permission, a manual build might require it.
```

::::{tab-set}
:::{tab-item} Apptainer (Download Prebuilt Image)

```bash
# Enter vOMIX-MEGA directory
# replace this with your native installation path
cd vOMIX-MEGA

# Pull Container Image 
VOMIX_VERSION="v0.1.0-beta.1"
apptainer pull --name workflow/apptainer/vomix_${VOMIX_VERSION}.sif oras://ghcr.io/erfanshekarriz/vomix:${VOMIX_VERSION}

# Dry run (test installation)
vomix viral-identify --sdm apptainer --fasta sample/contigs/contigs_simulated_viral_nonviral.fasta --outdir quick-run/results -j 64 --latency-wait 20 -n
```

:::
:::{tab-item} Apptainer (Build Local Image)

```bash
# Enter vOMIX-MEGA directory
# replace this with your native installation path
cd vOMIX-MEGA

# Built container Image 
VOMIX_VERSION="v0.1.0-beta.1"
apptainer build workflow/apptainer/vomix_${VERSION}.sif workflow/apptainer/vomix_${VERSION}.def

# Dry run (test installation)
vomix viral-identify --sdm apptainer --fasta sample/contigs/contigs_simulated_viral_nonviral.fasta --outdir quick-run/results -j 64 --latency-wait 20 -n
```

:::

::::

```{admonition} Apptainer Image Size
:class: note
The full vOMIX-MEGA apptainer build takes `9.5 GB` of storage, which includes all pre-installed conda packages mounted on a `condaforge/miniforge3:latest` base operating system (not including databases). Make sure that you have this space availabe on your local machine before running. Due to the filesize, your HPC administrator might limit the download bandwith, so you might need to contact them to set it up for you. Running `export GODEBUG="http2client=0"` before your pull might be a quick fix for this. 
```

```{admonition} Using Conda within a Container
:class: note
Each rule in vOMIX-MEGA's underlying snakemake files depends on a speicifc conda environment. To run what Snakemake calls ["Ad-hoc combination of Conda package management with containers"](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#ad-hoc-combination-of-conda-package-management-with-containers), which is essentially running apptainer containers with conda enviornments installed within them, you need to use the `--sdm conda --sdm apptainer` options in that specific format and order. This allows true full reproducibility.
```

```{admonition} Wrong Apptainer Flags for Snakemake
:class: warning
If you only use `--sdm apptainer`, Snakemake will not launch any conda environments and hence all jobs will fail. If you use `--sdm apptainer --use-conda` it will try and re-install conda enviornments in your local `.snakemake/conda` folder, which counteracts the purpose of containers. Using `--sdm apptainer --sdm conda` will also not work successfully.
```

## {octicon}`book;0.85em` Troubleshooting Guide

We have specific guidelines for troubleshooting vOMIX-MEGA so we can help you out in your analysis journey as efficiently as possible! If you run into any unexpected errors, warnings, etc. please visit our [Troubleshooting Guide](/troubleshoot.md).

## {octicon}`bug;0.85em` Report a bug to us

Have any questions or you've found a bug during your analysis? Please don't hesitate to report it to us by making an issue on our [{octicon}`mark-github;0.95em` GitHub repository](https://github.com/holab-hku/vOMIX-MEGA/issues/new).
