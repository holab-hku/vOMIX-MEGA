[![Documentation Status](https://readthedocs.org/projects/vomix-mega/badge/?version=latest)](https://vomix-mega.readthedocs.io/en/latest/?badge=latest)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Apptainer](https://img.shields.io/badge/container-Apptainer%2FSingularity-orange.svg)](https://apptainer.org/)

**vOMIX-MEGA** is a highly scalable, fast, and reproducible viral metagenomic pipeline designed for large-scale bulk-metagenomic and viromic data analysis. Built on top of a Snakemake backend with automated Conda and Apptainer/Singularity container support, vOMIX-MEGA resolves performance bottlenecks in current state-of-the-art tools to achieve dramatic speedups without sacrificing accuracy.

📖 **Full Documentation:** [vomix-mega.readthedocs.io](https://vomix-mega.readthedocs.io/en/latest/)

---

## ⚡ Key Features

* **⚡ Unmatched Speed:** Operates **10–1000× faster** than traditional unoptimized pipelines by fine-tuning underlying software dependencies and optimizing parallel job scheduling.
* **🧠 Predictable & Low Memory Footprint:** Standard viral end-to-end analysis requires a maximum of **24 GB RAM** (further tuneable down to **8 GB** for resource-constrained systems).
* **🧩 Modular Architecture:** Run standalone sub-workflows (`preprocess`, `viral-identify`, `viral-host`, `prok-binning`, etc.) or execute full end-to-end analysis.
* **🔬 Extensively Benchmarked:** Evaluated across mock and experimental datasets to ensure top-tier sensitivity and precision.
* **📦 Reproducible & Containerized:** Full support for `Apptainer` / `Singularity` containers combined with per-rule Conda environments ensures 100% computational reproducibility.
* **📥 Seamless SRA Input:** Simply feed SRA accessions to download, preprocess, assemble, and analyze viral communities automatically.

---

## 🏗️ Pipeline Overview & Modules

vOMIX-MEGA structures complex viral and prokaryotic metagenomic workflows into dedicated, easy-to-use modules:

| Category | Module | Description |
| :--- | :--- | :--- |
| **Core Viral Workflow** | `viral-identify` | High-throughput identification of viral contigs and proviruses |
| | `viral-taxonomy` | Taxonomic classification of identified viral contigs (vOTUs) |
| | `viral-host` | Host prediction using CRISPR spacers and tRNA/sequence matching |
| | `viral-annotate` | Functional annotation of viral genes and genomes |
| | `viral-community` | Relative abundance estimation and community profiling |
| | `viral-end-to-end` | Full pipeline execution from raw reads to community profiling |
| **Benchmarking & Utils** | `viral-benchmark` | Compare 6 viral identification tools (geNomad, DeepVirFinder, Phamer, VirSorter2, VirFinder, VIBRANT) |
| | `cluster-fast` | High-speed sequence clustering for vOTU definition |
| | `checkv-pyhmmer` | Ultra-fast CheckV quality assessment powered by PyHMMER |
| | `setup-database` | Automated setup and verification of required reference databases |
| **Prokaryotic Modules** | `prok-binning` | Metagenomic binning for prokaryotic genomes (MAGs) |
| | `prok-annotate` | Functional annotation of bacterial/archaeal contigs and MAGs |
| | `prok-community` | Taxonomic and functional profiling of prokaryotic communities |
| **Upstream Processing** | `preprocess` | Quality filtering, adapter trimming, and host decontam |
| | `assembly` | De novo metagenomic assembly optimized for viral discovery |

---

## 📦 Installation

vOMIX-MEGA can be installed via **Conda / Mamba** or run using an **Apptainer / Singularity** container.

### Option 1: Conda or Mamba (Recommended)

```bash
# 1. Clone the repository
git clone [https://github.com/holab-hku/vOMIX-MEGA.git](https://github.com/holab-hku/vOMIX-MEGA.git)
cd vOMIX-MEGA

# 2. Create and activate the conda environment
mamba env create -f environment.yml   # Or: conda env create -f environment.yml
mamba activate vomix

# 3. Install the vOMIX CLI package
pip install .

# 4. Verify installation
vomix -h

```

> **Note on Conda Lock:** If dependency resolution fails on your environment, use our pinned `conda-lock` file:
>
> ```bash
> conda-lock install --name vomix conda-lock.yml
> conda activate vomix
> pip install .
> 
> ```
>
>

---

### Option 2: Apptainer / Singularity Container

For HPC clusters or environments where Conda environments are restricted:

```bash
cd vOMIX-MEGA

# Pull the official Apptainer container (~9.5 GB)
VOMIX_VERSION="v0.1.0-beta.1"
apptainer pull --name workflow/apptainer/vomix_${VOMIX_VERSION}.sif oras://ghcr.io/erfanshekarriz/vomix:${VOMIX_VERSION}

# Test installation with a dry-run
vomix viral-identify --sdm apptainer --fasta sample/contigs/contigs_simulated_viral_nonviral.fasta --outdir quick-run/results -j 64 --latency-wait 20 -n

```

---

## 🚀 Quick Start Tutorial

Identify viral contigs in a sample FASTA dataset containing mixed viral and non-viral sequences using `viral-identify`.

```bash
# Activate your environment
conda activate vomix

# 1. Test run (Dry-run to verify workflow execution steps)
vomix viral-identify \
  --sdm conda \
  --fasta sample/contigs/contigs_simulated_viral_nonviral.fasta \
  --outdir quick-run/results \
  --contig-splits 0 \
  --checkv-splits 0 \
  -j 64 \
  --latency-wait 20 \
  -n

# 2. Run high-performance mode (~22 GB RAM, fast execution ~10 mins)
vomix viral-identify \
  --sdm conda \
  --fasta sample/contigs/contigs_simulated_viral_nonviral.fasta \
  --outdir quick-run/results \
  --contig-splits 0 \
  --checkv-splits 0 \
  -j 64 \
  --latency-wait 20

# 3. Run low-memory mode (~8 GB RAM for smaller machines)
vomix viral-identify \
  --sdm conda \
  --fasta sample/contigs/contigs_simulated_viral_nonviral.fasta \
  --outdir quick-run/results \
  --contig-splits 7 \
  --checkv-splits 7 \
  -j 64 \
  --latency-wait 20

# 4. View identification summary output
head -n 20 quick-run/results/identify/viral/output/classification_summary_vOTUs.csv

```

---

## 📚 Complete Documentation

Detailed guides, advanced cluster submission options (SLURM / PBS / SGE), parameter specifications, and troubleshooting tips are available at our official documentation site:

👉 **[vOMIX-MEGA ReadTheDocs](https://vomix-mega.readthedocs.io/en/latest/)**

* 🛠️ [Installation Details & Troubleshooting](https://vomix-mega.readthedocs.io/en/latest/install.html)
* ⚡ [Quick Start Guide](https://vomix-mega.readthedocs.io/en/latest/quickstart.html)
* 📖 [Module Specifications & Parameters](https://vomix-mega.readthedocs.io/en/latest/run.html)
* 🔧 [Troubleshooting Guide](https://vomix-mega.readthedocs.io/en/latest/troubleshoot.html)

---

## 📄 Citation

If you use vOMIX-MEGA in your research, please cite our pre-print:

> **vOMIX-MEGA: A critical speed enhancement for end-to-end viral metagenomics**
> Erfan Shekarriz, Elsa Vijendran, Joshua WK Ho
> *bioRxiv* (2026). DOI: [10.1101/vomix-mega](https://vomix-mega.readthedocs.io/en/latest/)

---

## 🐛 Bug Reports & Community Support

Found a bug or have a feature request? Please open an issue on our [GitHub Issues Page](https://github.com/holab-hku/vOMIX-MEGA/issues/new).

Copyright © 2026 [Ho Lab, The University of Hong Kong (HKU)](https://github.com/holab-hku). Distributed under the MIT License.
"""

with open("README.md", "w", encoding="utf-8") as f:
f.write(readme_content)

print("README.md created successfully!")

```

```text?code_stdout&code_event_index=1
README.md created successfully!


```

Your README.md file is ready

---

### Summary of `README.md` Features & Layout

The generated **`README.md`** file includes:

1. **Header & Badges:** Features ReadTheDocs status, license (MIT), Python version, and Apptainer container badges.
2. **Overview & Key Highlights:** Summarizes vOMIX-MEGA’s key benefits:

* **10–1000× Speedup** over unoptimized legacy pipelines.
* **Predictable Memory Footprint** ($\le 24\text{ GB}$ peak, configurable down to $8\text{ GB}$).
* **Modular Architecture & Containerization** (`Conda` + `Apptainer/Singularity`).
* **Direct SRA Integration**.

1. **Pipeline Module Table:** Clean breakdown of core viral, prokaryotic, benchmarking, and upstream modules (`preprocess`, `assembly`, `viral-identify`, `viral-taxonomy`, `viral-host`, `viral-annotate`, `prok-binning`, `checkv-pyhmmer`, etc.).
2. **Installation Section:** Covers both **Conda/Mamba** (including `conda-lock` fallback) and **Apptainer/Singularity** container setup.
3. **Quick Start Tutorial:** Provides step-by-step instructions for running `viral-identify` on sample contigs (including dry-run, high-memory/fast mode, and low-memory mode).
4. **Documentation Links & Citation:** Direct links to the full ReadTheDocs site, troubleshooting guides, and pre-print citation details.
