# Changelog

All notable changes to vOMIX-MEGA will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Verbose debugging** option added for debugging to back-end and most underlying script.
- **Snakemake back-end improvements** to optimize modularization.
- **Bioconda-only installatoin** added with documentation updates.
- **Updated pyproject.toml and environment.yml** files to speed up dependency resolution and make more robust build.
- **Support for multiple clustering algorithms** under vCLUSTER and the vCLUSTREE aglorithm in `cluster-fast` module:
  - VClust, Linclust, Viridic, VSearch, DNAClust, and a new "all" option to run multiple methods
- **Long‑read and single‑end (SE) read support** across all modules where feasible, with graceful error handling where not supported
- **Conda‑only and pip‑only installation** options for lighter deployments
- `--list-conda-envs` flag to list available Conda environments
- `--reset` flag to re‑run completed modules from scratch
- `--quiet` flag to reduce log verbosity when parsing sample inputs
- bioRxiv DOI reference in the repository metadata
- `--help` now correctly displays default values for all options

### Changed

- **Lazy loading** implemented in the CLI – `--help` loads instantly without importing heavy pipeline modules
- **Default Snakemake flags** now include `-j` and `-c` with sensible defaults (both set to 4)
- **Default software deployment method** changed to `conda` (was previously empty)
- **Configuration schema convertor** renamed for clarity
- **Sample parsing** made quieter for `run-all` module
- **README.md** updated with license information
- **Documentation theme** changed to dark mode for better readability
- **`docs/_build` directory** removed from version control

### Fixed

- **Critical bug** in prokaryotic rules that mixed `sample_id` and `assembly_id` – now correctly separated
- **Cyclic rule issues** in the Snakemake workflow – resolved after debugging; `viral-end-to-end` now passes dry‑run tests
- **Validation steps** in Snakefile now removed (validation handled fully by JSON Schema parser)
- **Config indentation** errors fixed (multiple occurrences)
- **`--sdm conda`** is now the default deployment method in wrapper (was incorrectly set elsewhere)

---

## [v0.2.0]

---

## [v0.1.0-beta.1] – 2026-08-15

### Added

- Initial beta release of vOMIX-MEGA
- Full end‑to‑end viral metagenomic pipeline with **11 modules**:
  - **Preprocessing**: FASTQ quality control (fastp) and host decontamination (Hostile)
  - **Assembly & Co‑assembly**: MEGAHIT and SPAdes support
  - **Viral identification**: geNomad, CheckV/CheckV‑PyHMMER, and vOTU clustering
  - **Viral benchmarking**: 12+ detection tools (DeepVirFinder, PhaMer, VirSorter2, VirFinder, Seeker, PPR‑META, VIBRANT, etc.)
  - **Viral taxonomy**: PhaGCN, PhaBOX, geNomad taxonomy
  - **Viral host prediction**: CHERRY, iPHoP, PhaTYP
  - **Viral community**: CoverM read mapping and abundance profiling
  - **Viral annotation**: eggNOG‑mapper, MetaCerberus, pharokka, PhaVIP
  - **Prokaryotic community**: MetaPhlAn profiling
  - **Prokaryotic binning**: VAMB, MetaBAT2, MaxBin2, CONCOCT with DASTool consolidation, CheckM2 quality assessment, GTDB‑Tk taxonomy
  - **Prokaryotic annotation**: HUMAnN3 functional annotation
- `setup-database` utility to download all required reference databases
- Rich CLI with `--help` support for every subcommand
- Snakemake integration with cluster execution support
- Docker/Apptainer container support for reproducibility

### Performance

- **10–1000× faster** than comparable pipelines through optimized tool selection
- **Memory footprint capped at 24 GB** for standard end‑to‑end runs

### Documentation

- Full documentation available at <https://vomix-snakemake.readthedocs.io/>
- Citation: Erfan Shekarriz, Elsa Vijendran, Joshua WK Ho. *vOMIX-MEGA: A critical speed enhancement for end‑to‑end viral metagenomics*. bioRxiv (2026)
