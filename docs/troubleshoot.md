# Troubleshooting Guide

## Quick guide to troubleshooting

If your run fails or behaves unexpectedly, follow these steps in order. This helps you identify the root cause and collect the information needed to fix the problem or report it clearly.

1. **Rerun with verbose output** – Add `--verbose` (or `-v`) to your vOMIX command to print detailed Snakemake logs to the terminal.  

   ```bash
   vomix preprocess --verbose ...   # or any module
   ```

2. **Capture the full terminal output** – Save everything printed to the console (including the initial log lines shown below) to a file:

   ```bash
   vomix preprocess ... 2>&1 | tee run.log
   ```

3. **Locate the generated Snakemake script** – The exact path is printed in the logs, e.g.:

   ```bash
   Running Script: /path/to/results/.vomix/log/vomix20260816_095407/snakemake.sh
   ```

   Keep this file – it contains the exact Snakemake command that was executed.
4. **Provide the metadata JSON files** – Inside the same timestamped directory (`${OUTDIR}/.vomix/log/<timestamp>/`) you will find:
   - `config.json` – all configuration values used for the run.
   - `assemblies.json` – assembly‑related metadata.
   - `samples.json` – sample and input file details.
   Save these files; they are essential for reproducing the environment.
5. **Reproduce with a dry‑run** – To isolate issues, try a dry run (add `-n` to the Snakemake options) or run a smaller subset of samples.
6. **Search existing resources** – Check the [documentation](https://vomix-snakemake.readthedocs.io/) and search the [GitHub issues](https://github.com/holab-hku/vOMIX-MEGA/issues) for similar problems.
7. **Inspect rule commands with `--printshellcmds`** – If a specific rule is failing, run your command with `--printshellcmds` (or `-p`) to see the exact shell commands Snakemake would execute. Then:
   - Activate the Conda environment used by that rule (see step 7 in the diagnostic section below).
   - Manually run the printed shell commands inside that environment to isolate whether the issue is with the command itself or the environment setup.

   ```bash
   vomix preprocess --printshellcmds ...   # shows all shell commands
   ```

8. **Prepare to submit a bug report** - After you've done all the tests, you can submit a bug report to our [{octicon}`mark-github;0.95em` GitHub repository](https://github.com/holab-hku/vomix-snakemake/issues/new).  Here is a quick guide on how to collect the diagnostic data:

:::{dropdown} Full Guide to Collecting Diagnostic Data for Bug Reports
:open: false

When you encounter an error and need to report it, please gather the following information. This will help us resolve the issue much faster.

### 1. vOMIX version

The version is logged at the start of every run:

```
INFO     vOMIX-MEGA v 0.1.0 initialized.
```

If you don't have the log, run:

```bash
vomix --version
```

### 2. Snakemake version

```bash
snakemake --version
```

### 3. The full Snakemake command/script

Copy the entire content of the generated `snakemake.sh` file (path is printed in the logs). This file contains the exact Snakemake command with all arguments.

### 4. Full terminal output

Include the complete stdout/stderr from the failed run. Use `--verbose` to get maximum detail.

### 5. Metadata JSON files

Provide the three JSON files from `${OUTDIR}/.vomix/log/<timestamp>/`:

- `config.json`
- `assemblies.json`
- `samples.json`

These contain all the pipeline configuration and input definitions.

### 6. Rule‑specific logs (if a rule fails)

The Snakemake output will print a block like this:

```
[Sun Aug 16 09:54:20 2026]
rule preprocessing.py download fastq from SRA:
    output: /path/to/SRR5898937_1.fastq.gz, /path/to/SRR5898937_2.fastq.gz
    log: /path/to/fastq/.log/SRR5898937.log
    jobid: 7
    ...
```

Attach the log file mentioned (e.g., `SRR5898937.log`). Also include any benchmark files if present.

### 7. Conda environment listing

Run your command with `--list-conda-envs` to see all environments used:

```bash
vomix preprocess --list-conda-envs ...
```

This outputs a table:

```
environment     container       location
workflow/envs/multiqc.yml               .snakemake/conda/ca1f39b1db43af...
workflow/envs/fastp.yml                 .snakemake/conda/f55b2991d0c2c...
...
```

For the failing rule, activate its environment and list all installed packages:

```bash
conda activate .snakemake/conda/<hash>
conda list
```

### 8. Cluster logs (if applicable)

If you ran on a cluster, provide the scheduler output files (e.g., `.out` and `.err` files, or job logs from `qacct`, `sacct`, etc.). Also include the exact submission command or cluster configuration used.

### 9. Sample data

If possible, provide a minimal example that reproduces the error – a small subset of reads or a single sample. This allows us to test and fix the issue on our end.

### 10. Any extra context

Include any custom modifications you made to the configuration, environment variables, or system‑specific details (e.g., filesystem type, network mounts).
:::

---

## Log files

vOMIX‑MEGA produces several log files that are critical for debugging:

- **Terminal output** – Everything printed by Snakemake to stdout/stderr.
- **Module‑specific logs** – Located in `results/<module>/logs/` (or `results/<module>/.log/`), these contain stdout/stderr of each tool.
- **Snakemake internal logs** – In `.snakemake/log/` (inside your working directory) – includes rule execution metadata and Conda environment creation logs.
- **Run metadata** – In `${OUTDIR}/.vomix/log/<timestamp>/` you will find:
  - `config.json`
  - `assemblies.json`
  - `samples.json`
  - `snakemake.sh` – the exact command used to invoke Snakemake.
- **Benchmark logs** – If benchmarking is enabled, rule‑specific benchmark files are stored in `results/<module>/benchmarks/`.

```{admonition} Note
:class: note
If you are using a cluster executor, scheduler logs (e.g., `qstat`, `squeue`, `qacct`) may also be available. Include those when reporting cluster‑related failures.
```

---

## FAQ & Common Issues

:::{dropdown} Conda / environment creation fails
:open: false

- Ensure you activated the correct environment: `conda activate vomix`.
- If Snakemake cannot create environments, verify that `conda`, `mamba`, or your container engine is installed and available.
- For Apptainer mode, use `snakemake --sdm apptainer --use-conda` and verify your Apptainer/Singularity setup.
- **Modern Snakemake requires conda ≥ 24.7.1**. If you see `CreateCondaEnvironmentException`, update conda (see next entry).
:::

:::{dropdown} CreateCondaEnvironmentException: Conda version too old
:open: false

Snakemake 8+ enforces a minimum conda version (24.7.1) for security and reproducibility. If your system has an older version (e.g., 4.10.3), you have several options:

1. **Update conda inside your active environment** (no admin needed):

   ```bash
   conda install conda>=24.7.1
   ```

2. **Update your global conda** (if you have admin rights):

   ```bash
   conda update -n base -c defaults conda
   ```

   (Use `-c conda-forge` if you use Miniforge/Mambaforge.)
3. **Bypass the check** (only as a last resort): Edit the Snakemake source file `site-packages/snakemake/deployment/conda.py` and comment out the version‑raising block. This is not recommended for production.
:::

:::{dropdown} Rule fails with cryptic error
:open: false

- Examine the rule‑specific log file (path is printed in the Snakemake output).
- Run the rule interactively: copy the command from `snakemake.sh` and execute it manually with `--dry-run` or `--printshellcmds` to see the exact shell commands.
- Verify input file existence and permissions.
- If the rule uses a Conda environment, list its contents (see diagnostic section above).
:::

---

## {octicon}`bug;0.85em` Report a bug to us

Have any questions or you've found a bug during your analysis? Please don't hesitate to report it to us by making an issue on our [{octicon}`mark-github;0.95em` GitHub repository](https://github.com/holab-hku/vomix-snakemake/issues/new).

```{admonition} Tip
:class: tip
Before reporting, search existing issues to avoid duplicates. Also check the [documentation](https://vomix-snakemake.readthedocs.io/) for any relevant notes.
```

We appreciate your contribution to making vOMIX‑MEGA better!
