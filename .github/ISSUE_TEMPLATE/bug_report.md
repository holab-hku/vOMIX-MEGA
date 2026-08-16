---
name: Bug report
about: Create a report to help us improve!
title: vOMIX-MEGA Bug
labels: bug
assignees: ''
type: Bug
---

**Before Submitting**

- [ ] I have searched existing issues and this is not a duplicate.
- [ ] I have read the [documentation](https://vomix-mega.readthedocs.io/).
- [ ] I have checked the [FAQ](https://vomix-mega.readthedocs.io/en/latest/troubleshoot.html#faq-common-issues) for common solutions.

---

**Brief Issue Description**
A clear and concise description of what the bug is.

---

**Version Information**

- **vOMIX-MEGA version**: (e.g., v0.1.0) – logged at startup: `INFO     vOMIX-MEGA v X.Y.Z initialized.`
- **Snakemake version**: (run `snakemake --version`)

---

**Wrapper‑Generated Snakemake Script**
Paste the **entire** content of the `.sh` script that the vOMIX wrapper generated.  
Location: `${OUTDIR}/.vomix/log/<timestamp>/snakemake.sh` (as logged in the vOMIX output, e.g., `Running Script: /path/to/.../snakemake.sh`).

```

[Paste your .sh script here]

```

---

**Full Terminal Output**
Paste the **complete stdout/stderr** from your terminal (including the initial vOMIX logs and the Snakemake output).  
Use the `--verbose` flag to increase verbosity.

```

[Paste terminal output here]

```

---

**Metadata JSON Files**
Attach or paste the contents of the three JSON files found in `${OUTDIR}/.vomix/log/<timestamp>/`:

- `config.json`
- `assemblies.json`
- `samples.json`

If you can provide them as separate files, that is ideal. Otherwise, paste their content here.

<details>
<summary>config.json</summary>

```

[Paste config.json content]

```

</details>

<details>
<summary>assemblies.json</summary>

```

[Paste assemblies.json content]

```

</details>

<details>
<summary>samples.json</summary>

```

[Paste samples.json content]

```

</details>

---

**Rule‑Specific Logs (if a rule fails)**
When a rule fails, Snakemake prints a block like:

```

[date]
rule <rule_name>:
    output: /path/to/output
    log: /path/to/rule.log
    jobid: N
    ...

```

Attach the log file mentioned (e.g., `rule.log`) and any benchmark files from `results/<module>/benchmarks/`.  
If you can, paste the relevant log content here.

```

[Paste rule log content]

```

---

**Conda Environment Listing**
Run your command with `--list-conda-envs` to list all Conda environments used:

```bash
vomix <module> --list-conda-envs ...
```

Paste the output table below, and for the failing rule, activate its environment (location from the table) and run `conda list`. Provide the list here.

<details>
<summary>--list-conda-envs output</summary>

```
[Paste table here]
```

</details>

<details>
<summary>conda list for the failing rule</summary>

```
[Paste conda list output]
```

</details>

---

**Input Files**
List your input files (samplesheet, FASTQ paths, etc.). Provide a minimal subset or a public link if possible.

---

**Faulty Output Files**
Describe or paste snippets from missing/corrupted output files.

---

**Cluster / HPC Information (if applicable)**

- Scheduler type (SLURM, SGE, PBS, etc.)
- Job submission command used
- Scheduler logs (e.g., `.out`/`.err` files, `sacct` output)

```
[Paste relevant cluster logs]
```

---

**Additional Context**
Any custom modifications to the configuration, environment variables, system‑specific details (e.g., filesystem type, network mounts), or steps you have already tried.

---

**Important Note on Network & Database Issues**

If this issue involves downloading databases, connecting to remote servers, or network timeouts, please be aware that this is most likely due to your institution's HPC firewall, proxy settings, or internet restrictions.

vOMIX‑MEGA **cannot override your system administration policies**. Please contact your local IT support or sysadmin to resolve connectivity issues (e.g., whitelisting domains, setting up proxies) **before** reporting network‑related errors here.
