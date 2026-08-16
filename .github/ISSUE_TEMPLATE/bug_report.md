---
name: Bug report
about: Create a report to help us improve!
title: vOMIX-MEGA Bug
labels: bug
assignees: ''
type: Bug
---

### **Before Submitting**

- [ ] I have searched existing issues and this is not a duplicate.
- [ ] I have read the [documentation](https://vomix-mega.readthedocs.io/).
- [ ] I have checked the [FAQ](https://vomix-mega.readthedocs.io/en/latest/troubleshoot.html#faq-common-issues) for common solutions.

### **Required Diagnostic Information**

📖 **For detailed instructions on how to collect each of the diagnostic items below, please see our [Troubleshooting Guide](https://vomix-snakemake.readthedocs.io/en/latest/troubleshoot.html).** It walks you through every step, including how to find log files, generate the environment listing, and capture the full terminal output.

---

**`1. Brief Issue Description`**

A clear and concise description of what the bug is.

---

**`1. Version Information`**

- **`vOMIX-MEGA version`**: (e.g., v0.1.0) – logged at startup
- **`Snakemake version`**: (run `snakemake --version`)

---

**`2. Full Terminal Output`** 

(use `vomix --verbose`) Includes your system information and working paths. Recommend to hide your personal information and user paths. 

```
[Paste the complete stdout/stderr here]
```

---

**`3. Wrapper-Generated .sh Script`** 

Location  printed on start up as `Running Script: ${OUTDIR}/.vomix/log/<timestamp>/snakemake.sh`

```
[Paste the entire snakemake.sh script here]
```

---

**`4. Metadata JSON Files`**  

Found in `${OUTDIR}/.vomix/log/<timestamp>/`

<details>
<summary>config.json</summary>

```
[Paste config.json here]
```

</details>

<details>
<summary>assemblies.json</summary>

```
[Paste assemblies.json here]
```

</details>

<details>
<summary>samples.json</summary>

```
[Paste samples.json here]
```

</details>

---

**`4. Rule-Specific Logs`** 

Found under Snakemake output `log: ${OUTDIR}path/to/log` if a specific rule fails.

```
[Paste the rule log file mentioned in the Snakemake output here]
```

---

**`5. Conda Environment Details`**

Run `vomix <module> --list-conda-envs ...` and paste the output table below.  Then activate the environment for the rule that is failing using `conda activate path/to/envs/<full SHA>` and run `conda list`.

<details>
<summary>--list-conda-envs output</summary>

```
[Paste table here]
```

</details>

<details>
<summary>conda list for the failing rule</summary>

```
[Paste conda list output here]
```

</details>

---

**`6. Input Files`**

List your input files (samplesheet, FASTQ paths, etc.). Provide a minimal subset or a public link if possible.

```
[Paste or describe input files here]
```

---

**`7. Faulty Output Files`**

Describe or paste snippets from missing/corrupted output files.

```
[Paste or describe faulty outputs here]
```

---

**`8. Cluster / HPC Logs`** (if applicable)

Include scheduler type, submission command, and any `.out`/`.err` logs.

```
[Paste cluster logs here]
```

---

**`9. Additional Context`**

Any custom modifications, environment variables, system-specific details, or steps you have already tried.

```
[Paste any extra context here]
```

---

### **📝 Important Note on Network & Database Issues**

If this issue involves downloading databases, connecting to remote servers, or network timeouts, please be aware that this is most likely due to your institution's HPC firewall, proxy settings, or internet restrictions.

vOMIX‑MEGA **cannot override your system administration policies**. Please contact your local IT support or sysadmin to resolve connectivity issues (e.g., whitelisting domains, setting up proxies) **before** reporting network‑related errors here.

---

### **👾 Thank you for submitting an issue to vOMIX-MEGA! It helps us imporove the pipeline for our community.**
