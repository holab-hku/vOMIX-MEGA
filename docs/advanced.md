# Advanced Usage

## {octicon}`cache;0.85em` HPC Job Scheduling

A great perk about having a Snakemake backend on vOMIX‑MEGA is that Snakemake can automatically schedule jobs for you if you use a cluster system. To do that, you will need to download and install a few extra steps through the [Snakemake Plugin Catalog](https://snakemake.github.io/snakemake-plugin-catalog/). Here we will take you through a few common systems, but Snakemake has a general cluster manager that will allow virtually any method to be used.

### {octicon}`note;0.85em` Before You Start: Snakemake Executor Plugins

All cluster and cloud execution in Snakemake is powered by **executor plugins**. The vOMIX CLI directly supports the `--cluster-generic-submit-cmd` and `--executor` flags (see sections below). For cloud executors arguments (AWS Batch, Google Life Sciences, Azure Batch, etc.), you'll pass native Snakemake options via the `--snakemake-args` flag via double quotes (see more examples bellow).

**Important steps before using any executor:**

1. **Check your Snakemake version** – Run `snakemake --version` (or `snakemake -h`) to see which version you have.
2. **Find the right plugin** – Visit the [Snakemake Executor Plugin Catalog](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor) and look for the executor that matches your environment.
3. **Match versions** – Ensure the plugin version is compatible with your Snakemake version (each plugin page states the required Snakemake version, e.g., `>=8.6`).
4. **Install the plugin** – Use `conda install -c bioconda <plugin-package>` or `pip install <plugin-package>`.
5. **Read the plugin’s documentation** – Each plugin has its own set of configuration options, environment variables, and required permissions.

```{admonition} vOMIX CLI support for executors
:class: note
The vOMIX CLI currently **only directly supports the cluster-generic executor** via the `--executor` and `--cluster-generic-submit-cmd` flags. For **cloud executors** (AWS, GCP, Azure, etc.), you must use the `--snakemake-args` flag to pass native Snakemake options. See the Cloud Execution section below for examples.
```

---

### SLURM

First, ensure you have the `cluster-generic` executor plugin installed in your conda environment:

```bash
conda activate vomix
conda install -c bioconda snakemake-executor-plugin-cluster-generic=1.0.9
```

To run your command with SLURM, you need to add additional arguments to your normal commands:

```bash
# Local Machine Conda Run
vomix viral-benchmark --sdm conda --fasta sample/contigs/contigs_simulated_viral_nonviral.fasta --outdir quick-run/results --contig-splits 0 -j 64 --latency-wait 20

# Cluster Execution Conda Run
EMAIL="your.email@example.com"
PARTITION="your_partition_name"
ACCOUNT="your_account"  # if required

vomix viral-benchmark --sdm conda --fasta sample/contigs/contigs_simulated_viral_nonviral.fasta --outdir quick-run/results --contig-splits 0 -j 64 --latency-wait 20 --executor cluster-generic --cluster-generic-submit-cmd "sbatch --parsable --job-name={rule} --output=slurm-%j.out --error=slurm-%j.err --time=120:00:00 --mem={resources.mem_mb} --cpus-per-task={threads} --partition=$PARTITION --account=$ACCOUNT --mail-user=$EMAIL --mail-type=FAIL,END"

# Check Job Scheduling
squeue -u $USER
```

```{admonition} SLURM note
:class: note
Make sure to change the partition, account (if required), and your email when running the command above. The `--parsable` flag ensures Snakemake can correctly track job IDs. Adjust `--time` and `--mem` based on your resource requirements.
```

**Key SLURM directives explained:**

| Directive | Purpose |
| ----------- | --------- |
| `--parsable` | Returns only the job ID for Snakemake to track |
| `--job-name={rule}` | Names the job after the Snakemake rule being executed |
| `--output=slurm-%j.out` | Saves stdout to a file with the job ID |
| `--error=slurm-%j.err` | Saves stderr to a file with the job ID |
| `--time=120:00:00` | Maximum runtime (adjust as needed) |
| `--mem={resources.mem_mb}` | Allocates memory based on Snakemake's resource calculation |
| `--cpus-per-task={threads}` | Allocates the number of threads requested by Snakemake |

---

### PBS / Torque

First, ensure you have the `cluster-generic` executor plugin installed in your conda environment:

```bash
conda activate vomix
conda install -c bioconda snakemake-executor-plugin-cluster-generic=1.0.9
```

To run your command with PBS, you need to add additional arguments to your normal commands:

```bash
# Local Machine Conda Run
vomix viral-benchmark --sdm conda --sdm apptainer --fasta sample/contigs/contigs_simulated_viral_nonviral.fasta --outdir quick-run/results --contig-splits 0 -j 64 --latency-wait 20

# Cluster Execution Conda Run
EMAIL="your.email@example.com"
QUEUE="cluster_queue_name"

vomix viral-benchmark --sdm conda --fasta sample/contigs/contigs_simulated_viral_nonviral.fasta --outdir quick-run/results --contig-splits 0 -j 64 --latency-wait 20 --executor cluster-generic --cluster-generic-submit-cmd "qsub -N {log} -l nodes=1:ppn={threads} -l mem={resources.mem_mb}m -l walltime=120:00:00 -M $EMAIL -q $QUEUE -o qsub.log -e qsub.log -m a"

# Check Job Scheduling
qstat
```

```{admonition} PBS note
:class: note
Make sure to change the queue name and your email when running the command above.
```

---

### General Cluster (SGE, LSF, or Custom Schedulers)

If your cluster uses a scheduler not covered above (e.g., Sun Grid Engine, LSF, or a custom system), Snakemake's `cluster-generic` executor can still handle it. The key is providing the correct submission command.

First, install the `cluster-generic` plugin:

```bash
conda activate vomix
conda install -c bioconda snakemake-executor-plugin-cluster-generic=1.0.9
```

Then, run your command with the appropriate submission command for your scheduler:

#### Sun Grid Engine (SGE)

```bash
EMAIL="your.email@example.com"
QUEUE="your_queue"

vomix viral-benchmark --sdm conda --fasta sample/contigs/contigs_simulated_viral_nonviral.fasta --outdir quick-run/results --contig-splits 0 -j 64 --latency-wait 20 --executor cluster-generic --cluster-generic-submit-cmd "qsub -V -N {rule} -o sge.out -e sge.err -pe smp {threads} -l mem_free={resources.mem_mb}m -l h_rt=120:00:00 -q $QUEUE -M $EMAIL -m abe"

# Check Job Scheduling
qstat
```

#### LSF (IBM Spectrum LSF)

```bash
EMAIL="your.email@example.com"
QUEUE="your_queue"

vomix viral-benchmark --sdm conda --fasta sample/contigs/contigs_simulated_viral_nonviral.fasta --outdir quick-run/results --contig-splits 0 -j 64 --latency-wait 20 --executor cluster-generic --cluster-generic-submit-cmd "bsub -J {rule} -o lsf.%J.out -e lsf.%J.err -n {threads} -R 'rusage[mem={resources.mem_mb}]' -W 120:00 -q $QUEUE -u $EMAIL"

# Check Job Scheduling
bjobs
```

#### Custom Script Submission

For complex or custom schedulers, you can also write a wrapper script. Create a file called `submit.sh`:

```bash
#!/bin/bash
# submit.sh - Wrapper script for custom cluster submission

# Parse arguments passed by Snakemake
# Snakemake passes: {rule}, {threads}, {resources.mem_mb}, {log}, etc.
RULE=$1
THREADS=$2
MEM=$3
LOG=$4

# Custom submission logic here
# Example: submission to a REST API, custom queue system, etc.
echo "Submitting job for rule: $RULE with $THREADS threads and ${MEM}MB memory"
```

Then, make it executable and use it:

```bash
chmod +x submit.sh

vomix viral-benchmark --sdm conda --fasta sample/contigs/contigs_simulated_viral_nonviral.fasta --outdir quick-run/results --contig-splits 0 -j 64 --latency-wait 20 --executor cluster-generic --cluster-generic-submit-cmd "./submit.sh {rule} {threads} {resources.mem_mb} {log}"
```

```{admonition} General Cluster note
:class: note
For custom schedulers, you may need to adjust the submission command to fit your system's requirements. The placeholders `{rule}`, `{threads}`, `{resources.mem_mb}`, and `{log}` are automatically populated by Snakemake.
```

---

## {octicon}`cloud;0.85em` Cloud Execution

vOMIX‑MEGA can also run natively on cloud platforms through Snakemake's cloud executor plugins. This allows you to scale your analysis to hundreds of nodes without managing a local cluster.

````{admonition} Cloud execution via vOMIX CLI
:class: important
The vOMIX CLI **does not** have built‑in flags for cloud executors. Instead, you pass native Snakemake options using the `--snakemake-args` flag. All cloud‑specific arguments (e.g., `--executor`, `--default-remote-prefix`, `--aws-batch-job-queue`, etc.) go inside the double‑quoted string of `--snakemake-args`.
Example structure:
```bash
vomix viral-benchmark ... --snakemake-args "--executor aws-batch --default-remote-prefix s3://my-bucket --aws-batch-job-queue my-queue ..."
```
````

### AWS Batch

First, install the AWS Batch executor plugin:

```bash
conda activate vomix
conda install -c bioconda snakemake-executor-plugin-aws-batch
```

Before running, ensure you have:

- AWS CLI configured with appropriate credentials (`aws configure`)
- An AWS Batch compute environment and job queue set up
- Proper IAM permissions to submit jobs

**Set your variables:**

```bash
# AWS Batch variables
JOB_QUEUE="your-batch-queue-name"
JOB_DEFINITION="your-job-definition-name"
S3_BUCKET="s3://your-bucket/vomix-output"
AWS_REGION="us-east-1"  # optional, defaults to AWS CLI configured region
```

Then run your pipeline using `--snakemake-args`:

```bash
vomix viral-benchmark --sdm conda --fasta sample/contigs/contigs_simulated_viral_nonviral.fasta --outdir quick-run/results --contig-splits 0 -j 64 --snakemake-args "--executor aws-batch --default-remote-prefix $S3_BUCKET --aws-batch-job-queue $JOB_QUEUE --aws-batch-job-definition $JOB_DEFINITION --aws-batch-region $AWS_REGION"
```

### Google Cloud Life Sciences

Install the Google Cloud Life Sciences executor plugin:

```bash
conda activate vomix
conda install -c bioconda snakemake-executor-plugin-google-life-sciences
```

**Set your variables:**

```bash
# Google Cloud variables
export GOOGLE_APPLICATION_CREDENTIALS="/path/to/your/service-account-key.json"
export GCP_PROJECT="your-project-id"
export GCP_REGION="us-central1"
export GCP_ZONE="us-central1-a"
GCS_BUCKET="gs://your-bucket/vomix-output"
```

Then run your pipeline:

```bash
vomix viral-benchmark --sdm conda --fasta sample/contigs/contigs_simulated_viral_nonviral.fasta --outdir quick-run/results --contig-splits 0 -j 64 --snakemake-args "--executor google-life-sciences --default-remote-prefix $GCS_BUCKET --google-life-sciences-location $GCP_REGION --google-life-sciences-project $GCP_PROJECT"
```

### Azure Batch

Install the Azure Batch executor plugin:

```bash
conda activate vomix
conda install -c bioconda snakemake-executor-plugin-azure-batch
```

**Set your variables:**

```bash
# Azure Batch variables
export AZURE_BATCH_ACCOUNT="your-batch-account"
export AZURE_BATCH_ACCESS_KEY="your-access-key"
export AZURE_BATCH_ACCOUNT_URL="https://your-account.region.batch.azure.com"
AZURE_CONTAINER="azure://your-container/vomix-output"
AZURE_POOL_ID="your-pool-id"
```

Then run your pipeline:

```bash
vomix viral-benchmark --sdm conda --fasta sample/contigs/contigs_simulated_viral_nonviral.fasta --outdir quick-run/results --contig-splits 0 -j 64 --snakemake-args "--executor azure-batch --default-remote-prefix $AZURE_CONTAINER --azure-batch-pool-id $AZURE_POOL_ID --azure-batch-account $AZURE_BATCH_ACCOUNT --azure-batch-account-url $AZURE_BATCH_ACCOUNT_URL"
```

```{admonition} Cloud Execution note
:class: note
Cloud execution requires additional setup (credentials, compute environments, and storage) that varies by provider. For detailed configuration, see the [Snakemake Cloud Plugin Documentation](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor) and the specific plugin page for your chosen executor.
```

---

### {octicon}`light-bulb;0.85em` Finding the Right Executor Plugin

| Provider | Executor Name | Plugin Package | Documentation |
| ---------- | --------------- | ---------------- | --------------- |
| SLURM / PBS / SGE / LSF | `cluster-generic` | `snakemake-executor-plugin-cluster-generic` | [Cluster Generic Plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/cluster-generic.html) |
| AWS Batch | `aws-batch` | `snakemake-executor-plugin-aws-batch` | [AWS Batch Plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/aws-batch.html) |
| Google Cloud Life Sciences | `google-life-sciences` | `snakemake-executor-plugin-google-life-sciences` | [Google Life Sciences Plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/google-life-sciences.html) |
| Azure Batch | `azure-batch` | `snakemake-executor-plugin-azure-batch` | [Azure Batch Plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/azure-batch.html) |
| Kubernetes | `kubernetes` | `snakemake-executor-plugin-kubernetes` | [Kubernetes Plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/kubernetes.html) |
| Google Batch | `google-batch` | `snakemake-executor-plugin-googlebatch` | [Google Batch Plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/google-batch.html) |

```{admonition} Always check version compatibility!
:class: warning
Before installing any executor plugin, always:
1. Check your Snakemake version: `snakemake --version`
2. Visit the plugin's page in the [Snakemake Plugin Catalog](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor)
3. Verify the required Snakemake version (most plugins require `>=8.6`)
4. Read the plugin's specific configuration options, environment variables, and permission requirements
```

---

### {octicon}`light-bulb;0.85em` General Tips for Cluster Execution

| Tip | Description |
| ----- | ------------- |
| **Test with a small job** | Start with `-n` (dry-run) to verify your submission command is formatted correctly. Then run with `-j 2` on a small dataset before scaling up. |
| **Monitor your jobs** | Use the scheduler's monitoring tools (`squeue`, `qstat`, `bjobs`, or cloud dashboards) to track job status. |
| **Check resource requests** | Ensure `--mem` and `--cpus-per-task` match your cluster's available resources. Over-requesting can lead to long queue times. |
| **Use `--latency-wait`** | Network filesystems can have latency – increase `--latency-wait` if jobs fail due to missing files. |
| **Store logs** | Include `-o` and `-e` directives to save scheduler logs for debugging. |
| **Set email notifications** | Most schedulers support email alerts for job completion or failure – useful for long-running jobs. |
| **Consider `--rerun-incomplete`** | If a job fails due to a node issue, this flag will automatically re-run incomplete jobs. |

```{admonition} Need help with your cluster?
:class: tip
If you're unsure about the correct submission command for your cluster, check your system's documentation or contact your system administrator. The `{rule}`, `{threads}`, and `{resources.mem_mb}` placeholders are the most commonly used variables in job submission commands.
```

## Quick Updating vOMIX-MEGA

While we're developing a stable version of vOMIX‑MEGA, we've made it easy to update the development version to facilitate quick bug fixes for your analysis.

```bash
# 1) Enter your vOMIX‑MEGA directory
cd vOMIX-MEGA
conda activate vomix

# 2) Copy update script to the environment bin
cp workflow/scripts/vomix_update.sh $CONDA_PREFIX/bin/
chmod +x $CONDA_PREFIX/bin/vomix_update.sh

# 3) Check the usage guide and if the script is in $PATH
vomix_update.sh -h

# 4) Run script to update current directory 
vomix_update.sh . 

# 5) Rebuild packages using build
pip install . 
```

```{admonition} Update script behavior
:class: note
The `vomix_update.sh` command will ONLY update the i) Snakefile ii) config.yml file iii) rules iv) environments v) scripts IF they have changed since your current version. It will not affect any other file in your directory including analysis.
```

## {octicon}`book;0.85em` Troubleshooting Guide

We have specific guidelines for troubleshooting vOMIX-MEGA --use-conda so we can help you out in your analysis journey as efficiently as possible! If you run into any unexpected errors, warnings, etc. please visit our [Troubleshooting Guide](/troubleshoot.md).

## {octicon}`bug;0.85em` Report a bug to us

Have any questions or you've found a bug during your analysis? Please don't hesitate to report it to us by making an issue on our [{octicon}`mark-github;0.95em` GitHub repository](https://github.com/holab-hku/vOMIX-MEGA/issues/new).
