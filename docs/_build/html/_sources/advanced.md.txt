# Advanced Usage

## HPC Job Scheduling

One great thing about vomix-snakemake is that is can automatically schedule jobs for you if you use a cluster system. To do that, you will need to download and install a few extra steps through the [Snakemake Plug-in Catalouge](https://snakemake.github.io/snakemake-plugin-catalog/). Here we will take you through a few common systems, but Snakemake has a general cluster manager that will allow virtually any method to be used.

### SLURM

### PBS

```bash
# install into your pre-existing conda environment
conda activate vomix
conda install -c bioconda snakemake --use-conda-executor-plugin-cluster-generic=1.0.9
```

To run your command with PBS, you need to add additional arguments to your normal commands:

```bash
# Local Machine Run
snakemake --config module="viral-taxonomy" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
# Cluster Execution
snakemake --config module="viral-taxonomy" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20 --retries 3 --executor cluster-generic --cluster-generic-submit-cmd "qsub -N {log} -l nodes=1:ppn={threads} -l mem={resources.mem_mb}m -l walltime=120:00:00 -M youremail@gmail.com -q your_queu_name -o qsub.log -e qsub.log -m a"
# Check Job Scheudling
qstat
```

```{admonition} PBS note
:class: note
Make sure to change the queue name and your email when running the command above.
```

### General Cluster

## Cloud Execution

## Quick Updating

While we're developing a stable version of vomix-snakemake --use-conda, we've made it easy to update the development version to facilitate quick bug fixes for your analysis.

```bash
# 1) Enter your vomix-snakemake --use-conda directory
cd vomix-snakemake --use-conda
conda activate vomix

# 2) Copy update script to the environment bin
cp workflow/scripts/vomix_update.sh $CONDA_PREFIX/bin/
chmod +x $CONDA_PREFIX/bin/vomix_update.sh

# 3) Check the usage guide and if the script is in $PATH
vomix_update.sh -h

# 4) Run script to update current directory 
vomix_update.sh . 
```

```{admonition} Update script behavior
:class: note
The `vomix_update.sh` command will ONLY update the i) Snakefile ii) config.yml file iii) rules iv) environments v) scripts IF they have changed since your current version. It will not affect any other file in your directory including analysis.
```
