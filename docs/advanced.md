# Advanced Usage
## HPC Job Scheduling
One great thing about vOMIX-MEGA is that is can automatically schedule jobs for you if you use a cluster system. To do that, you will need to download and install a few extra steps through the [Snakemake Plug-in Catalouge](https://snakemake.github.io/snakemake-plugin-catalog/). Here we will take you through a few common systems, but Snakemake has a general cluster manager that will allow virtually any method to be used. 

#### SLURM

#### PBS
```bash
# install into your pre-existing conda environment
conda activate vomix
conda install -c bioconda snakemake-executor-plugin-cluster-generic=1.0.9
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

> Note: Make sure to change the queue name and your email when running the command above

#### General Cluster

## Cloud Execution 

## Docker & Containerization

## Snakemake Back-end

## Quick Updating vOMIX-MEGA

While we're working on a stable version of vOMIX-MEGA, we've made it easy to update the development version to facilitate quick bug fixes for your analysis. 

```bash
# 1) Enter your vOMIX-MEGA directory
cd vOMIX-MEGA
conda activate vomix

# 2) Copy update script to the environment bin
cp workflow/scripts/vomix_update.sh $CONDA_PREFIX/bin/
chmod +x $CONDA_PREFIX/bin/vomix_update.sh

# 3) Check the usage guide and if the script is in $PATH
vomix_update.sh -h

# 4) Run script to update current directory 
vomix_update.sh . 
```

> NOTE: The `vomix_update.sh` command will ONLY update the i) Snakefile ii) config.yml file iii) rules iv) environments v) scripts IF they have changed since your current version. It will not affect any other file in your directory including analysis. 
 

## Sample List CSV
```bash
snakemake --config module="end-to-end" datadir="sample/fastq" samplelist="sample/sample_list.csv" outdir="sample/results" --use-conda -j 4 -c 4 --latency-wait 30
```
_Options_:
```bash
--config 
   datadir 
   Points to data directory where pair end files exist or will be downloaded
   samplelist  
   Points to valid sample_list.csv file
```
The `sample_list.csv` file maps paired-end sample data to assemblies or co-assemblies. It is a comma-delimited file with four main columns:

- **sample_id:** The name of the sample
- **accession:** The SRA accession to the sample for downloading (if not locally available)
- **assembly:** The assembly name (if using co-assembly or mix-assembly)
- **R1:** The path to the local forward read of raw FASTQ files (if locally available)
- **R2:** The path to the local reverse read of raw FASTQ files (if locally available)
 
You can use `--config datadir="path/to/fastqdir"` to specify the path where the data will be downloaded or already resides. The default directory is `./fastq`.

>_1. Example A: Remotely Downloaded Samples_. vOMIX-MEGA will automatically download paired-end files and retrieve all metadata from NCBI. This is the preferred approach for rapid analysis.



```csv
sample_id,accession,assembly,R1,R2
,SRR5898936,,,
,SRR5898937,,,
,SRR5898934,,,
```

>_2. Example B: Remotely Downloaded Co-assemblies_. To co-assemble samples (currently supported by `megahit` assembler), assign the same assembly name to different sample names. This will be used consistently in downstream analysis.

```csv
sample_id,accession,assembly,R1,R2
,DRR093002,mouse-A,,
,DRR093003,mouse-A,,
,SRR7716469,mouse-B,,
,SRR7716465,mouse-B,,
,SRR7716471,mouse-B,,
,SRR5716301,mouse-C,,
,SRR5716302,mouse-C,,
```


>_3. Example C: Locally Stored Files_. To configure local files you can either directly provide the path to your reads or use the `<sample_id>_{1,2}.fastq.gz` format in your `--config datadir="./fastq"` configuration. If using the latter approach, remote files will be automatically downloaded and adopt the same naming system.
```csv
sample_id,accession,assembly,R1,R2
Sample-A,,,,
Sample-B,,,,
Sample-C,,,,
WhateverNameSamp,,,Sample-D_whatevername8.fastq.gz,Sample-D_whatevername2.fastq.gz
,SRR7716469,,,
```


>For local files, you can either:
>1. Provide full file paths in the R1 and R2 columns of `sample_list.csv`.
>2. Place the FASTQ files in the `config['datadir']` path with the `<sample_id>_{1,2}.fastq.gz` naming format.
