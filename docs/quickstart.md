# Quick Start

vOMIX-MEGA is a big pipeline, and we know it can feel intimidating. Still, we'd like to show you otherwise! vOMIX-MEGA is here to simplify everything, remove tedious ad-hoc coding, and speed up viral metagenomic analysis by miles, all while maintaining full freedom to tune parameters to your liking. To have a quick start on your analysis, you need to: 

1. **Install vOMIX-MEGA**
2. **Configure a `sample_list.csv` file**
3. **Start your end-to-end viral metagenomic analysis!**

> **NOTE:** The `sample_list.csv` takes SRA accessions and automatically downloads and processes the data for you!

> **NOTE:** vOMIX-MEGA has been designed to use a maximum of 24Gb of memory regardless of the size of data and number of threads used. That's what makes vOMIX-MEGA awesome. You can further reduce this number by introducing the `--config splits=8` which will reduce memory use at the expense of computation time. More on that later.

***

## 1. Install vOMIX-MEGA

> **NOTE:** Before installation, please make sure that your conda is up to date for snakemake compatibility.

_1.1 Install the vOMIX-MEGA base environment:_

```bash
# Make sure conda is updated for snakemake compatibility [IMPORTANT]
# Set channel priority to strict before running vOMIX-MEGA to ensure reproducibility [IMPORTANT]
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# Install base environment
conda create -n vomix -c conda-forge snakemake=8.25.5 biopython=1.84 -y
conda activate vomix

# Verify the two essential base tools are running
snakemake -v

```

*1.2 Download the GitHub repository*

```bash
# clone from GitHub
git clone [https://github.com/holab-hku/vOMIX-MEGA](https://github.com/holab-hku/vOMIX-MEGA)
cd vOMIX-MEGA

```

*1.3 Test Viral Contig Identification using Sample Data*

```bash
# Use more memory (22 GB) but run faster (~10 mins)
snakemake --use-conda --config module="viral-identify" outdir="test_res" splits=0 fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" -j 4 --latency-wait 20     

# Use less memory (8 GB) but run slower (~40 mins)
snakemake --use-conda --config module="viral-identify" outdir="test_res" splits=8 fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" -j 4 --latency-wait 20                                 

```

*Options:*

```text
--use-conda
    Prompts snakemake to use conda to run environments [Required]

--sdm conda
    Choose conda as the environment manager [Required]

--latency-wait
    Wait-time to check for files when they are created. Min 20 seconds for CheckV-PyHMMER

--cores 
   Use at most N CPU cores/jobs in parallel. If N is omitted or 'all', the limit is set to the number of available CPU cores.

--config
    module
        Chooses the vOMIX-MEGA module to run || default: "end-to-end"
    fasta
        Single-sample fasta input for analysis. Available only in certain modules.
    outdir
        Select the output directory for hierarchal results formatting || default: "./results"
    splits
        Splits data into N chunks to reduce memory usage wherever possible || default: 0

```

> **NOTE:** To see all snakemake parameters run `snakemake -h`. Configurations nested after `--config` are all specific to vOMIX-MEGA and can only be accessed via the `config/config.yml` file. Please read more in our Run & Configuration documentation page.

---

## 2. Configure a `sample_list.csv` file

```tsv
sample_id       accession       assembly        R1      R2
        SRR5898936
        SRR5898937
        SRR5898934

```

The sample_list.csv is a comma-delimited table noting sample names and their assemblies, and pointing to the location of the paired-end sequencing files if available. The easiest way to run vOMIX-MEGA is by writing NCBI SRA accession names as noted above. For the full tutorial, please visit our Advanced Usage documentation page.

---

## 3. Start your end-to-end viral metagenomic analysis!

```bash
# Dry run to check which analysis steps will be performed
snakemake --use-conda --configfile config/config.yml --config module="end-to-end" decontam-host=False outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" -j 4 --dry-run

# Run your pipeline
snakemake --use-conda --configfile config/config.yml --config module="end-to-end" decontam-host=False outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" -j 4 --latency-wait 20

```

*Options:*

```text
--use-conda
    Prompts snakemake to use conda to run environments [Required]

--sdm conda
    Choose conda as the environment manager [Required]

--latency-wait
    Wait-time to check for files when they are created. Min 20 seconds for CheckV-PyHMMER

--cores 
   Use at most N CPU cores/jobs in parallel.

--jobs 
   Use at most N CPU cluster/cloud jobs in parallel.

--dry-run
    Dry-run flag. Show all analysis to be done without running anything

--configfile
    Specify or overwrite the config file of the workflow. 

--config 
    decontam-host=False
        Flag to decontaminate host during pre-processing. If set to True it will automatically download and use the human-t2t-hla index
    outdir
        Directory in which results will be written
    datadir
        Directory in which the paired-end fasta files will be downloaded or already reside
    samplelist
        Specify the sample_list.csv file with sample information 

```

The pipeline handles all processes. To run on a cluster, you can visit Advanced Usage for a detailed setup guide. To change the pipeline configuration, you may change the `config.yml` file or visit Run & Configuration for a detailed run and configuration explanation.

> **Note:** Conda environments for each module will be installed automatically at the beginning of the command. Databases will be downloaded for each software automatically before the command is run.

# ⚙️ Troubleshooting Guide

We have specific guidelines for troubleshooting vOMIX-MEGA so we can help you out in your analysis journey as efficiently as possible! If you run into any unexpected errors, warnings, etc. please visit our Troubleshooting Guide.

