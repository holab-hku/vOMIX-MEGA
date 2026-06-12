# Run & Configuration
## Table of Contents
* 🛠️ [Preprocessing](#pre-processing)
* 🛠️ [Assembly & Co-assembly](#assembly-co-assembly) 
* 🦠 [Viral Identify](##viral-identify)
* 🦠 [Viral Taxonomy](##viral-taxonomy) 
* 🦠 [Viral Host Analysiss](#viral-host)
* 🦠 [Viral Annotate](##viral-annotate) 
* 🦠 [Viral Community](#-viral-community)
* 🧫 [Prokaryotic Community](#prokaryotic-community)
* 🧫 [Prokaryotic Binning](#-prokaryotic-binning)
* 🧫 [Prokaryotic Annotate](##prokaryotic-annotate)
* 💻 [End-to-end](#end-to-end) 
* ⚙️ [Clustering Fast (vOTU)](#clustering-fast-votu)
* ⚙️ [CheckV PyHMMER](#checkv-pyhmmer)
* ⚙️ [Setup Database](#setup-database)
* ⚙️ [Viral Benchmark Tools](#viral-benchmark-tools)

## 📑 General Configuration 
All the parameters of vOMIX-snakemake are configured using the `config/config.yml` file. You can create a custom `config.yml` by downloading the template from our GitHub page, altering the parameters, and passing it to the `--configfile config.yml` parameter. 

> NOTE: The config.yml file needs to maintain the correct formatting. We check the format of the config file using JSON schemas and will throw warnings if the proper formatting is not maintained. You can always reference the original correct formatting of [config.yml here](https://github.com/holab-hku/vOMIX-snakemake/blob/main/config/config.yml).  

Alternatively, you can pass individual parameters to snakemake using the `--config` flag, which we will demonstrate in this tutorial per module. Even though vOMIX-snakemake is rigorously benchmarked and we use the best standard settings for all modules, we've made it entirely adjustable to your use. 

_Universal Configuration:_
```
--config
    module
        Chooses the vOMIX-snakemake module to run || default: "end-to-end"
    workdir
        Set the working directory for Snakefile (We recommend not changing this)
    outdir
        Select the output directory for hierarchal results formatting || default: "./results"
    intermediate
        Flag to keep LARGE intermediate files generated during analysis || default: False
    splits
        Splits data into N chunks to reduce memory usage wherever possible || default: 0
```
_To check the full snakemake run options run:_

```h
snakemake -h
```
## 📑 Command Line Format 
Running vOMIX-snakemake with snakemake on the command line is simple. The general structure of running an analysis with vOMIX-snakemake is made up of three different components

```bash
# 1) Snakemake command
snakemake 

# 2) --config to pass vOMIX-snakemake commands 
--config module="preprocess" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" 

# 3) Snakemake execution parameters 
--use-conda -j 4 --latency-wait 20
```

_Together they make:_
```bash
snakemake --config module="preprocess" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
```
> NOTE: To view all configurations you can pass to vOMIX-snakemake using `--config`, you can check the `config/config.yml` file. To view all Snakemake execution parameters, you can run `snakemake -h`

You can use three different input types to pass on data to vOMIX-snakemake following the format:
```bash
# sample_list.csv format with datadir for fastq files
snakemake --config module="viral-identify" samplelist="sample/sample_list.csv" datadir="sample/fastq" outdir="sample/results" --use-conda -j 4 --latency-wait 20

# single fasta file input accepted in certain modules
snakemake --config module="viral-identify" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results" --use-conda -j 4 --latency-wait 20

# directory of fasta files accepted in certain modules
snakemake --config module="viral-identify" fastadir="sample/contigs/" outdir="sample/results" --use-conda -j 4 --latency-wait 20

```

> NOTE: Not all modules accept all three inputs. Please check this page for each module to see what inputs are accepted. In general, `sample_list.csv` is what you will need for comprehensive analysis, while `fasta` and `fastadir` allow each module to be used independently without the need for running previous vOMIX-snakemake steps.  

## 📑 Module-Based Analysis 

### 🛠️ Pre-processing

_Quick Run:_
```bash
snakemake --config module="preprocess" decontam-host=False outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
```
_Options:_
```
--config
    dwnldparams
        Parameters for fasterq-dump for downloading from sra tools https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump || default: ""
    pigzparams
        Parameters of pigz for compressing downloaded fastq files https://github.com/madler/pigz || default: ""
    fastpparams
        Parameters to pass on fastp software https://github.com/OpenGene/fastp || default: ""
    hostileparams
        Parameters for hostile decontamination https://github.com/bede/hostile|| default: ""
    hostilealigner
        Which mapper to use for host decontamination- bowtie2 or minimap2 (recommended) || default: "minimap2"
    alignerparams
        PLEASE DO NOT change the -x sr for minimap2 to make sure it can accurately map short reads || default: "-x sr"
    indexpath
        Path to host contamination || default: "./workflow/database/hostile/human-t2t-hla.fa.gz"
```
_Accepted Inputs:_
1. Sample List CSV ✅
2. Single Fasta ❌
3. Multi Fasta Directory ❌

### 🛠️ Assembly & Co-assembly
_Quick Run:_
```bash
snakemake --config module="assembly" assembler="megahit" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
```
_Options:_
```
--config
    megahit-minlen
        Minimum length for snakemakeHIT to use for contig building || default: 300
    megahit-params
        Extra parameters to hand off to snakemakeHIT software https://github.com/voutcn/megahit || default: "--prune-level 3"
    spades-params
        Parameters to pass on fastp software https://github.com/OpenGene/fastp || default: "--meta"
    spades-memory
        Parameters for hostile decontamination https://github.com/bede/hostile|| default: 250 NUM
```
_Accepted Inputs:_
1. Sample List CSV ✅
2. Single Fasta ❌
3. Multi Fasta Directory ❌

### 🦠 Viral Identify
_Quick Run:_
```bash
snakemake --config module="viral-identify" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
```
_Options:_
```
--config
    contig-minlen
        Minimum contig length to filter BEFORE viral identification || default: 0 [INT]
    genomad-db
        Path to geNomad's databases || default: "workflow/database/genomad" [STR]
    genomad-minlen
        Minimum viral contig length for the geNomad cutoff || default: 1500 [INT]
    genomad-params
        Additional parameters to hand off to geNomad's analysis || default: "" [STR]
    genomad-cutoff
        Parameters for hostile decontamination https://github.com/bede/hostile|| default: 0.7 [INT]
    checkv-original
        Flag to use CheckV original instead of the much faster version in vOMIX-snakemake, CheckV-PyHMMER. || default: False [True or False]
    checkv-params
        Additional parameters to pass on to CheckV. Read more at https://bitbucket.org/berkeleylab/CheckV/src || default: "" [STR]
    checkv-database
        Path to CheckV's database || default: "workflow/database/checkv" [STR]
    clustering-fast
        Flag to run fast clustering using CheckV's snakemakeBLAST approach. If set to False, CD-HIT will be used. Proceed with caution as it can be extremely slow at large sequence numbers. || default: True [True or False]
    cdhit-params
        Additional parameters to pass on to CD-HIT if clustering-fast is set to False. Read more at https://github.com/weizhongli/cdhit/blob/master/doc/cdhit-user-guide.wiki || default: "-c 0.95 -aS 0.85 -d 400 -M 0 -n 5" [STR]
    vOTU-ani
        Minimum average nucleotide identity for fast clustering algorithm of viral contigs || default: 95 [INT]
    vOTU-targetcov
        Minimum target coverage for fast clustering algorithm of viral contigs || default: 85 [NUM]
    vOTU-querycov
        Minimum query coverage for fast clustering algorithm of viral contigs || default: 0 [NUM]
```
_Accepted Inputs:_
1. Sample List CSV ✅
2. Single Fasta ✅
3. Multi Fasta Directory ✅

### 🦠 Viral Taxonomy
_Quick Run:_
```
snakemake --config module="viral-taxonomy" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
```
_Options:_
```
--config
    viphogs-hmmeval
        Minimum e value for ViPhogs hmms to be considered a hit || default: 0.01 [NUM]
    viphogs-prop
        Minimum proportion of annotated genes required for taxonomic assignment || default: 0.6 [NUM]
    PhaBox2-db
        Path to phabox database directory || default: "workflow/database/phabox_db_v2"
    phagcn-minlen
        Minimum contig length to filter before PaGCN taxonomy annotation || default: 1500 [INT]
    phagcn-params
        Additional parameters to pass on to PhaGCN || default: "" [STR]
    diamond-params
        Parameters for taxonomic classification using diamond || default: "--query-cover 50 --subject-cover 50 --evalue 1e-5 --max-target-seqs 1000" [INT] || WARNING: strongly recommend not changing this as it has been extensively tested by https://doi.org/10.1038/s41564-021-00928-6
    genomad-db
        Path to geNomad database directory || default: "workflow/database/genomad" [STR]
    genomad-params
        Additional parameters to pass on to geNomad || default: "--enable-score-calibration --relaxed" [INT]
```
_Accepted Inputs:_
1. Sample List CSV ✅
2. Single Fasta ✅
3. Multi Fasta Directory ❌

### 🦠 Viral Host
_Quick Run:_
```
snakemake --config module="viral-host" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
```
_Options:_
```
--config
    CHERRY-params
        Parameters to pass on to CHERRY for host identification. Read more at https://phage.ee.cityu.edu.hk/wiki || default: "" [STR]
    PhaTYP-params
        Parameters to pass on to MaxBin2 for lifestyle identification. Read more at https://phage.ee.cityu.edu.hk/wiki || default: "" [STR]
    iphop-cutoff
        The number of correct host predictions was evaluated for 3 different score cutoffs corresponding to 20%, 10%, and 5% estimated FDR || default: 90 [NUM]
    iphop-params
        Parameters to pass on to iPhOp for consensus host analysis. Read more at https://bitbucket.org/srouxjgi/iphop/src || default: "" [STR]

```
_Accepted Inputs:_
1. Sample List CSV ✅
2. Single Fasta ✅
3. Multi Fasta Directory ❌

### 🦠 Viral Community
_Quick Run:_
```
snakemake --config module="viral-community" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 -c 4 --latency-wait 20
```
_Options:_
```
--config
    mpa-indexv
        The version of the MetaPhlAn4 database to download || default: "mpa_vOct22_CHOCOPhlAnSGB_202212" [STR]
    mpa-params
        Additional parameters to pass on to metaphlan function. See https://huttenhower.sph.harvard.edu/metaphlan/ for more. || default: "--ignore_eukaryotes" [STR]

```
_Accepted Inputs:_
1. Sample List CSV ✅
2. Single Fasta ❌
3. Multi Fasta Directory ❌


### 🦠 Viral Annotate
_Quick Run:_
```
snakemake --config module="viral-annotate" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
```
_Options:_
```
--config
    eggNOG-params
        Parameters for running eggNOG-mapper v2. See more at https://github.com/eggnogdb/eggnog-mapper/wiki || default: "-m diamond --hmm_evalue 0.001 --hmm_score 60 --query-cover 20 --subject-cover 20 --tax_scope auto --target_orthologs all --go_evidence non-electronic --report_orthologs" [INT]
    PhaVIP-params
        Minimum contig length to filter BEFORE viral identification || default: "" [STR]

```
_Accepted Inputs:_
1. Sample List CSV ✅
2. Single Fasta ✅
3. Multi Fasta Directory ❌

### 🧫 Prokaryotic Community
_Quick Run:_
```
snakemake --config module="prok-community" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 -c 4 --latency-wait 20
```
_Options:_
```
--config
    mpa-params
        Parameters for metaphlan function. See more at https://huttenhower.sph.harvard.edu/metaphlan/ || default: "--ignore_eukaryotes" [STR]
    mpa-indexv
        Database version for metaphlan to use. See more at https://huttenhower.sph.harvard.edu/metaphlan/ || default: "mpa_vOct22_CHOCOPhlAnSGB_202212" [STR]

```
_Accepted Inputs:_
1. Sample List CSV ✅
2. Single Fasta ❌
3. Multi Fasta Directory ❌

### 🧫 Prokaryotic Binning 
_Quick Run:_
```bash
snakemake --config module="prok-binning" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
```
_Options:_
```
--config
    binning-consensus
        Flag to indicate whether multitool binning should be performed with CONCOCT, MetaBAT2, and MaxBin2 and consensus derived from DASTool OR only single-tool binning is performed with VAMB || default: True [True or False]
    strobealign-params
        Parameters to pass on to strobealign for short-read alignment. Read more at https://github.com/ksahlin/strobealign || default: "" [STR]
    MetaBAT2-params
        Parameters to pass on to MetaBAT2 for consensus binning. Read more at https://gensoft.pasteur.fr/docs/MetaBAT/2.15/ || default: "-m 1500" [STR]
    MaxBin2-params
        Parameters to pass on to MaxBin2 for consensus binning. Read more at https://flowcraft.readthedocs.io/en/latest/user/components/maxbin2.html || default: "-min_contig_length 1500 -max_iteration 50 -prob_threshold 0.9" [STR]
    CONCOCT-params
        Parameters to pass on to CONCOCT for consensus binning. Read more at https://concoct.readthedocs.io/en/stable/|| default: "" [STR]
    jgi-summarize-params
        Parameters to pass on to the jgi_summarize_bam_contig_depths function for MetaBAT2 for transforming depth samples from MaxBin2 for MetaBAT2.|| default: "--percentIdentity 97" [STR]
    DASTool-params
        Parameters to pass on to DASTool for consensus binning. Read more at https://github.com/cmks/DAS_Tool || default: "" [STR]
    drep-params
        Parameters to pass on to to dRep for dereplication of MAGs. Read more at https://drep.readthedocs.io/en/latest/overview.html  || default: "" [STR]
    VAMB-params
        Parameters to pass on to VAMB for single-tool binning. Read more at https://github.com/RasmussenLab/vamb || default: "-m 100" [STR]

```
_Accepted Inputs:_
1. Sample List CSV ✅
2. Single Fasta ❌
3. Multi Fasta Directory ❌

### 🧫 Prokaryotic Annotate
_Quick Run:_
```bash
snakemake --config module="prok-annotate" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
```
_Options:_
```
--config
```
_Accepted Inputs:_
1. Sample List CSV ✅
2. Single Fasta ✅
3. Multi Fasta Directory ❌

### 💻 End-to-end
_Quick Run:_
```
snakemake --config module="end-to-end" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 -c 4
```
_Options are same as listed above_

_Accepted Inputs:_
1. Sample List CSV ✅
2. Single Fasta ❌
3. Multi Fasta Directory ❌

### ⚙️ Cluster Fast
_Quick Run:_
```
snakemake --config module="cluster-fast" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
```
_Options:_
```
--config
    clustering-fast
        Flag to run fast clustering using CheckV's snakemakeBLAST approach. If set to False, CD-HIT will be used. Proceed with caution as it can be extremely slow at large sequence numbers. || default: True [True or False]
    cdhit-params
        Additional parameters to pass on to CD-HIT if clustering-fast is set to False. Read more at https://github.com/weizhongli/cdhit/blob/master/doc/cdhit-user-guide.wiki || default: "-c 0.95 -aS 0.85 -d 400 -M 0 -n 5" [STR]
    vOTU-ani
        Minimum average nucleotide identity for fast clustering algorithm of viral contigs || default: 95 [INT]
    vOTU-targetcov
        Minimum target coverage for fast clustering algorithm of viral contigs || default: 85 [NUM]
    vOTU-querycov
        Minimum query coverage for fast clustering algorithm of viral contigs || default: 0 [NUM]
```
_Accepted Inputs:_
1. Sample List CSV ❌
2. Single Fasta ✅
3. Multi Fasta Directory ❌

### ⚙️ CheckV-PyHMMER
_Quick Run:_
```
snakemake --config module="checkv-pyhmmer" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
```
_Options:_
```
--config
    checkv-original
        Flag to use CheckV original instead of the much faster version in vOMIX-snakemake, CheckV-PyHMMER. || default: False [True or False]
    checkv-params
        Additional parameters to pass on to CheckV. Read more at https://bitbucket.org/berkeleylab/CheckV/src || default: "" [STR]
    checkv-database
        Path to CheckV's database || default: "workflow/database/checkv" [STR]
```
_Accepted Inputs:_
1. Sample List CSV ❌
2. Single Fasta ✅
3. Multi Fasta Directory ❌

### ⚙️ Setup-Database
_Quick Run:_
```
snakemake --config module="setup-database" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
```
_Options:_
```
--config
    PhaBox2-db
        Path to PhaBox2 database for download || default: "workflow/database/phabox_db_v2" [STR]
    genomad-db
        Path to geNomad database for download || default: "workflow/database/genomad" [STR]
    checkv-db
        Path to CheckV database for download || default: "workflow/database/phabox_db_v2" [STR]
    eggNOG-db
        Path to eggNOG v2 database for download || default: "workflow/database/eggNOGv2" [STR]
    eggNOG-db-params
        Parameters for downloading eggNOG v2 database || default: "" [STR]
    virsorter2-db
        Path to VirSorter2 database for download || default: "workflow/database/virsorter2" [STR]
    iphop-db
        Path to iPHoP database for download || default: "workflow/database/iphop/Aug_2023_pub_rw" [STR]
    humann-db
        Path to HUMAnN3 databases for download || default: "workflow/database/humann" [STR]

```
_Accepted Inputs:_
1. Sample List CSV ❌
2. Single Fasta ✅
3. Multi Fasta Directory ❌


### ⚙️ Viral Benchmark


## ⚙️ Troubleshooting Guide

We have specific guidelines for troubleshooting vOMIX-snakemake so we can help you out in your analysis journey as efficiently as possible! If you run into any unexpected errors, warnings, etc. please visit our [Troubleshooting Guide](/troubleshoot.md).

