from vomix.module import Module


class PreProcessingModule(Module):
    # snakemake --config module="preprocess" decontam-host=False outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
    name = "preprocess"

    def __init__(
        self,
        decontam_host=False,
        hasOptions=False,
        dwnld_params=None,
        pigz_paramss=None,
        fastp_params=None,
        hostile_params=None,
        hostile_aligner=None,
        hostile_aligner_params=None,
        hostile_index_name=None,
        dwnld_only=False,
    ):
        self.decontam_host = decontam_host
        self.hasOptions = hasOptions
        self.dwnld_params = dwnld_params
        self.pigz_params = pigz_paramss
        self.fastp_params = fastp_params
        self.hostile_params = hostile_params
        self.hostile_aligner = hostile_aligner
        self.hostile_aligner_params = hostile_aligner_params
        self.hostile_index_name = hostile_index_name
        self.dwnld_only = dwnld_only


class AssemblyCoAssemblyModule(Module):
    # snakemake --config module="assembly" assembler="megahit" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
    name = "assembly"

    def __init__(
        self,
        assembler=False,
        hasOptions=False,
        megahit_min_len=300,
        megahit_params=None,
        spades_params=None,
        spades_memory=250,
    ):
        self.assembler = assembler
        self.hasOptions = hasOptions
        self.megahit_min_len = megahit_min_len
        self.megahit_params = megahit_params
        self.spades_params = spades_params
        self.spades_memory = spades_memory


class ViralIdentifyModule(Module):
    # snakemake --config module="viral-identify" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
    name = "viral-identify"

    def __init__(
        self,
        hasOptions=False,
        contig_min_len=0,
        genomad_db=None,
        genomad_min_len=1500,
        genomad_params=None,
        genomad_cutoff=0.7,
        genomad_cutoff_s=0,
        checkv_original=False,
        checkv_params=None,
        checkv_database=None,
        clustering_fast=True,
        cdhit_params=None,
        vOTU_ani=95,
        vOTU_targetcov=85,
        vOTU_querycov=0,
    ):
        self.hasOptions = hasOptions
        self.contig_min_len = contig_min_len
        self.genomad_db = genomad_db
        self.genomad_min_len = genomad_min_len
        self.genomad_params = genomad_params
        self.genomad_cutoff = genomad_cutoff
        self.genomad_cutoff_s = genomad_cutoff_s
        self.checkv_original = checkv_original
        self.checkv_params = checkv_params
        self.checkv_database = checkv_database
        self.clustering_fast = clustering_fast
        self.cdhit_params = cdhit_params
        self.vOTU_ani = vOTU_ani
        self.vOTU_targetcov = vOTU_targetcov
        self.vOTU_querycov = vOTU_querycov


class ViralBenchmarkModule(Module):
    # snakemake --config module="viral-benchmark" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
    name = "viral-benchmark"

    def __init__(
        self,
        hasOptions=False,
        PhaBox2_db=None,
        genomad_db=None,
        virsorter2_db=None,
        vibrant_db=None,
        contig_min_len=0,
        genomad_min_len=1000,
        genomad_params=None,
        genomad_cutoff=0.7,
        genomad_cutoff_s=0,
        checkv_original=False,
        checkv_params=None,
        clustering_fast=True,
        cdhit_params=None,
        vOTU_ani=95,
        vOTU_targetcov=85,
        vOTU_querycov=0,
        dvf_min_len=1500,
        phamer_min_len=2000,
        dvf_params=None,
        phamer_params=None,
        virsorter2_params=None,
        vf_params=None,
        vibrant_params=None,
        seeker_params=None,
        ppr_params=None,
        dvf_cutoff=0.7,
        dvf_pval=0.05,
        phamer_pred=None,
        phamer_cutoff=0,
        vf_cutoff=0,
        virsorter2_cutoff=0,
        seeker_cutoff=0,
        ppr_cutoff=0,
        vibrant_cutoff=0,
        phabox2_db_name=None,
        phabox2_db_baselink=None,
    ):
        self.hasOptions = hasOptions
        self.PhaBox2_db = PhaBox2_db
        self.genomad_db = genomad_db
        self.virsorter2_db = virsorter2_db
        self.vibrant_db = vibrant_db
        self.contig_min_len = contig_min_len
        self.genomad_min_len = genomad_min_len
        self.genomad_params = genomad_params
        self.genomad_cutoff = genomad_cutoff
        self.genomad_cutoff_s = genomad_cutoff_s
        self.checkv_original = checkv_original
        self.checkv_params = checkv_params
        self.clustering_fast = clustering_fast
        self.cdhit_params = cdhit_params
        self.vOTU_ani = vOTU_ani
        self.vOTU_targetcov = vOTU_targetcov
        self.vOTU_querycov = vOTU_querycov
        self.dvf_min_len = dvf_min_len
        self.phamer_min_len = phamer_min_len
        self.dvf_params = dvf_params
        self.phamer_params = phamer_params
        self.virsorter2_params = virsorter2_params
        self.vf_params = vf_params
        self.vibrant_params = vibrant_params
        self.seeker_params = seeker_params
        self.ppr_params = ppr_params
        self.dvf_cutoff = dvf_cutoff
        self.dvf_pval = dvf_pval
        self.phamer_pred = phamer_pred
        self.phamer_cutoff = phamer_cutoff
        self.vf_cutoff = vf_cutoff
        self.virsorter2_cutoff = virsorter2_cutoff
        self.seeker_cutoff = seeker_cutoff
        self.ppr_cutoff = ppr_cutoff
        self.vibrant_cutoff = vibrant_cutoff
        self.phabox2_db_name = phabox2_db_name
        self.phabox2_db_baselink = phabox2_db_baselink


class ViralTaxonomyModule(Module):
    # snakemake --config module="viral-taxonomy" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
    name = "viral-taxonomy"

    def __init__(
        self,
        hasOptions=False,
        PhaBox2_db=None,
        phabox2_db_name=None,
        phabox2_db_baselink=None,
        phagcn_min_len=1500,
        phagcn_params=None,
        diamond_params=None,
        genomad_db=None,
        genomad_tax_params=None,
    ):
        self.hasOptions = hasOptions
        self.PhaBox2_db = PhaBox2_db
        self.phabox2_db_name = phabox2_db_name
        self.phabox2_db_baselink = phabox2_db_baselink
        self.phagcn_min_len = phagcn_min_len
        self.phagcn_params = phagcn_params
        self.diamond_params = diamond_params
        self.genomad_db = genomad_db
        self.genomad_params = genomad_tax_params


class ViralHostModule(Module):
    # snakemake --config module="viral-host" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
    name = "viral-host"

    def __init__(
        self,
        hasOptions=False,
        PhaBox2_db=None,
        phabox2_db_name=None,
        phabox2_db_baselink=None,
        CHERRY_params=None,
        PhaTYP_params=None,
        iphop_host=False,
        iphop_cutoff=90,
        iphop_params=None,
        iphop_db=None,
        iphop_db_version=None,
        iphop_db_basename=None,
    ):
        self.hasOptions = hasOptions
        self.PhaBox2_db = PhaBox2_db
        self.phabox2_db_name = phabox2_db_name
        self.phabox2_db_baselink = phabox2_db_baselink
        self.CHERRY_params = CHERRY_params
        self.PhaTYP_params = PhaTYP_params
        self.iphop_host = iphop_host
        self.iphop_cutoff = iphop_cutoff
        self.iphop_params = iphop_params
        self.iphop_db = iphop_db
        self.iphop_db_version = iphop_db_version
        self.iphop_db_basename = iphop_db_basename


class ViralCommunityModule(Module):
    # snakemake --config module="viral-community" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 -c 4 --latency-wait 20
    name = "viral-community"

    def __init__(self, hasOptions=False, coverm_params=None, coverm_methods=None):
        self.hasOptions = hasOptions
        self.coverm_params = coverm_params
        self.coverm_methods = coverm_methods


class ViralAnnotateModule(Module):
    # snakemake --config module="viral-annotate" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
    name = "viral-annotate"

    def __init__(
        self,
        hasOptions=False,
        PhaBox2_db=None,
        eggNOG_db=None,
        eggNOG_db_params=None,
        eggNOG_params=None,
        PhaVIP_params=None,
        virsorter2_annotate_params=None,
        dram_v_annotate_params=None,
        dram_v_distill_params=None,
        dram_db=None,
        dram_setup_params=None,
        metacerberus_db=None,
        metacerberus_setup_params=None,
        metacerberus_params=None,
        pharokka_db=None,
        pharokka_params=None,
        phabox2_db_name=None,
        phabox2_db_baselink=None,
    ):
        self.hasOptions = hasOptions
        self.PhaBox2_db = PhaBox2_db
        self.eggNOG_db = eggNOG_db
        self.eggNOG_db_params = eggNOG_db_params
        self.eggNOG_params = eggNOG_params
        self.PhaVIP_params = PhaVIP_params
        self.virsorter2_annotate_params = virsorter2_annotate_params
        self.dram_v_annotate_params = dram_v_annotate_params
        self.dram_v_distill_params = dram_v_distill_params
        self.dram_db = dram_db
        self.pharokka_params = pharokka_params
        self.phabox2_db_name = phabox2_db_name
        self.phabox2_db_baselink = phabox2_db_baselink
        self.dram_setup_params = dram_setup_params
        self.metacerberus_db = metacerberus_db
        self.metacerberus_setup_params = metacerberus_setup_params
        self.metacerberus_params = metacerberus_params
        self.pharokka_db = pharokka_db
        self.pharokka_params = pharokka_params


class ProkaryoticCommunityModule(Module):
    # snakemake --config module="prok-community" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 -c 4 --latency-wait 20
    name = "prok-community"

    def __init__(self, hasOptions=False, mpa_params=None, mpa_indexv=None):
        self.hasOptions = hasOptions
        self.mpa_params = mpa_params
        self.mpa_indexv = mpa_indexv


class ProkaryoticBinningModule(Module):
    # snakemake --config module="prok-binning" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
    name = "prok-binning"

    def __init__(
        self,
        hasOptions=False,
        checkm2_db=None,
        GTDBTk_db=None,
        GTDBTk_db_version=None,
        binning_consensus=True,
        strobealign_params=None,
        MetaBAT2_params=None,
        MaxBin2_params=None,
        CONCOCT_params=None,
        jgi_summarize_params=None,
        DASTool_params=None,
        checkm2_params=None,
        galah_params=None,
        GTDBTk_identify_params=None,
        GTDBTk_align_params=None,
        GTDBTk_classify_params=None,
        VAMB_params=None,
    ):
        self.hasOptions = hasOptions
        self.checkm2_db = checkm2_db
        self.GTDBTk_db = GTDBTk_db
        self.GTDBTk_db_version = GTDBTk_db_version
        self.binning_consensus = binning_consensus
        self.strobealign_params = strobealign_params
        self.MetaBAT2_params = MetaBAT2_params
        self.MaxBin2_params = MaxBin2_params
        self.CONCOCT_params = CONCOCT_params
        self.jgi_summarize_params = jgi_summarize_params
        self.DASTool_params = DASTool_params
        self.checkm2_params = checkm2_params
        self.galah_params = galah_params
        self.GTDBTk_identify_params = GTDBTk_identify_params
        self.GTDBTk_align_params = GTDBTk_align_params
        self.GTDBTk_classify_params = GTDBTk_classify_params
        self.VAMB_params = VAMB_params


class ProkaryoticAnnotateModule(Module):
    # snakemake --config module="prok-annotate" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 --latency-wait 20
    name = "prok-annotate"

    def __init__(self, hasOptions=False, humann_params=None, humann_db=None):
        self.hasOptions = hasOptions
        self.humann_params = humann_params
        self.humann_db = humann_db


class ViralEndToEndModule(Module):
    # snakemake --config module="viral-end-to-end" outdir="sample/results" datadir="sample/fastq" samplelist="sample/sample_list.csv" --use-conda -j 4 -c 4
    name = "viral-end-to-end"

    def __init__(self, hasOptions=False):
        self.hasOptions = hasOptions


class ClusterFastModule(Module):
    # snakemake --config module="cluster-fast" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
    name = "cluster-fast"

    def __init__(
        self,
        hasOptions=False,
        clustering_fast=True,
        cdhit_params=None,
        vOTU_ani=95,
        vOTU_targetcov=85,
        vOTU_querycov=0,
    ):
        self.hasOptions = hasOptions
        self.clustering_fast = clustering_fast
        self.cdhit_params = cdhit_params
        self.vOTU_ani = vOTU_ani
        self.vOTU_targetcov = vOTU_targetcov
        self.vOTU_querycov = vOTU_querycov


class CheckVPyHMMERModule(Module):
    # snakemake --config module="checkv-pyhmmer" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
    name = "checkv-pyhmmer"

    def __init__(
        self,
        hasOptions=False,
        checkv_original=False,
        checkv_params=None,
        checkv_database=None,
    ):
        self.hasOptions = hasOptions
        self.checkv_original = checkv_original
        self.checkv_params = checkv_params
        self.checkv_database = checkv_database


class SetupDatabaseModule(Module):
    # snakemake --config module="setup-database" fasta="sample/contigs/contigs_simulated_viral_nonviral.fasta" outdir="sample/results"  --use-conda -j 4 --latency-wait 20
    name = "setup-database"

    def __init__(
        self,
        hasOptions=False,
        hostile_index_db=None,
        PhaBox2_db=None,
        genomad_db=None,
        virsorter2_db=None,
        vibrant_db=None,
        checkv_db=None,
        iphop_db=None,
        iphop_db_version=None,
        iphop_db_basename=None,
        eggNOG_db=None,
        eggNOG_db_params=None,
        humann_db=None,
        checkm2_db=None,
        GTDBTk_db=None,
        GTDBTk_db_version=None,
        dram_db=None,
        dram_setup_params=None,
        metacerberus_db=None,
        metacerberus_setup_params=None,
        pharokka_db=None,
        phabox2_db_name=None,
        phabox2_db_baselink=None,
    ):
        self.hasOptions = hasOptions
        self.hostile_index_db = hostile_index_db
        self.PhaBox2_db = PhaBox2_db
        self.genomad_db = genomad_db
        self.virsorter2_db = virsorter2_db
        self.vibrant_db = vibrant_db
        self.checkv_db = checkv_db
        self.iphop_db = iphop_db
        self.iphop_db_version = iphop_db_version
        self.iphop_db_basename = iphop_db_basename
        self.eggNOG_db = eggNOG_db
        self.eggNOG_db_params = eggNOG_db_params
        self.humann_db = humann_db
        self.checkm2_db = checkm2_db
        self.GTDBTk_db = GTDBTk_db
        self.GTDBTk_db_version = GTDBTk_db_version
        self.dram_db = dram_db
        self.dram_setup_params = dram_setup_params
        self.metacerberus_db = metacerberus_db
        self.metacerberus_setup_params = metacerberus_setup_params
        self.pharokka_db = pharokka_db
        self.phabox2_db_name = phabox2_db_name
        self.phabox2_db_baselink = phabox2_db_baselink
