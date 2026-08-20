from vomix.module import Module


class BaseModule(Module):
    """Optimized base module that automatically assigns defaults and keyword arguments."""

    def __init__(self, hasOptions=False, **kwargs):
        self.hasOptions = hasOptions
        defaults = getattr(self, "DEFAULTS", {})

        # Assign defaults, allowing kwargs to override them
        for key, default_val in defaults.items():
            setattr(self, key, kwargs.get(key, default_val))

        # Catch any extra kwargs not explicitly declared in DEFAULTS
        for key, val in kwargs.items():
            if key not in defaults:
                setattr(self, key, val)


class PreProcessingModule(BaseModule):
    name = "preprocess"
    DEFAULTS = {
        "decontam_host": None,
        "dwnld_params": None,
        "pigz_params": None,
        "fastp_params": None,
        "hostile_params": None,
        "hostile_aligner": None,
        "hostile_aligner_params": None,
        "hostile_index_name": None,
        "dwnld_only": None,
    }


class AssemblyCoAssemblyModule(BaseModule):
    name = "assembly"
    DEFAULTS = {
        "short_read_assembler": None,
        "metamdbg_params": None,
        "nanomdbg_params": None,
        "megahit_min_len": None,
        "megahit_params": None,
        "spades_params": None,
        "spades_memory": None,
    }


class ViralIdentifyModule(BaseModule):
    name = "viral-identify"
    DEFAULTS = {
        "contig_min_len": None,
        "contig_splits": None,
        "genomad_db": None,
        "genomad_min_len": None,
        "genomad_params": None,
        "genomad_cutoff": None,
        "genomad_cutoff_s": None,
        "checkv_original": None,
        "checkv_splits": None,
        "checkv_params": None,
        "checkv_database": None,
        "clustering_fast": None,
        "cluster_iter": None,
        "cdhit_params": None,
        "vOTU_ani": None,
        "vOTU_targetcov": None,
        "vOTU_querycov": None,
    }


class ViralBenchmarkModule(BaseModule):
    name = "viral-benchmark"
    DEFAULTS = {
        "PhaBox2_db": None,
        "genomad_db": None,
        "virsorter2_db": None,
        "vibrant_db": None,
        "contig_min_len": None,
        "contig_splits": None,
        "genomad_min_len": None,
        "genomad_params": None,
        "genomad_cutoff": None,
        "dvf_min_len": None,
        "phamer_min_len": None,
        "dvf_params": None,
        "phamer_params": None,
        "virsorter2_params": None,
        "vf_params": None,
        "vibrant_params": None,
        "seeker_params": None,
        "ppr_params": None,
        "dvf_cutoff": None,
        "dvf_pval": None,
        "phamer_pred": None,
        "phamer_cutoff": None,
        "vf_cutoff": None,
        "virsorter2_cutoff": None,
        "seeker_cutoff": None,
        "ppr_cutoff": None,
        "vibrant_cutoff": None,
        "phabox2_db_name": None,
        "phabox2_db_baselink": None,
    }


class ViralTaxonomyModule(BaseModule):
    name = "viral-taxonomy"
    DEFAULTS = {
        "genomad_db": None,
        "PhaBox2_db": None,
        "phabox2_db_name": None,
        "phabox2_db_baselink": None,
        "phagcn_min_len": None,
        "phagcn_params": None,
        "diamond_params": None,
        "genomad_params_tax": None,
    }


class ViralHostModule(BaseModule):
    name = "viral-host"
    DEFAULTS = {
        "PhaBox2_db": None,
        "phabox2_db_name": None,
        "phabox2_db_baselink": None,
        "CHERRY_params": None,
        "PhaTYP_params": None,
        "iphop_host": None,
        "iphop_cutoff": None,
        "iphop_params": None,
        "iphop_db": None,
        "iphop_db_version": None,
        "iphop_db_basename": None,
    }


class ViralCommunityModule(BaseModule):
    name = "viral-community"
    DEFAULTS = {
        "coverm_sr_mapper": None,
        "coverm_pacbio_mapper": None,
        "coverm_nanopore_mapper": None,
        "coverm_params": None,
        "coverm_methods": None
    }


class ViralAnnotateModule(BaseModule):
    name = "viral-annotate"
    DEFAULTS = {
        "PhaBox2_db": None,
        "phabox2_db_name": None,
        "phabox2_db_baselink": None,
        "eggNOG_db": None,
        "eggNOG_db_params": None,
        "eggNOG_params": None,
        "PhaVIP_params": None,
        "virsorter2_annotate_params": None,
        "dram_v_annotate_params": None,
        "dram_v_distill_params": None,
        "dram_db": None,
        "dram_setup_params": None,
        "metacerberus_db": None,
        "metacerberus_setup_params": None,
        "metacerberus_params": None,
        "pharokka_db": None,
        "pharokka_params": None,
    }


class ProkaryoticCommunityModule(BaseModule):
    name = "prok-community"
    DEFAULTS = {
        "mpa_params": None,
        "mpa_indexv": None,
    }


class ProkaryoticBinningModule(BaseModule):
    name = "prok-binning"
    DEFAULTS = {
        "checkm2_db": None,
        "GTDBTk_db": None,
        "GTDBTk_db_version": None,
        "binning_consensus": None,
        "strobealign_params": None,
        "MetaBAT2_params": None,
        "MaxBin2_params": None,
        "CONCOCT_params": None,
        "jgi_summarize_params": None,
        "DASTool_params": None,
        "checkm2_params": None,
        "galah_params": None,
        "GTDBTk_identify_params": None,
        "GTDBTk_align_params": None,
        "GTDBTk_classify_params": None,
        "VAMB_params": None,
    }


class ProkaryoticAnnotateModule(BaseModule):
    name = "prok-annotate"
    DEFAULTS = {
        "humann_params": None,
        "humann_db": None,
    }


class ViralEndToEndModule(BaseModule):
    name = "viral-end-to-end"
    DEFAULTS = {
        "decontam_host": None,
        "dwnld_params": None,
        "pigz_params": None,
        "fastp_params": None,
        "hostile_params": None,
        "hostile_aligner": None,
        "hostile_aligner_params": None,
        "hostile_index_name": None,
        "dwnld_only": None,
        "megahit_min_len": None,
        "megahit_params": None,
        "spades_params": None,
        "assembler": None,
        "spades_memory": None,
        "contig_min_len": None,
        "contig_splits": None,
        "genomad_db": None,
        "genomad_min_len": None,
        "genomad_params": None,
        "genomad_cutoff": None,
        "genomad_cutoff_s": None,
        "checkv_original": None,
        "checkv_splits": None,
        "checkv_params": None,
        "checkv_database": None,
        "clustering_fast": None,
        "cluster_iter": None,
        "cdhit_params": None,
        "vOTU_ani": None,
        "vOTU_targetcov": None,
        "vOTU_querycov": None,
        "PhaBox2_db": None,
        "phabox2_db_name": None,
        "phabox2_db_baselink": None,
        "phagcn_min_len": None,
        "phagcn_params": None,
        "diamond_params": None,
        "genomad_tax_params": None,
    }


class ClusterFastModule(BaseModule):
    name = "cluster-fast"
    DEFAULTS = {
        "cluster_method": None,
        "cluster_iter": None,
        "cdhit_params": None,
        "checkv_megablast_ani": None,
        "checkv_megablast_targetcov": None,
        "checkv_megablast_querycov": None,
        "vclust_prefilter_params": None,
        "vclust_align_params": None,
        "vclust_cluster_algorithm": None,
        "vclust_cluster_params": None,
        "dnaclust_similarity": None, 
        "dnaclust_params": None, 
        "vsearch_identity": None, 
        "vsearch_params": None,
        "viridic_threshold": None, 
        "viridic_params": None
    }


class CheckVPyHMMERModule(BaseModule):
    name = "checkv-pyhmmer"
    DEFAULTS = {
        "checkv_original": None,
        "checkv_splits": None,
        "checkv_params": None,
        "checkv_database": None,
    }


class SetupDatabaseModule(BaseModule):
    name = "setup-database"
    DEFAULTS = {
        "hostile_index_db": None,
        "PhaBox2_db": None,
        "genomad_db": None,
        "virsorter2_db": None,
        "vibrant_db": None,
        "checkv_db": None,
        "humann_db": None,
        "checkm2_db": None,
        "GTDBTk_db": None,
        "GTDBTk_db_version": None,
        "phabox2_db_name": None,
        "phabox2_db_baselink": None,
        "eggNOG_db": None,
        "eggNOG_db_params": None,
        "eggNOG_params": None,
        "PhaVIP_params": None,
        "virsorter2_annotate_params": None,
        "dram_v_annotate_params": None,
        "dram_v_distill_params": None,
        "dram_db": None,
        "dram_setup_params": None,
        "metacerberus_db": None,
        "metacerberus_setup_params": None,
        "metacerberus_params": None,
        "pharokka_db": None,
        "pharokka_params": None,
        "CHERRY_params": None,
        "PhaTYP_params": None,
        "iphop_host": None,
        "iphop_cutoff": None,
        "iphop_params": None,
        "iphop_db": None,
        "iphop_db_version": None,
        "iphop_db_basename": None,
        "coverm_params": None,
        "coverm_methods": None,
    }
