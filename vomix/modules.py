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
        "decontam_host": False,
        "dwnld_params": None,
        "pigz_params": None,
        "fastp_params": None,
        "hostile_params": None,
        "hostile_aligner": None,
        "hostile_aligner_params": None,
        "hostile_index_name": None,
        "dwnld_only": False,
    }


class AssemblyCoAssemblyModule(BaseModule):
    name = "assembly"
    DEFAULTS = {
        "megahit_min_len": 300,
        "megahit_params": None,
        "spades_params": None,
        "assembler": None,
        "spades_memory": None,
    }


class ViralIdentifyModule(BaseModule):
    name = "viral-identify"
    DEFAULTS = {
        "contig_min_len": 0,
        "contig_splits": 0,
        "genomad_db": None,
        "genomad_min_len": 1000,
        "genomad_params": None,
        "genomad_cutoff": 0.7,
        "genomad_cutoff_s": 0,
        "checkv_original": False,
        "checkv_splits": 0,
        "checkv_params": None,
        "checkv_database": None,
        "clustering_fast": True,
        "cluster_iter": 1,
        "cdhit_params": None,
        "vOTU_ani": 95,
        "vOTU_targetcov": 85,
        "vOTU_querycov": 0,
    }


class ViralBenchmarkModule(BaseModule):
    name = "viral-benchmark"
    DEFAULTS = {
        "PhaBox2_db": None,
        "genomad_db": None,
        "virsorter2_db": None,
        "vibrant_db": None,
        "contig_min_len": 0,
        "contig_splits": 0,
        "genomad_min_len": 1000,
        "genomad_params": None,
        "genomad_cutoff": 0.7,
        "dvf_min_len": 1500,
        "phamer_min_len": 2000,
        "dvf_params": None,
        "phamer_params": None,
        "virsorter2_params": None,
        "vf_params": None,
        "vibrant_params": None,
        "seeker_params": None,
        "ppr_params": None,
        "dvf_cutoff": 0.7,
        "dvf_pval": 0.05,
        "phamer_pred": None,
        "phamer_cutoff": 0,
        "vf_cutoff": 0,
        "virsorter2_cutoff": 0,
        "seeker_cutoff": 0,
        "ppr_cutoff": 0,
        "vibrant_cutoff": 0,
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
        "phagcn_min_len": 1000,
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
        "iphop_host": False,
        "iphop_cutoff": 90,
        "iphop_params": None,
        "iphop_db": None,
        "iphop_db_version": None,
        "iphop_db_basename": None,
    }


class ViralCommunityModule(BaseModule):
    name = "viral-community"
    DEFAULTS = {
        "coverm_params": None,
        "coverm_methods": None,
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
        "binning_consensus": True,
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
        "decontam_host": False,
        "dwnld_params": None,
        "pigz_params": None,
        "fastp_params": None,
        "hostile_params": None,
        "hostile_aligner": None,
        "hostile_aligner_params": None,
        "hostile_index_name": None,
        "dwnld_only": False,
        "megahit_min_len": 300,
        "megahit_params": None,
        "spades_params": None,
        "assembler": None,
        "spades_memory": 250,
        "contig_min_len": 0,
        "contig_splits": 0,
        "genomad_db": None,
        "genomad_min_len": 1000,
        "genomad_params": None,
        "genomad_cutoff": 0.7,
        "genomad_cutoff_s": 0,
        "checkv_original": False,
        "checkv_splits": 0,
        "checkv_params": None,
        "checkv_database": None,
        "clustering_fast": True,
        "cluster_iter": 1,
        "cdhit_params": None,
        "vOTU_ani": 95,
        "vOTU_targetcov": 85,
        "vOTU_querycov": 0,
        "PhaBox2_db": None,
        "phabox2_db_name": None,
        "phabox2_db_baselink": None,
        "phagcn_min_len": 1000,
        "phagcn_params": None,
        "diamond_params": None,
        "genomad_tax_params": None,
    }


class ClusterFastModule(BaseModule):
    name = "cluster-fast"
    DEFAULTS = {
        "clustering_fast": True,
        "cluster_iter": 1,
        "cdhit_params": None,
        "vOTU_ani": 95,
        "vOTU_targetcov": 85,
        "vOTU_querycov": 0,
    }


class CheckVPyHMMERModule(BaseModule):
    name = "checkv-pyhmmer"
    DEFAULTS = {
        "checkv_original": False,
        "checkv_splits": 0,
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
        "iphop_db": None,
        "iphop_db_version": None,
        "iphop_db_basename": None,
        "eggNOG_db": None,
        "eggNOG_db_params": None,
        "humann_db": None,
        "checkm2_db": None,
        "GTDBTk_db": None,
        "GTDBTk_db_version": None,
        "dram_db": None,
        "dram_setup_params": None,
        "metacerberus_db": None,
        "metacerberus_setup_params": None,
        "pharokka_db": None,
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
