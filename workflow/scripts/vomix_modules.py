#!/usr/bin/env python
"""Module definitions for vOMIX-MEGA workflow."""

import os
import sys
from typing import List, Union

import vomix_utils
from vomix_parse_samples import parse_sample_list

# ----------------------------------------------------------------------
# Global objects (set from the Snakefile)
# ----------------------------------------------------------------------
console = None
Panel = None
config = None  # set to the Snakemake config dict
outdir = None  # set to the output directory
targets = []  # global list of target logfiles
current_module = None


class Module:
    """Base class for all workflow modules."""

    name: str = ""
    base: str = ""
    snakemake_files: List[str] = []
    target_logs: Union[str, List[str]] = []
    add_targets_condition: bool = True
    verbose_log_bool = (isinstance(config, dict)) and (config.get("verbose", False))

    def __init__(self):
        # Context attributes (set during setup)
        self.logdir = None
        self.benchmarks = None
        self.tmpd = None
        self.samples = None
        self.assemblies = None
        self.fastap = None
        self.sample_id = None
        self.assembly_ids = None

    def should_run(self) -> bool:
        raise NotImplementedError

    def setup(self) -> "Module":
        """Prepare module context: create dirs, parse inputs, store itself globally."""
        if self.verbose_log_bool:
            console.log(f"[bold cyan]Setting up module: {self.name}[/]")
            console.log(f"  base = {self.base}")
            console.log(f"  outdir = {outdir}")

        base_path = os.path.join(outdir, self.base) if self.base else outdir
        self.logdir = os.path.join(base_path, "logs")
        self.benchmarks = os.path.join(base_path, "benchmarks")
        self.tmpd = os.path.join(base_path, "tmp")

        if self.verbose_log_bool:
            console.log(f"  logdir = {self.logdir}")
            console.log(f"  benchmarks = {self.benchmarks}")
            console.log(f"  tmpd = {self.tmpd}")

        for d in [self.logdir, self.benchmarks, self.tmpd]:
            os.makedirs(d, exist_ok=True)

        # Parse inputs (subclass-specific)
        self.parse_inputs()

        # Expose the module as the current context for .smk files
        vomix_utils.current_module = self
        if self.verbose_log_bool:
            console.log(f"  current_module set to {self.name}")
        return self

    def parse_inputs(self) -> None:
        """Hook for subclasses to parse inputs. Default does nothing."""
        if self.verbose_log_bool:
            console.log("  No specific input parsing for this module.")

    def add_targets(self) -> None:
        """Add the module's target logfile(s) to the global target list."""
        if not self.add_targets_condition or not self.target_logs:
            if self.verbose_log_bool:
                console.log(f"  No targets to add for {self.name}")
            return
        if isinstance(self.target_logs, str):
            targets = [vomix_utils.relpath(self.target_logs)]
            if self.verbose_log_bool:
                console.log(f"  target (str) added = {targets}")
        else:
            targets = [vomix_utils.relpath(t) for t in self.target_logs]
            if self.verbose_log_bool:
                console.log(f"  target (list) added = {targets}")
        if self.verbose_log_bool:
            console.log(f"  Adding targets: {targets}")
        vomix_utils.set_module_target(targets)


class PreprocessModule(Module):
    name = "preprocess"
    base = "preprocess"
    snakemake_files = ["rules/preprocess.smk"]
    target_logs = "preprocess/logs/done.log"

    def should_run(self) -> bool:
        run = config.get("module") in ["preprocess", "viral-end-to-end", "run-all"]
        if run and self.verbose_log_bool:
            console.log(f"[green]Module {self.name} will run[/]")
        return run

    def parse_inputs(self) -> None:
        if self.verbose_log_bool:
            console.log(f"Parsing sample list for {self.name}...")
        parse_quiet = config.get("module") in ["viral-end-to-end", "run-all"]
        parse_verbose = config.get("verbose", False)
        self.samples, self.assemblies = parse_sample_list(
            config["samplelist"],
            config["datadir"],
            config["outdir"],
            config["NCBI-email"],
            config["NCBI-API-key"],
            config["latest-run"],
            quiet=parse_quiet,
            verbose=parse_verbose,
        )
        if self.verbose_log_bool:
            console.log(f"  samples keys: {list(self.samples.keys())[:5]}...")
            console.log(f"  assemblies keys: {list(self.assemblies.keys())[:5]}...")


class AssemblyModule(Module):
    name = "assembly"
    base = "assembly"
    snakemake_files = ["rules/assembly.smk"]
    target_logs = "assembly/logs/done.log"

    def should_run(self) -> bool:
        run = config.get("module") in ["assembly", "viral-end-to-end", "run-all"]
        if run and self.verbose_log_bool:
            console.log(f"[green]Module {self.name} will run[/]")
        return run

    def parse_inputs(self) -> None:
        if self.verbose_log_bool:
            console.log(f"Parsing sample list for {self.name}...")
        parse_quiet = config.get("module") in ["viral-end-to-end", "run-all"]
        parse_verbose = config.get("verbose", False)
        self.samples, self.assemblies = parse_sample_list(
            config["samplelist"],
            config["datadir"],
            config["outdir"],
            config["NCBI-email"],
            config["NCBI-API-key"],
            config["latest-run"],
            quiet=parse_quiet,
            verbose=parse_verbose,
        )
        self.sr_assembler = config.get("short-read-assembler", "megahit")
        if self.verbose_log_bool:
            console.log(f"  samples keys: {list(self.samples.keys())[:5]}...")
            console.log(f"  assemblies keys: {list(self.assemblies.keys())[:5]}...")


class ViralIdentifyModule(Module):
    name = "viral-identify"
    base = "identify/viral"
    snakemake_files = ["rules/viral-identify.smk"]
    target_logs = "identify/viral/logs/done.log"

    def should_run(self) -> bool:
        run = config.get("module") in ["viral-identify", "viral-end-to-end", "run-all"]
        if run and self.verbose_log_bool:
            console.log(f"[green]Module {self.name} will run[/]")
        return run

    def parse_inputs(self) -> None:
        if self.verbose_log_bool:
            console.log(f"Parsing input for {self.name}...")
        if config.get("fasta", "") != "":
            self.fastap = vomix_utils.readfasta(config["fasta"])
            self.sample_id = config.get("sample-name", "")
            self.assembly_ids = [self.sample_id]
            if self.verbose_log_bool:
                console.log(f"  fasta mode: sample_id = {self.sample_id}")
        elif config.get("fastadir", "") != "":
            self.fastap = vomix_utils.readfastadir(config["fastadir"])
            self.assembly_ids = config.get("assembly-ids", [])
            if self.verbose_log_bool:
                console.log(f"  fastadir mode: assembly_ids = {self.assembly_ids}")
        else:
            parse_quiet = config.get("module") in ["viral-end-to-end", "run-all"]
            parse_verbose = config.get("verbose", False)
            self.samples, self.assemblies = parse_sample_list(
                config["samplelist"],
                config["datadir"],
                config["outdir"],
                config["NCBI-email"],
                config["NCBI-API-key"],
                config["latest-run"],
                quiet=parse_quiet,
                verbose=parse_verbose,
            )
            self.fastap = os.path.join(
                outdir,
                "assembly",
                "samples",
                "{sample_id}",
                "output",
                "final.contigs.fa",
            )
            self.assembly_ids = list(self.assemblies.keys())
            if self.verbose_log_bool:
                console.log(f"  sample list mode: assembly_ids = {self.assembly_ids}")


class ViralAnnotateModule(Module):
    name = "viral-annotate"
    base = "annotate/viral"
    snakemake_files = ["rules/viral-annotate.smk"]
    target_logs = "annotate/viral/logs/done.log"

    def should_run(self) -> bool:
        run = config.get("module") in ["viral-annotate", "viral-end-to-end", "run-all"]
        if run and self.verbose_log_bool:
            console.log(f"[green]Module {self.name} will run[/]")
        return run

    def parse_inputs(self) -> None:
        if self.verbose_log_bool:
            console.log(f"Parsing input for {self.name}...")
        if config.get("fasta", "") != "":
            self.fastap = vomix_utils.readfasta(config["fasta"])
            self.sample_id = config.get("sample-name", "")
            self.assembly_ids = [self.sample_id]
            if self.verbose_log_bool:
                console.log(f"  fasta mode: sample_id = {self.sample_id}")
        elif config.get("fastadir", "") != "":
            console.print(
                Panel.fit(
                    f"The {self.name} module does not accept a fasta directory. "
                    "Please provide a single FASTA file via config['fasta'] or run the full viral-end-to-end workflow.",
                    title="Input Error",
                    border_style="red",
                )
            )
            sys.exit(1)
        else:
            self.fastap = os.path.join(
                outdir, "identify/viral/output/combined.final.vOTUs.fa"
            )
            self.sample_id = "combined.final.vOTUs"
            self.assembly_ids = [self.sample_id]
            if self.verbose_log_bool:
                console.log(f"  default mode: using fixed vOTUs file {self.fastap}")


class ViralTaxonomyModule(Module):
    name = "viral-taxonomy"
    base = "taxonomy/viral"
    snakemake_files = ["rules/viral-taxonomy.smk"]
    target_logs = "taxonomy/viral/logs/done.log"

    def should_run(self) -> bool:
        run = config.get("module") in ["viral-taxonomy", "viral-end-to-end", "run-all"]
        if run and self.verbose_log_bool:
            console.log(f"[green]Module {self.name} will run[/]")
        return run

    def parse_inputs(self) -> None:
        if self.verbose_log_bool:
            console.log(f"Parsing input for {self.name}...")
        if config.get("fasta", "") != "":
            self.fastap = vomix_utils.readfasta(config["fasta"])
            self.sample_id = config.get("sample-name", "")
            self.assembly_ids = [self.sample_id]
            if self.verbose_log_bool:
                console.log(f"  fasta mode: sample_id = {self.sample_id}")
        elif config.get("fastadir", "") != "":
            console.print(
                Panel.fit(
                    f"The {self.name} module does not accept a fasta directory. "
                    "Please provide a single FASTA file via config['fasta'] or run the full viral-end-to-end workflow.",
                    title="Input Error",
                    border_style="red",
                )
            )
            sys.exit(1)
        else:
            self.fastap = os.path.join(
                outdir, "identify/viral/output/combined.final.vOTUs.fa"
            )
            self.sample_id = "combined.final.vOTUs"
            self.assembly_ids = [self.sample_id]
            if self.verbose_log_bool:
                console.log(f"  default mode: using fixed vOTUs file {self.fastap}")


class ViralHostModule(Module):
    name = "viral-host"
    base = "host"
    snakemake_files = ["rules/viral-host.smk"]
    target_logs = "host/logs/done.log"

    def should_run(self) -> bool:
        run = config.get("module") in ["viral-host", "viral-end-to-end", "run-all"]
        if run and self.verbose_log_bool:
            console.log(f"[green]Module {self.name} will run[/]")
        return run

    def parse_inputs(self) -> None:
        if self.verbose_log_bool:
            console.log(f"Parsing input for {self.name}...")
        if config.get("fasta", "") != "":
            self.fastap = vomix_utils.readfasta(config["fasta"])
            self.sample_id = config.get("sample-name", "")
            self.assembly_ids = [self.sample_id]
            if self.verbose_log_bool:
                console.log(f"  fasta mode: sample_id = {self.sample_id}")
        elif config.get("fastadir", "") != "":
            console.print(
                Panel.fit(
                    f"The {self.name} module does not accept a fasta directory. "
                    "Please provide a single FASTA file via config['fasta'] or run the full viral-end-to-end workflow.",
                    title="Input Error",
                    border_style="red",
                )
            )
            sys.exit(1)
        else:
            self.fastap = os.path.join(
                outdir, "identify/viral/output/combined.final.vOTUs.fa"
            )
            self.sample_id = "combined.final.vOTUs"
            self.assembly_ids = [self.sample_id]
            if self.verbose_log_bool:
                console.log(f"  default mode: using fixed vOTUs file {self.fastap}")


class ViralCommunityModule(Module):
    name = "viral-community"
    base = "community/viral"
    snakemake_files = ["rules/viral-community.smk"]
    target_logs = "community/viral/logs/done.log"

    def should_run(self) -> bool:
        run = config.get("module") in [
            "viral-community",
            "viral-end-to-end",
            "run-all",
        ]
        if run and self.verbose_log_bool:
            console.log(f"[green]Module {self.name} will run[/]")
        return run

    def parse_inputs(self) -> None:
        if self.verbose_log_bool:
            console.log(f"Parsing sample list for {self.name}...")
        parse_quiet = config.get("module") in ["viral-end-to-end", "run-all"]
        parse_verbose = config.get("verbose", False)
        self.samples, self.assemblies = parse_sample_list(
            config["samplelist"],
            config["datadir"],
            config["outdir"],
            config["NCBI-email"],
            config["NCBI-API-key"],
            config["latest-run"],
            quiet=parse_quiet,
            verbose=parse_verbose,
        )
        self.fastap = os.path.join(
            outdir, "identify/viral/output/combined.final.vOTUs.fa"
        )
        self.sample_id = "combined.final.vOTUs"
        self.assembly_ids = [self.sample_id]
        self.sr_assembler = config.get("short-read-assembler", "megahit")
        self.methodslist = config.get("coverm-methods", "tpm rpkm").split()
        if self.verbose_log_bool:
            console.log(f"  samples keys: {list(self.samples.keys())[:5]}...")
            console.log(f"  assemblies keys: {list(self.assemblies.keys())[:5]}...")


class ViralBenchmarkModule(Module):
    name = "viral-benchmark"
    base = "identify/viral"
    snakemake_files = ["rules/viral-benchmark.smk"]
    target_logs = "identify/viral/logs/done_benchmarks.log"

    def should_run(self) -> bool:
        run = config.get("module") in ["viral-benchmark", "run-all"]
        if run and self.verbose_log_bool:
            console.log(f"[green]Module {self.name} will run[/]")
        return run

    def parse_inputs(self) -> None:
        if self.verbose_log_bool:
            console.log(f"Parsing input for {self.name}...")
        if config.get("fasta", "") != "":
            self.fastap = vomix_utils.readfasta(config["fasta"])
            self.sample_id = config.get("sample-name", "")
            self.assembly_ids = [self.sample_id]
            if self.verbose_log_bool:
                console.log(f"  fasta mode: sample_id = {self.sample_id}")
        elif config.get("fastadir", "") != "":
            self.fastap = vomix_utils.readfastadir(config["fastadir"])
            self.assembly_ids = config.get("assembly-ids", [])
            if self.verbose_log_bool:
                console.log(f"  fastadir mode: assembly_ids = {self.assembly_ids}")
        else:
            parse_quiet = config.get("module") in ["viral-end-to-end", "run-all"]
            parse_verbose = config.get("verbose", False)
            self.samples, self.assemblies = parse_sample_list(
                config["samplelist"],
                config["datadir"],
                config["outdir"],
                config["NCBI-email"],
                config["NCBI-API-key"],
                config["latest-run"],
                quiet=parse_quiet,
                verbose=parse_verbose,
            )
            self.fastap = os.path.join(
                outdir,
                "assembly",
                "samples",
                "{sample_id}",
                "output",
                "final.contigs.fa",
            )
            self.assembly_ids = list(self.assemblies.keys())
            if self.verbose_log_bool:
                console.log(f"  sample list mode: assembly_ids = {self.assembly_ids}")


class ProkBinningModule(Module):
    name = "prok-binning"
    base = "binning/prok"
    snakemake_files = ["rules/prok-binning.smk"]
    target_logs = "binning/prok/logs/done.log"

    def should_run(self) -> bool:
        run = config.get("module") in ["prok-binning", "viral-end-to-end", "run-all"]
        if run and self.verbose_log_bool:
            console.log(f"[green]Module {self.name} will run[/]")
        return run

    def parse_inputs(self) -> None:
        if self.verbose_log_bool:
            console.log(f"Parsing sample list for {self.name}...")
        parse_quiet = config.get("module") in ["viral-end-to-end", "run-all"]
        parse_verbose = config.get("verbose", False)
        self.samples, self.assemblies = parse_sample_list(
            config["samplelist"],
            config["datadir"],
            config["outdir"],
            config["NCBI-email"],
            config["NCBI-API-key"],
            config["latest-run"],
            quiet=parse_quiet,
            verbose=parse_verbose,
        )
        self.sr_assembler = config.get("short-read-assembler", "megahit")
        if self.verbose_log_bool:
            console.log(f"  samples keys: {list(self.samples.keys())[:5]}...")
            console.log(f"  assemblies keys: {list(self.assemblies.keys())[:5]}...")


class ProkCommunityModule(Module):
    name = "prok-community"
    base = "community/metaphlan"
    snakemake_files = ["rules/prok-community.smk"]
    target_logs = "community/metaphlan/logs/done.log"

    def should_run(self) -> bool:
        run = config.get("module") in ["prok-community", "run-all"]
        if run and self.verbose_log_bool:
            console.log(f"[green]Module {self.name} will run[/]")
        return run

    def parse_inputs(self) -> None:
        if self.verbose_log_bool:
            console.log("Parsing sample list for {self.name}...")
        parse_quiet = config.get("module") in ["viral-end-to-end", "run-all"]
        parse_verbose = config.get("verbose", False)
        self.samples, self.assemblies = parse_sample_list(
            config["samplelist"],
            config["datadir"],
            config["outdir"],
            config["NCBI-email"],
            config["NCBI-API-key"],
            config["latest-run"],
            quiet=parse_quiet,
            verbose=parse_verbose,
        )
        self.sr_assembler = config.get("short-read-assembler", "megahit")
        if self.verbose_log_bool:
            console.log(f"  samples keys: {list(self.samples.keys())[:5]}...")
            console.log(f"  assemblies keys: {list(self.assemblies.keys())[:5]}...")


class ProkAnnotateModule(Module):
    name = "prok-annotate"
    base = "annotate/prok"
    snakemake_files = ["rules/prok-annotate.smk"]
    target_logs = "annotate/prok/logs/done.log"

    def should_run(self) -> bool:
        run = config.get("module") in ["prok-annotate", "run-all"]
        if run and self.verbose_log_bool:
            console.log(f"[green]Module {self.name} will run[/]")
        return run

    def parse_inputs(self) -> None:
        if self.verbose_log_bool:
            console.log("Parsing sample list for {self.name}...")
        parse_quiet = config.get("module") in ["viral-end-to-end", "run-all"]
        parse_verbose = config.get("verbose", False)
        self.samples, self.assemblies = parse_sample_list(
            config["samplelist"],
            config["datadir"],
            config["outdir"],
            config["NCBI-email"],
            config["NCBI-API-key"],
            config["latest-run"],
            quiet=parse_quiet,
            verbose=parse_verbose,
        )
        self.sr_assembler = config.get("short-read-assembler", "megahit")
        if self.verbose_log_bool:
            console.log(f"  samples keys: {list(self.samples.keys())[:5]}...")
            console.log(f"  assemblies keys: {list(self.assemblies.keys())[:5]}...")


class CheckvPyhmmerModule(Module):
    name = "checkv-pyhmmer"
    base = "identify/viral"
    snakemake_files = ["rules/checkv-pyhmmer.smk"]
    target_logs = "identify/viral/logs/checkv-done.log"

    def should_run(self) -> bool:
        run = (config.get("module")) in [
            "checkv-pyhmmer",
            "viral-identify",
            "viral-end-to-end",
            "run-all",
        ]
        if run and self.verbose_log_bool:
            console.log(f"[green]Module {self.name} will run[/]")
        return run

    def parse_inputs(self) -> None:
        if self.verbose_log_bool:
            console.log(f"Parsing input for {self.name}...")
        if config.get("fasta", "") != "":
            self.fastap = vomix_utils.readfasta(config["fasta"])
            self.sample_id = config.get("sample-name", "")
            self.assembly_ids = [self.sample_id]
            if self.verbose_log_bool:
                console.log(f"  fasta mode: sample_id = {self.sample_id}")
        elif config.get("fastadir", "") != "":
            console.print(
                Panel.fit(
                    f"The {self.name} module does not accept a fasta directory. "
                    "Please provide a single FASTA file via config['fasta'] or run the full viral-end-to-end workflow.",
                    title="Input Error",
                    border_style="red",
                )
            )
            sys.exit(1)
        else:
            self.fastap = os.path.join(
                outdir, "identify/viral/output/derep/combined.viralcontigs.derep.fa"
            )
            self.sample_id = "combined.viralcontigs.derep"
            self.assembly_ids = [self.sample_id]
            if self.verbose_log_bool:
                console.log(f"  default mode: using fixed vOTUs file {self.fastap}")


class ClusterFastModule(Module):
    name = "cluster-fast"
    base = "identify/viral"
    snakemake_files = ["rules/cluster-fast.smk"]
    target_logs = "identify/viral/logs/clustering-done.log"

    def should_run(self) -> bool:
        run = (config.get("module")) in [
            "cluster-fast",
            "viral-identify",
            "viral-end-to-end",
            "run-all",
        ]
        if run and self.verbose_log_bool:
            console.log(f"[green]Module {self.name} will run[/]")
        return run

    def parse_inputs(self) -> None:
        if self.verbose_log_bool:
            console.log(f"Parsing input for {self.name}...")
        if config.get("fasta", "") != "":
            self.fastap = vomix_utils.readfasta(config["fasta"])
            self.sample_id = config.get("sample-name", "")
            self.assembly_ids = [self.sample_id]
            if self.verbose_log_bool:
                console.log(f"  fasta mode: sample_id = {self.sample_id}")
        elif config.get("fastadir", "") != "":
            console.print(
                Panel.fit(
                    f"The {self.name} module does not accept a fasta directory. "
                    "Please provide a single FASTA file via config['fasta'] or run the full viral-end-to-end workflow.",
                    title="Input Error",
                    border_style="red",
                )
            )
            sys.exit(1)
        else:
            self.fastap = os.path.join(
                outdir, "identify/viral/intermediate/scores/combined.viralcontigs.fa"
            )
            self.sample_id = "combined.viralcontigs"
            self.assembly_ids = [self.sample_id]
            if self.verbose_log_bool:
                console.log(f"  default mode: using fixed vOTUs file {self.fastap}")


class ClusterBenchmarkModule(Module):
    name = "cluster-benchmark"
    base = "cluster-benchmark"
    snakemake_files = ["rules/cluster-benchmark.smk"]
    target_logs = "cluster-benchmark/logs/done.log"

    def should_run(self) -> bool:
        run = config.get("module") == "cluster-benchmark"
        if run and self.verbose_log_bool:
            console.log(f"[green]Module {self.name} will run[/]")
        return run


class SetupDatabaseModule(Module):
    name = "setup-database"
    base = "database"
    snakemake_files = ["rules/setup-database.smk"]
    target_logs = []

    def should_run(self) -> bool:
        run = config.get("setup-database", False) or config.get("module") in [
            "setup-database",
            "run-all",
        ]
        if run and self.verbose_log_bool:
            console.log(f"[green]Module {self.name} will run[/]")
        return run

    def setup(self) -> "Module":
        """Custom setup for setup-database: use basedir and dot-directories."""
        if self.verbose_log_bool:
            console.log(f"[bold cyan]Setting up module: {self.name}[/]")
        base_path = os.path.join(config["basedir"], self.base)
        self.logdir = os.path.join(base_path, ".logs")
        self.benchmarks = os.path.join(base_path, ".benchmarks")
        self.tmpd = os.path.join(base_path, ".tmp")

        if self.verbose_log_bool:
            console.log(f"  base_path = {base_path}")
            console.log(f"  logdir = {self.logdir}")
            console.log(f"  benchmarks = {self.benchmarks}")
            console.log(f"  tmpd = {self.tmpd}")

        for d in [self.logdir, self.benchmarks, self.tmpd]:
            os.makedirs(d, exist_ok=True)

        vomix_utils.current_module = self
        return self

    def add_targets(self) -> None:
        if vomix_utils.config.get("module") == "setup-database":
            vomix_utils.set_module_target(
                os.path.join(
                    vomix_utils.config["basedir"], "database/.benchmarks/done.log"
                )
            )


class SymlinkModule(Module):
    name = "symlink"
    base = "symlink"
    snakemake_files = ["rules/symlink.smk"]
    target_logs = ".vomix/log/symlink_done.log"

    def should_run(self) -> bool:
        run = config.get("module") == "symlink"
        if run and self.verbose_log_bool:
            console.log(f"[green]Module {self.name} will run[/]")
        return run


# List of all module instances in order
ALL_MODULES = [
    PreprocessModule(),
    AssemblyModule(),
    ViralIdentifyModule(),
    ViralTaxonomyModule(),
    ViralCommunityModule(),
    ViralHostModule(),
    ViralAnnotateModule(),
    ViralBenchmarkModule(),
    CheckvPyhmmerModule(),
    ClusterBenchmarkModule(),
    ClusterFastModule(),
    ProkBinningModule(),
    ProkCommunityModule(),
    ProkAnnotateModule(),
    SetupDatabaseModule(),
    SymlinkModule(),
]
