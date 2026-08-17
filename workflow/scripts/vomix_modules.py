#!/usr/bin/env python
"""Module definitions for vOMIX-MEGA workflow."""

import os
import sys
from typing import List, Union

import vomix_utils
from vomix_parse_samples import parse_sample_list


class Module:
    """Base class for all workflow modules."""

    name: str = ""
    base: str = ""
    snakemake_files: List[str] = []
    target_logs: Union[str, List[str]] = []
    add_targets_condition: bool = True

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

    def _debug(self, message: str, config: dict = None) -> None:
        """Log a message if verbose is True."""
        cfg = config if config is not None else vomix_utils.config
        if cfg is not None and cfg.get("verbose", False):
            vomix_utils.console.log(message)

    def should_run(self, config: dict) -> bool:
        raise NotImplementedError

    def setup(self, config: dict, outdir: str) -> "Module":
        """Prepare module context: create dirs, parse inputs, store itself globally."""
        vomix_utils.outdir = outdir

        self._debug(f"[bold cyan]Setting up module: {self.name}[/]", config)
        self._debug(f"  base = {self.base}", config)
        self._debug(f"  outdir = {outdir}", config)

        base_path = os.path.join(outdir, self.base) if self.base else outdir
        self.logdir = os.path.join(base_path, "logs")
        self.benchmarks = os.path.join(base_path, "benchmarks")
        self.tmpd = os.path.join(base_path, "tmp")

        self._debug(f"  logdir = {self.logdir}", config)
        self._debug(f"  benchmarks = {self.benchmarks}", config)
        self._debug(f"  tmpd = {self.tmpd}", config)

        for d in [self.logdir, self.benchmarks, self.tmpd]:
            os.makedirs(d, exist_ok=True)

        # Parse inputs (subclass-specific)
        self.parse_inputs(config, outdir)

        # Expose the module as the current context for .smk files
        vomix_utils.current_module = self
        self._debug(f"  current_module set to {self.name}", config)
        return self

    def parse_inputs(self, config: dict, outdir: str) -> None:
        """Hook for subclasses to parse inputs. Default does nothing."""
        self._debug("  No specific input parsing for this module.", config)

    def add_targets(self) -> None:
        """Add the module's target logfile(s) to the global target list."""
        cfg = vomix_utils.config
        if not self.add_targets_condition or not self.target_logs:
            self._debug(f"  No targets to add for {self.name}", cfg)
            return
        if isinstance(self.target_logs, str):
            targets = [vomix_utils.relpath(self.target_logs)]
            self._debug(f"  target (str) added = {targets}", cfg)
        else:
            targets = [vomix_utils.relpath(t) for t in self.target_logs]
            self._debug(f"  target (list) added = {targets}", cfg)
        self._debug(f"  Adding targets: {targets}", cfg)
        vomix_utils.set_module_target(targets)


class PreprocessModule(Module):
    name = "preprocess"
    base = "preprocess"
    snakemake_files = ["rules/preprocess.smk"]
    target_logs = "preprocess/logs/done.log"

    def should_run(self, config: dict) -> bool:
        run = config.get("module") in ["preprocess", "viral-end-to-end", "run-all"]
        if run:
            self._debug(f"[green]Module {self.name} will run[/]", config)
        return run

    def parse_inputs(self, config: dict, outdir: str) -> None:
        self._debug("Parsing sample list for {self.name}...", config)
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
        self._debug(f"  samples keys: {list(self.samples.keys())[:5]}...", config)
        self._debug(f"  assemblies keys: {list(self.assemblies.keys())[:5]}...", config)


class AssemblyModule(Module):
    name = "assembly"
    base = "assembly"
    snakemake_files = ["rules/assembly.smk"]
    target_logs = "assembly/logs/done.log"

    def should_run(self, config: dict) -> bool:
        run = config.get("module") in ["assembly", "viral-end-to-end", "run-all"]
        if run:
            self._debug(f"[green]Module {self.name} will run[/]", config)
        return run

    def parse_inputs(self, config: dict, outdir: str) -> None:
        self._debug("Parsing sample list for {self.name}...", config)
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
        self._debug(f"  samples keys: {list(self.samples.keys())[:5]}...", config)
        self._debug(f"  assemblies keys: {list(self.assemblies.keys())[:5]}...", config)


class ViralIdentifyModule(Module):
    name = "viral-identify"
    base = "identify/viral"
    snakemake_files = [
        "rules/viral-identify.smk",
        "rules/checkv-pyhmmer.smk",
        "rules/cluster-fast.smk",
    ]
    target_logs = "identify/viral/logs/done.log"

    def should_run(self, config: dict) -> bool:
        run = config.get("module") in ["viral-identify", "viral-end-to-end", "run-all"]
        if run:
            self._debug(f"[green]Module {self.name} will run[/]", config)
        return run

    def parse_inputs(self, config: dict, outdir: str) -> None:
        self._debug(f"Parsing input for {self.name}...", config)
        if config.get("fasta", "") != "":
            self.fastap = vomix_utils.readfasta(config["fasta"])
            self.sample_id = config.get("sample-name", "")
            self.assembly_ids = [self.sample_id]
            self._debug(f"  fasta mode: sample_id = {self.sample_id}", config)
        elif config.get("fastadir", "") != "":
            self.fastap = vomix_utils.readfastadir(config["fastadir"])
            self.assembly_ids = config.get("assembly-ids", [])
            self._debug(f"  fastadir mode: assembly_ids = {self.assembly_ids}", config)
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
            self._debug(
                f"  sample list mode: assembly_ids = {self.assembly_ids}", config
            )


class ViralAnnotateModule(Module):
    name = "viral-annotate"
    base = "annotate/viral"
    snakemake_files = ["rules/viral-annotate.smk"]
    target_logs = "annotate/viral/logs/done.log"

    def should_run(self, config: dict) -> bool:
        run = config.get("module") in ["viral-annotate", "viral-end-to-end", "run-all"]
        if run:
            self._debug(f"[green]Module {self.name} will run[/]", config)
        return run

    def parse_inputs(self, config: dict, outdir: str) -> None:
        self._debug(f"Parsing input for {self.name}...", config)
        if config.get("fasta", "") != "":
            self.fastap = vomix_utils.readfasta(config["fasta"])
            self.sample_id = config.get("sample-name", "")
            self.assembly_ids = [self.sample_id]
            self._debug(f"  fasta mode: sample_id = {self.sample_id}", config)
        elif config.get("fastadir", "") != "":
            vomix_utils.console.print(
                vomix_utils.Panel.fit(
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
            self._debug(f"  default mode: using fixed vOTUs file {self.fastap}", config)


class ViralTaxonomyModule(Module):
    name = "viral-taxonomy"
    base = "taxonomy/viral"
    snakemake_files = ["rules/viral-taxonomy.smk"]
    target_logs = "taxonomy/viral/logs/done.log"

    def should_run(self, config: dict) -> bool:
        run = config.get("module") in ["viral-taxonomy", "viral-end-to-end", "run-all"]
        if run:
            self._debug(f"[green]Module {self.name} will run[/]", config)
        return run

    def parse_inputs(self, config: dict, outdir: str) -> None:
        self._debug(f"Parsing input for {self.name}...", config)
        if config.get("fasta", "") != "":
            self.fastap = vomix_utils.readfasta(config["fasta"])
            self.sample_id = config.get("sample-name", "")
            self.assembly_ids = [self.sample_id]
            self._debug(f"  fasta mode: sample_id = {self.sample_id}", config)
        elif config.get("fastadir", "") != "":
            vomix_utils.console.print(
                vomix_utils.Panel.fit(
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
            self._debug(f"  default mode: using fixed vOTUs file {self.fastap}", config)


class ViralHostModule(Module):
    name = "viral-host"
    base = "host"
    snakemake_files = ["rules/viral-host.smk"]
    target_logs = "host/logs/done.log"

    def should_run(self, config: dict) -> bool:
        run = config.get("module") in ["viral-host", "viral-end-to-end", "run-all"]
        if run:
            self._debug(f"[green]Module {self.name} will run[/]", config)
        return run

    def parse_inputs(self, config: dict, outdir: str) -> None:
        self._debug(f"Parsing input for {self.name}...", config)
        if config.get("fasta", "") != "":
            self.fastap = vomix_utils.readfasta(config["fasta"])
            self.sample_id = config.get("sample-name", "")
            self.assembly_ids = [self.sample_id]
            self._debug(f"  fasta mode: sample_id = {self.sample_id}", config)
        elif config.get("fastadir", "") != "":
            vomix_utils.console.print(
                vomix_utils.Panel.fit(
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
            self._debug(f"  default mode: using fixed vOTUs file {self.fastap}", config)


class ViralCommunityModule(Module):
    name = "viral-community"
    base = "community/viral"
    snakemake_files = ["rules/viral-community.smk"]
    target_logs = "community/viral/logs/done.log"

    def should_run(self, config: dict) -> bool:
        run = config.get("module") in [
            "viral-community",
            "viral-end-to-end",
            "run-all",
        ]
        if run:
            self._debug(f"[green]Module {self.name} will run[/]", config)
        return run

    def parse_inputs(self, config: dict, outdir: str) -> None:
        self._debug("Parsing sample list for {self.name}...", config)
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
        self._debug(f"  samples keys: {list(self.samples.keys())[:5]}...", config)
        self._debug(f"  assemblies keys: {list(self.assemblies.keys())[:5]}...", config)


class ViralBenchmarkModule(Module):
    name = "viral-benchmark"
    base = "identify/viral"
    snakemake_files = ["rules/viral-benchmark.smk"]
    target_logs = "identify/viral/logs/done_benchmarks.log"

    def should_run(self, config: dict) -> bool:
        run = config.get("module") in ["viral-benchmark", "run-all"]
        if run:
            self._debug(f"[green]Module {self.name} will run[/]", config)
        return run

    def parse_inputs(self, config: dict, outdir: str) -> None:
        self._debug(f"Parsing input for {self.name}...", config)
        if config.get("fasta", "") != "":
            self.fastap = vomix_utils.readfasta(config["fasta"])
            self.sample_id = config.get("sample-name", "")
            self.assembly_ids = [self.sample_id]
            self._debug(f"  fasta mode: sample_id = {self.sample_id}", config)
        elif config.get("fastadir", "") != "":
            self.fastap = vomix_utils.readfastadir(config["fastadir"])
            self.assembly_ids = config.get("assembly-ids", [])
            self._debug(f"  fastadir mode: assembly_ids = {self.assembly_ids}", config)
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
            self._debug(
                f"  sample list mode: assembly_ids = {self.assembly_ids}", config
            )


class ViralRefilterModule(Module):
    name = "viral-refilter"
    base = "identify/viral"
    snakemake_files = [
        "rules/refilter-genomad.smk",
        "rules/checkv-pyhmmer.smk",
        "rules/cluster-fast.smk",
    ]
    target_logs = [
        "identify/viral/logs/done.log",
        "identify/viral/logs/clustering-done.log",
        "identify/viral/logs/checkv-done.log",
    ]

    def should_run(self, config: dict) -> bool:
        run = config.get("module") == "viral-refilter"
        if run:
            self._debug(f"[green]Module {self.name} will run[/]", config)
        return run

    def parse_inputs(self, config: dict, outdir: str) -> None:
        self._debug(f"Parsing input for {self.name}...", config)
        if config.get("fasta", "") != "":
            self.fastap = vomix_utils.readfasta(config["fasta"])
            self.sample_id = config.get("sample-name", "")
            self.assembly_ids = [self.sample_id]
            self._debug(f"  fasta mode: sample_id = {self.sample_id}", config)
        elif config.get("fastadir", "") != "":
            self.fastap = vomix_utils.readfastadir(config["fastadir"])
            self.assembly_ids = config.get("assembly-ids", [])
            self._debug(f"  fastadir mode: assembly_ids = {self.assembly_ids}", config)
        else:
            self.fastap = os.path.join(
                outdir, "identify/viral/output/derep/combined.viralcontigs.derep.fa"
            )
            self.sample_id = "combined.viralcontigs.derep"
            self.assembly_ids = [self.sample_id]
            self._debug(
                f"  sample list mode: assembly_ids = {self.assembly_ids}", config
            )


class ProkBinningModule(Module):
    name = "prok-binning"
    base = "binning/prok"
    snakemake_files = ["rules/prok-binning.smk"]
    target_logs = "binning/prok/logs/done.log"

    def should_run(self, config: dict) -> bool:
        run = config.get("module") in ["prok-binning", "viral-end-to-end", "run-all"]
        if run:
            self._debug(f"[green]Module {self.name} will run[/]", config)
        return run

    def parse_inputs(self, config: dict, outdir: str) -> None:
        self._debug("Parsing sample list for {self.name}...", config)
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
        self._debug(f"  samples keys: {list(self.samples.keys())[:5]}...", config)
        self._debug(f"  assemblies keys: {list(self.assemblies.keys())[:5]}...", config)


class ProkCommunityModule(Module):
    name = "prok-community"
    base = "community/metaphlan"
    snakemake_files = ["rules/prok-community.smk"]
    target_logs = "community/metaphlan/logs/done.log"

    def should_run(self, config: dict) -> bool:
        run = config.get("module") in ["prok-community", "run-all"]
        if run:
            self._debug(f"[green]Module {self.name} will run[/]", config)
        return run

    def parse_inputs(self, config: dict, outdir: str) -> None:
        self._debug("Parsing sample list for {self.name}...", config)
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
        self._debug(f"  samples keys: {list(self.samples.keys())[:5]}...", config)
        self._debug(f"  assemblies keys: {list(self.assemblies.keys())[:5]}...", config)


class ProkAnnotateModule(Module):
    name = "prok-annotate"
    base = "annotate/prok"
    snakemake_files = ["rules/prok-annotate.smk"]
    target_logs = "annotate/prok/logs/done.log"

    def should_run(self, config: dict) -> bool:
        run = config.get("module") in ["prok-annotate", "run-all"]
        if run:
            self._debug(f"[green]Module {self.name} will run[/]", config)
        return run

    def parse_inputs(self, config: dict, outdir: str) -> None:
        self._debug("Parsing sample list for {self.name}...", config)
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
        self._debug(f"  samples keys: {list(self.samples.keys())[:5]}...", config)
        self._debug(f"  assemblies keys: {list(self.assemblies.keys())[:5]}...", config)


class CheckvPyhmmerModule(Module):
    name = "checkv-pyhmmer"
    base = "identify/viral"
    snakemake_files = ["rules/checkv-pyhmmer.smk"]
    target_logs = "identify/viral/logs/checkv-done.log"

    def should_run(self, config: dict) -> bool:
        run = config.get("module") == "checkv-pyhmmer"
        if run:
            self._debug(f"[green]Module {self.name} will run[/]", config)
        return run

    def parse_inputs(self, config: dict, outdir: str) -> None:
        self._debug(f"Parsing input for {self.name}...", config)
        if config.get("fasta", "") != "":
            self.fastap = vomix_utils.readfasta(config["fasta"])
            self.sample_id = config.get("sample-name", "")
            self.assembly_ids = [self.sample_id]
            self._debug(f"  fasta mode: sample_id = {self.sample_id}", config)
        elif config.get("fastadir", "") != "":
            vomix_utils.console.print(
                vomix_utils.Panel.fit(
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
            self._debug(f"  default mode: using fixed vOTUs file {self.fastap}", config)


class ClusterFastModule(Module):
    name = "cluster-fast"
    base = "identify/viral"
    snakemake_files = ["rules/cluster-fast.smk"]
    target_logs = "identify/viral/logs/clustering-done.log"

    def should_run(self, config: dict) -> bool:
        run = config.get("module") == "cluster-fast"
        if run:
            self._debug(f"[green]Module {self.name} will run[/]", config)
        return run

    def parse_inputs(self, config: dict, outdir: str) -> None:
        self._debug(f"Parsing input for {self.name}...", config)
        if config.get("fasta", "") != "":
            self.fastap = vomix_utils.readfasta(config["fasta"])
            self.sample_id = config.get("sample-name", "")
            self.assembly_ids = [self.sample_id]
            self._debug(f"  fasta mode: sample_id = {self.sample_id}", config)
        elif config.get("fastadir", "") != "":
            vomix_utils.console.print(
                vomix_utils.Panel.fit(
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
            self._debug(f"  default mode: using fixed vOTUs file {self.fastap}", config)


class ClusterBenchmarkModule(Module):
    name = "cluster-benchmark"
    base = "cluster-benchmark"
    snakemake_files = ["rules/cluster-benchmark.smk"]
    target_logs = "cluster-benchmark/logs/done.log"

    def should_run(self, config: dict) -> bool:
        run = config.get("module") == "cluster-benchmark"
        if run:
            self._debug(f"[green]Module {self.name} will run[/]", config)
        return run


class SetupDatabaseModule(Module):
    name = "setup-database"
    base = "database"
    snakemake_files = ["rules/setup-database.smk"]
    target_logs = []

    def should_run(self, config: dict) -> bool:
        run = config.get("setup-database", False) or config.get("module") in [
            "setup-database",
            "run-all",
        ]
        if run:
            self._debug(f"[green]Module {self.name} will run[/]", config)
        return run

    def setup(self, config: dict, outdir: str) -> "Module":
        """Custom setup for setup-database: use basedir and dot-directories."""
        self._debug(f"[bold cyan]Setting up module: {self.name}[/]", config)
        base_path = os.path.join(config["basedir"], self.base)
        self.logdir = os.path.join(base_path, ".logs")
        self.benchmarks = os.path.join(base_path, ".benchmarks")
        self.tmpd = os.path.join(base_path, ".tmp")

        self._debug(f"  base_path = {base_path}", config)
        self._debug(f"  logdir = {self.logdir}", config)
        self._debug(f"  benchmarks = {self.benchmarks}", config)
        self._debug(f"  tmpd = {self.tmpd}", config)

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

    def should_run(self, config: dict) -> bool:
        run = config.get("module") == "symlink"
        if run:
            self._debug(f"[green]Module {self.name} will run[/]", config)
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
    ViralRefilterModule(),
    CheckvPyhmmerModule(),
    ClusterBenchmarkModule(),
    ClusterFastModule(),
    ProkBinningModule(),
    ProkCommunityModule(),
    ProkAnnotateModule(),
    SetupDatabaseModule(),
    SymlinkModule(),
]
