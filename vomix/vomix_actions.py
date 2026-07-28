import subprocess
import pty
import os
import sys
from subprocess import Popen, PIPE, CalledProcessError
import datetime
import json
import yaml
import shutil
import logging
from inspect import getsourcefile
from os.path import abspath
import inspect
import re
from rich.logging import RichHandler
from rich.console import Console

# ---------------------------------------------------------
# Rich Logging Configuration (Matching vomix.py)
# ---------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, markup=True, show_path=False)],
)

log = logging.getLogger("rich")
console = Console()


class vomix_actions:
    def __init__(self):
        self.name = "vomix"
        self.version = "1.0.0"
        self.description = "vomix is a tool for viral metagenomics analysis."

    def __repr__(self):
        return f"vomix(name={self.name}, version={self.version}, description={self.description})"

    def __str__(self):
        return f"{self.name} v{self.version}: {self.description}"

    @staticmethod
    def get_snakefile():
        filename = inspect.getframeinfo(inspect.currentframe()).filename
        base_dir = os.path.dirname(os.path.abspath(filename))
        sf = base_dir.replace(
            "/.venv/lib/python3.9/site-packages/vomix",
            "/vomix/workflow/Snakefile",
        )

        if not os.path.exists(sf):
            log.error(
                f"Unable to locate the Snakemake file; tried [bold red]{sf}[/bold red]"
            )
            sys.exit(1)
        return sf

    @staticmethod
    def env_setup_script() -> str:
        script_path = os.path.realpath(os.path.join("vomix", "env_setup.sh"))
        log.info(f"Running script: [cyan]{script_path}[/cyan]")

        cmd = ["bash", script_path]
        try:
            result = subprocess.run(
                cmd,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            return result.stdout
        except subprocess.CalledProcessError as e:
            log.error(f"Error executing setup script: [red]{e.stderr}[/red]")
            return f"Error: {e.stderr}"

    def createScript(self, module, module_obj, snakemake_obj):
        cwd = os.path.abspath(os.getcwd())

        script = f'snakemake --config module="{module}" '

        for attr, value in module_obj.__dict__.items():
            attr = attr.replace("_", "-")
            if str(value) == "True" or str(value) == "False":
                script += f"{attr}={value} "
            elif value is not None and attr != "custom-config" and attr != "name":
                if attr in ["samplelist", "datadir", "outdir", "fasta", "fastadir"]:
                    script += f'{attr}="{os.path.join(cwd, str(value))}" '
                else:
                    script += f'{attr}="{value}" '

        for attr, value in snakemake_obj.__dict__.items():
            if value is not None and attr != "snakemake_args":
                if isinstance(value, (list, tuple)) and len(value) == 0:
                    continue
                attr = attr.replace("_", "-")
                if isinstance(value, (list, tuple)):
                    val_str = " ".join(str(v) for v in value)
                    script += f"--{attr} {val_str} "
                elif str(value) == "True":
                    script += f"--{attr} "
                elif str(value) != "False":
                    if " " in str(value):
                        script += f'--{attr} "{value}" '
                    else:
                        script += f"--{attr} {value} "
            if attr == "snakemake_args" and value is not None and value != "":
                script += f"{value} "

        return script

    def _get_working_path(self):
        """Helper method to robustly find the working path by looking for the vomix parent."""
        current_vomix_dir = os.path.dirname(os.path.abspath(__file__))
        search_dir = current_vomix_dir
        currentWorkingPath = None

        # Walk up the directory tree one folder at a time
        while True:
            parent_dir, current_folder = os.path.split(search_dir)

            if current_folder == "vomix":
                currentWorkingPath = parent_dir
                break

            # If we reach the root directory without finding "vomix", stop looking
            if parent_dir == search_dir:
                break

            # Move up one level for the next iteration
            search_dir = parent_dir

        if not currentWorkingPath:
            log.error("Could not determine the path to the current working directory.")
            raise FileNotFoundError(
                "Could not determine the path to the current working directory."
            )

        return currentWorkingPath

    def createFoldersAndUpdateConfig(self, module_obj):
        currentWorkingPath = self._get_working_path()

        # Based on your previous logic replacing "/vomix" with "/config/config.yml"
        fullConfigPath = os.path.join(currentWorkingPath, "config", "config.yml")
        # Fallback if config is inside vomix instead of adjacent to it
        if not os.path.exists(fullConfigPath):
            fullConfigPath = os.path.join(
                currentWorkingPath, "vomix", "config", "config.yml"
            )

        log.info(f"Using Template config: [cyan]{fullConfigPath}[/cyan]")

        # get workdir
        workdir = module_obj.workdir
        if workdir is None:
            # get workdir from custom config
            if module_obj.custom_config is not None:
                with open(module_obj.custom_config) as f:
                    list_doc = yaml.safe_load(f)
                    workdir = list_doc["workdir"]
            else:
                # get workdir from template config.yml
                with open(fullConfigPath) as f:
                    list_doc = yaml.safe_load(f)
                    workdir = list_doc["workdir"]

        log.info(f"Working Directory: [green]{workdir}[/green]")

        outdir = os.path.join(workdir, module_obj.outdir)

        if module_obj.datadir is not None:
            datadir = os.path.join(workdir, module_obj.datadir)
            os.makedirs(datadir, exist_ok=True)

        if module_obj.fastadir is not None:
            fastadir = os.path.join(workdir, module_obj.fastadir)
            os.makedirs(fastadir, exist_ok=True)

        if not (
            os.path.exists(outdir) and os.path.exists(os.path.join(outdir, ".vomix"))
        ):
            os.makedirs(os.path.join(outdir, ".vomix"), exist_ok=True)

        now = datetime.datetime.now()
        latest_run = now.strftime("%Y%m%d_%H%M%S")
        outdir_folder = os.path.join(outdir, ".vomix", "log", f"vomix{latest_run}")

        os.makedirs(outdir_folder, exist_ok=True)

        # if custom config is specified
        if module_obj.custom_config is not None:
            log.info(f"Using custom config: [cyan]{module_obj.custom_config}[/cyan]")
            log.warning(
                "[yellow]REMINDER - Any command line flags specified will override those options in your custom config.[/yellow]"
            )
            shutil.copy(os.path.realpath(module_obj.custom_config), outdir_folder)

            # Use os.path.basename so it doesn't crash if custom_config was a nested path (e.g., 'configs/my_config.yml')
            copied_config_name = os.path.basename(module_obj.custom_config)
            os.rename(
                os.path.join(outdir_folder, copied_config_name),
                os.path.join(outdir_folder, "config.yml"),
            )
        else:
            shutil.copy(fullConfigPath, outdir_folder)

        # edit new config with user options + latest_run
        config_file_path = os.path.join(outdir_folder, "config.yml")
        with open(config_file_path) as f:
            list_doc = yaml.safe_load(f)
            list_doc["latest-run"] = latest_run

            for module, value in module_obj.__dict__.items():
                if value is not None:
                    module_formatted = module.replace("_", "-")
                    list_doc[module_formatted] = value
                    log.debug(f"Setting config option: {module_formatted} = {value}")

        with open(config_file_path, "w") as f:
            yaml.dump(list_doc, f)

        return outdir_folder

    def run_module(self, module, module_obj, snakemake_obj):
        outdir_folder = self.createFoldersAndUpdateConfig(module_obj)

        # create the script to run the module
        script_path = os.path.realpath(os.path.join(outdir_folder, "snakemake.sh"))
        script = self.createScript(module, module_obj, snakemake_obj)

        # save the script
        with open(script_path, "w") as f:
            f.write(script)

        log.info(f"Running Script: [cyan]{script_path}[/cyan]")
        cmd = ["bash", script_path]

        # Use the newly abstracted helper method instead of the old regex block
        currentWorkingPath = self._get_working_path()

        log.info(f"Working Path: [green]{currentWorkingPath}[/green]")
        log.info("Delegating execution to Snakemake backend...\n")

        try:
            # Create a pseudo-terminal pair to preserve Snakemake's color output
            master_fd, slave_fd = pty.openpty()

            with Popen(
                cmd,
                stdout=slave_fd,
                stderr=slave_fd,
                close_fds=True,
                cwd=currentWorkingPath,
            ) as p:
                # Close the slave fd in the parent process so EOF reads correctly
                os.close(slave_fd)

                # Read line by line and guarantee a trailing newline
                with os.fdopen(
                    master_fd, "r", encoding="utf-8", errors="replace"
                ) as master:
                    last_line = ""
                    for line in master:
                        print(line, end="", flush=True)
                        last_line = line

                    # If the stream ended without a newline, force one so the next log starts clean
                    if last_line and not last_line.endswith("\n"):
                        print()

            if p.returncode != 0:
                raise CalledProcessError(p.returncode, p.args)

        except (OSError, subprocess.CalledProcessError) as e:
            if isinstance(e, subprocess.CalledProcessError):
                log.error(f"Process failed with error code [red]{e.returncode}[/red]")
                return f"Error: {e.stderr}"
