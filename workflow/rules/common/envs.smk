# Per-tool conda environments generated from config["versions"], plus the
# pipeline's own version/commit for the report.
#
# Absolute path: rules are included from workflow/Snakefile and Snakemake
# resolves relative paths in an included file against that file's directory,
# while the env files are written relative to the process CWD.

import yaml

GENERATED_ENV_DIR = os.path.abspath("workflow/envs/generated")
os.makedirs(GENERATED_ENV_DIR, exist_ok=True)


def _write_env(name, dependencies):
    path = f"{GENERATED_ENV_DIR}/{name}.yaml"
    with open(path, "w") as fh:
        yaml.safe_dump(
            {"channels": ["conda-forge", "bioconda"], "dependencies": list(dependencies)},
            fh,
            sort_keys=False,
        )
    return path


STRINGTIE_ENV = _write_env("stringtie", [f"stringtie={V['stringtie']}"])
GFFCOMPARE_ENV = _write_env("gffcompare", [f"gffcompare={V['gffcompare']}"])
GFFREAD_ENV = _write_env("gffread", [f"gffread={V['gffread']}"])
MULTIQC_ENV = _write_env("multiqc", [f"multiqc={V['multiqc']}"])
PYTHON_ENV = _write_env("python", ["python>=3.9"])

# realpath: run directories link workflow/ to a shared checkout.
_REPO_ROOT = os.path.dirname(os.path.realpath(workflow.basedir))
with open(os.path.join(_REPO_ROOT, "VERSION")) as _fh:
    PIPELINE_VERSION = _fh.read().strip()


def _git_state():
    """Commit (and dirty flag) of the checkout the workflow runs from."""
    try:
        commit = subprocess.run(
            ["git", "-C", _REPO_ROOT, "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "-C", _REPO_ROOT, "status", "--porcelain", "--untracked-files=no"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
        return commit + ("-dirty" if dirty else "")
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


PIPELINE_COMMIT = _git_state()
SCRIPTS_DIR = os.path.join(_REPO_ROOT, "workflow", "scripts")
