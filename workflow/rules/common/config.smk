# Imports, config validation, per-rule resource lookup.

import math
import os
import subprocess

import pandas as pd
from snakemake.exceptions import WorkflowError
from snakemake.logging import logger
from snakemake.utils import validate

validate(config, schema="../../schemas/config.schema.yaml")

V = config["versions"]

# Per-rule threads/mem_mb/runtime from workflow/default-config/resources.yaml,
# overridable in input/resources.yaml. Unlisted rules get small defaults.
RESOURCES = config.get("resources", {})
_RESOURCE_DEFAULTS = {"threads": 1, "mem_mb": 4000, "runtime": 60}


def get_resources(rule_name):
    """Return {threads, mem_mb, runtime} for a rule name."""
    res = {**_RESOURCE_DEFAULTS, **RESOURCES.get(rule_name, {})}
    res.pop("mem_per_sample", None)
    return res


def get_scaled_mem_mb(rule_name, n_samples):
    """Base mem_mb plus mem_per_sample for each of *n_samples* samples."""
    per_sample = RESOURCES.get(rule_name, {}).get("mem_per_sample", 0)
    return get_resources(rule_name)["mem_mb"] + per_sample * n_samples
