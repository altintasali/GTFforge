# Sample sheet: load, validate, groups, per-sample inputs.
#
# columns: sample, group, gtf | bam, strandedness (optional)

if not os.path.exists(config["samples"]):
    raise WorkflowError(
        f"samples: '{config['samples']}' does not exist.\n"
        f"Relative paths resolve against the directory you run snakemake "
        f"from ({os.getcwd()}), not workflow/ -- cd into the repo/run directory."
    )


def _nonempty(value):
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    text = str(value).strip()
    return text or None


_raw = pd.read_csv(config["samples"], dtype=str, comment="#")
for _col in ("gtf", "bam", "strandedness"):
    if _col not in _raw.columns:
        _raw[_col] = None
# jsonschema's "null" only matches None, not NaN.
_raw = _raw.astype(object).where(pd.notnull(_raw), None)
validate(_raw, schema="../../schemas/samples.schema.yaml")

_rows = []
for _, row in _raw.iterrows():
    sample, group = _nonempty(row["sample"]), _nonempty(row["group"])
    gtf, bam = _nonempty(row["gtf"]), _nonempty(row["bam"])
    if sample is None or group is None:
        raise WorkflowError(f"{config['samples']}: every row needs a sample and a group.")
    if (gtf is None) == (bam is None):
        raise WorkflowError(
            f"sample '{sample}': give exactly one of 'gtf' (an existing assembly) "
            f"or 'bam' (assembled by GTFforge) in {config['samples']}."
        )
    _rows.append(
        {
            "sample": sample,
            "group": group,
            "gtf": gtf,
            "bam": bam,
            "strandedness": _nonempty(row["strandedness"]) or "unstranded",
        }
    )

samples = pd.DataFrame(_rows)
_dups = sorted(samples["sample"][samples["sample"].duplicated()].unique())
if _dups:
    raise WorkflowError(f"duplicated sample names in {config['samples']}: {', '.join(_dups)}")
samples = samples.set_index("sample", drop=False).sort_index()

SAMPLES = list(samples["sample"])
BAM_SAMPLES = list(samples.loc[samples["bam"].notna(), "sample"])
# Sorted with Python's ordering (not the shell's locale): the support_map rule
# and finalize.py both rely on this exact order.
GROUPS = sorted(samples["group"].unique())
GROUP_SAMPLES = {g: sorted(samples.loc[samples["group"] == g, "sample"]) for g in GROUPS}

_small = {g: len(s) for g, s in GROUP_SAMPLES.items() if len(s) < 2}
if _small:
    raise WorkflowError(
        "every group needs at least 2 samples -- replicate support is the "
        "point of GTFforge. Groups with one sample: "
        + ", ".join(f"{g} ({n})" for g, n in _small.items())
        + "\nMerge them into a larger group or drop them from the sample sheet."
    )


def min_support(n):
    s = config["support"]
    return max(s["min_samples"], math.ceil(n * s["min_fraction"] - 1e-9))


_impossible = {g: len(s) for g, s in GROUP_SAMPLES.items() if min_support(len(s)) > len(s)}
if _impossible:
    logger.warning(
        "support.min_samples is larger than these groups, so none of their "
        "transcripts can pass the support filter: "
        + ", ".join(f"{g} (n={n}, needs {min_support(n)})" for g, n in _impossible.items())
    )


def group_single_exon(group):
    """Resolve support.single_exon for a group: 'drop' or 'keep_stranded'."""
    mode = config["support"]["single_exon"]
    if mode != "auto":
        return mode
    stranded = all(samples.loc[s, "strandedness"] != "unstranded" for s in GROUP_SAMPLES[group])
    return "keep_stranded" if stranded else "drop"


def sample_gtf(sample):
    """The assembly GTF of a sample: as given, or GTFforge's own for BAM rows."""
    gtf = samples.loc[sample, "gtf"]
    return gtf if gtf else f"results/assembly/{sample}.gtf"


_missing = [
    f"{s}: {p}"
    for s in SAMPLES
    for p in (samples.loc[s, "gtf"], samples.loc[s, "bam"])
    if p and not os.path.exists(p)
]
if _missing:
    logger.warning(
        f"{len(_missing)} input file(s) in the sample sheet do not exist yet:\n  "
        + "\n  ".join(_missing[:20])
        + ("\n  ..." if len(_missing) > 20 else "")
    )
