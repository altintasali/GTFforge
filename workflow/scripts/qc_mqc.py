#!/usr/bin/env python3
"""Turn GTFforge's JSON summaries into MultiQC custom-content sections.

Writes one *_mqc.json per section into --outdir:
  overview              headline numbers (table)
  support_histogram     per group: how many replicates assembled each chain
  support_kept          per group: chains kept vs dropped, by reason
  novel_classes         kept novel transcripts by gffcompare class code
  novel_by_group        kept novel transcripts supported by each group
  dropped               every filter and how much it removed
  accuracy              gffcompare sensitivity/precision of the merge vs reference
  te_families           TE families at novel transcript starts (with ref.te_gtf)
"""
import argparse
import json
import os
import re

CLASS_LABELS = {
    "=": "= identical intron chain",
    "j": "j new isoform (shares a junction)",
    "k": "k contains a reference transcript",
    "x": "x antisense",
    "u": "u intergenic",
    "i": "i inside an intron",
    "o": "o other same-strand overlap",
    "y": "y contains a reference in its intron",
}


def section(id_, name, description, plot_type, data, pconfig=None, **extra):
    out = {
        "id": f"gtfforge_{id_}",
        "parent_id": "gtfforge",
        "parent_name": "GTFforge",
        "parent_description": "Replicate-supported, reference-anchored custom transcriptome.",
        "section_name": name,
        "description": description,
        "plot_type": plot_type,
        "data": data,
    }
    if pconfig:
        out["pconfig"] = {"id": f"gtfforge_{id_}_plot", "title": f"GTFforge: {name}", **pconfig}
    out.update(extra)
    return out


def parse_gffcompare_stats(path):
    """{level: {"Sensitivity": x, "Precision": y}} from a gffcompare .stats file."""
    out = {}
    pat = re.compile(r"^\s*([A-Za-z ]+level):\s+([\d.]+)\s+\|\s+([\d.]+)\s+\|")
    with open(path) as fh:
        for line in fh:
            m = pat.match(line)
            if m:
                out[m.group(1).strip()] = {
                    "Sensitivity": float(m.group(2)),
                    "Precision": float(m.group(3)),
                }
    return out


def build(support_stats, finalize, validation, accuracy, te=None):
    """Return {filename: section dict}."""
    secs = {}
    v = validation["stats"]
    secs["overview"] = section(
        "overview",
        "Overview",
        "What went into the final GTF. Novel isoforms keep their reference "
        "gene_id, so gene-level counts of known genes stay comparable with a "
        "reference-only run.",
        "table",
        {
            "GTFforge": {
                "Reference transcripts": finalize["reference_transcripts_kept"],
                "Novel transcripts": finalize["novel_transcripts_kept"],
                "Novel isoforms of ref. genes": finalize["novel_isoforms_of_reference_genes"],
                "Novel genes": finalize["novel_genes"],
                "Transcripts in GTF": v["transcripts"],
                "Genes in GTF": v["genes"],
                "Validation": "passed" if validation["ok"] else "FAILED",
            }
        },
        {"col1_header": "Run"},
    )
    secs["support_histogram"] = section(
        "support_histogram",
        "Replicate support",
        "For each group, the number of unique intron chains assembled in exactly "
        "N of its samples. A long tail at N=1 is normal: most single-sample "
        "chains are noise, which is what the support filter removes.",
        "linegraph",
        {
            s["group"]: {int(k): n for k, n in s["support_histogram"].items()}
            for s in support_stats
        },
        {"xlab": "Samples in the group that assembled the chain", "ylab": "Intron chains",
         "logswitch": True},
    )
    secs["support_kept"] = section(
        "support_kept",
        "Support filter",
        "Intron chains kept per group, and why the rest were dropped. "
        "Threshold: max(support.min_samples, ceil(n * support.min_fraction)).",
        "bargraph",
        {
            f"{s['group']} (n={s['n_samples']}, >={s['min_support']})": {
                "kept": s["n_kept"],
                **{k.replace("_", " "): n for k, n in s["dropped"].items()},
            }
            for s in support_stats
        },
        {"ylab": "Intron chains", "cpswitch_c_active": True},
    )
    secs["novel_classes"] = section(
        "novel_classes",
        "Novel transcripts by class",
        "gffcompare class code of each novel transcript relative to the reference.",
        "bargraph",
        {"novel transcripts": {CLASS_LABELS.get(c, c): n for c, n in finalize["kept_by_class"].items()}},
        {"ylab": "Transcripts"},
    )
    secs["novel_by_group"] = section(
        "novel_by_group",
        "Novel transcripts by group",
        "Novel transcripts whose exact intron chain was supported in each group "
        "(a transcript can be supported by several groups).",
        "bargraph",
        {g: {"novel transcripts": n} for g, n in finalize["kept_by_group"].items()},
        {"ylab": "Transcripts"},
    )
    secs["dropped"] = section(
        "dropped",
        "Filters after merging",
        "Everything finalize removed, and why: reference transcripts reused on "
        "several loci or on excluded contigs, and novel transcripts with an "
        "excluded class code, no strand, an excluded contig, or no exact "
        "group-level support.",
        "table",
        {k.replace("_", " "): {"count": n} for k, n in finalize["dropped"].items()},
        {"col1_header": "Filter"},
    )
    if accuracy:
        secs["accuracy"] = section(
            "accuracy",
            "Merged transcriptome vs reference",
            "gffcompare of the merged transcriptome against the reference. "
            "Sensitivity near 100 means the reference was kept; lower precision "
            "reflects the novel transcripts added.",
            "table",
            accuracy,
            {"col1_header": "Level"},
        )
    if te is not None:
        secs["te_families"] = section(
            "te_families",
            "TEs at novel transcript starts",
            f"{te['with_te_at_tss']} of {te['novel_transcripts']} novel transcripts "
            "start inside a transposable element. Families shown by number of "
            "transcripts; see te_at_tss in results/novel_transcripts.tsv for the "
            "subfamily and orientation.",
            "bargraph",
            {"novel transcripts": dict(list(te["te_families_at_tss"].items())[:25])},
            {"ylab": "Transcripts"},
        )
    return secs


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--support-stats", nargs="+", required=True)
    ap.add_argument("--finalize-stats", required=True)
    ap.add_argument("--validation", required=True)
    ap.add_argument("--classify-stats", required=True)
    ap.add_argument("--te-stats", default=None)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args(argv)

    def load(p):
        with open(p) as fh:
            return json.load(fh)

    secs = build(
        [load(p) for p in a.support_stats],
        load(a.finalize_stats),
        load(a.validation),
        parse_gffcompare_stats(a.classify_stats),
        load(a.te_stats) if a.te_stats else None,
    )
    os.makedirs(a.outdir, exist_ok=True)
    for name, sec in secs.items():
        with open(os.path.join(a.outdir, f"{name}_mqc.json"), "w") as fh:
            json.dump(sec, fh, indent=2)


if __name__ == "__main__":
    main()
