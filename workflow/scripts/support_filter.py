#!/usr/bin/env python3
"""Keep a group's consensus transcripts that enough replicates agree on.

Input is one gffcompare run over every sample GTF of a group. Its .tracking
file has one row per unique intron chain (TCONS id) and one column per sample
("-" where that sample did not assemble the chain). A chain is kept when it
was seen in at least

    min_support = max(min_samples, ceil(n_samples * min_fraction))

samples. Single-exon chains (only present when gffcompare ran without -M)
are kept only if they carry a strand; unstranded single-exon transfrags
cannot be placed on a strand and would compete with genes in quantification.

Writes the kept transcripts (from <prefix>.combined.gtf, each tagged with a
`support_group` attribute) and a JSON summary for the report.
"""
import argparse
import json
import math
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gtf_utils import get_attr, open_read, open_write


def min_support(n_samples, min_samples, min_fraction):
    """Replicates needed to keep a chain. The epsilon keeps e.g. 12 * 0.25
    from rounding up to 4 through floating-point error."""
    return max(min_samples, math.ceil(n_samples * min_fraction - 1e-9))


def read_tracking(path):
    """Return ({tcons_id: n_samples_present}, n_samples)."""
    support, n = {}, None
    with open_read(path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            present = cols[4:]
            if n is None:
                n = len(present)
            support[cols[0]] = sum(q != "-" for q in present)
    return support, (n or 0)


def read_structure(path):
    """Return {transcript_id: (strand, n_exons)} from a gffcompare combined GTF."""
    strand, exons = {}, Counter()
    with open_read(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.split("\t")
            if len(p) < 9:
                continue
            tid = get_attr(p[8], "transcript_id")
            strand[tid] = p[6]
            if p[2] == "exon":
                exons[tid] += 1
    return {t: (strand[t], exons[t]) for t in strand}


def select(support, structure, needed):
    """Return (kept ids, Counter of drop reasons)."""
    kept, dropped = set(), Counter()
    for tid, n in support.items():
        if n < needed:
            dropped["below_min_support"] += 1
            continue
        strand, n_exons = structure.get(tid, (".", 0))
        if n_exons <= 1 and strand not in "+-":
            dropped["single_exon_unstranded"] += 1
            continue
        kept.add(tid)
    return kept, dropped


def write_kept(combined, kept, group, out):
    with open_read(combined) as fh, open_write(out) as oh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.split("\t")
            if len(p) < 9 or get_attr(p[8], "transcript_id") not in kept:
                continue
            oh.write(line.rstrip("\n").rstrip() + f' support_group "{group}";\n')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tracking", required=True)
    ap.add_argument("--combined", required=True)
    ap.add_argument("--group", required=True)
    ap.add_argument("--min-samples", type=int, required=True)
    ap.add_argument("--min-fraction", type=float, required=True)
    ap.add_argument("--out-gtf", required=True)
    ap.add_argument("--out-stats", required=True)
    a = ap.parse_args(argv)

    support, n = read_tracking(a.tracking)
    needed = min_support(n, a.min_samples, a.min_fraction)
    kept, dropped = select(support, read_structure(a.combined), needed)
    write_kept(a.combined, kept, a.group, a.out_gtf)

    hist = defaultdict(int)
    for v in support.values():
        hist[v] += 1
    stats = {
        "group": a.group,
        "n_samples": n,
        "min_support": needed,
        "n_chains": len(support),
        "n_kept": len(kept),
        "dropped": dict(dropped),
        "support_histogram": {str(k): hist[k] for k in sorted(hist)},
    }
    with open(a.out_stats, "w") as fh:
        json.dump(stats, fh, indent=2)
    print(
        f"[{a.group}] n_samples={n} min_support={needed} "
        f"kept {len(kept)} / {len(support)} chains",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
