#!/usr/bin/env python3
"""Check that a GTF is safe to hand to STAR / gffread / Salmon / featureCounts.

Every check here is something that has broken a real run:
  * a transcript_id used on two loci or by two genes (Salmon/tximport collapse
    or reject it);
  * a gene_id spread over chromosomes or strands (featureCounts/tximport);
  * transcripts without exons, overlapping exons, strand ".", start > end
    (gffread/STAR reject or silently mangle them);
  * transcript/gene lines that do not span their exons/transcripts.

Writes a JSON report and exits 1 if any check failed. With --transcripts-fa
(gffread -w output) it also checks every transcript was extracted.
"""
import argparse
import json
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gtf_utils import get_attr, iter_gtf, open_read

MAX_EXAMPLES = 10


class Problems:
    def __init__(self):
        self.items = defaultdict(list)

    def add(self, check, example):
        self.items[check].append(example)

    def summary(self):
        return {k: {"count": len(v), "examples": v[:MAX_EXAMPLES]} for k, v in self.items.items()}


def validate(path):
    """Return (stats, Problems) for the GTF at *path*."""
    bad = Problems()
    tx = {}  # tid -> dict(gene, chrom, strand, exons, line_span)
    gene_lines = defaultdict(list)
    gene_loci = defaultdict(set)

    for p in iter_gtf(path):
        chrom, feat, s, e, strand, attrs = p[0], p[2], p[3], p[4], p[6], p[8]
        if s > e:
            bad.add("start_after_end", f"{chrom}:{s}-{e} {feat}")
        if strand not in ("+", "-"):
            bad.add("no_strand", f"{chrom}:{s}-{e} {feat}")
        gid = get_attr(attrs, "gene_id")
        if gid is None:
            bad.add("missing_gene_id", f"{chrom}:{s}-{e} {feat}")
            continue
        if feat == "gene":
            gene_lines[gid].append((chrom, strand, s, e))
            continue
        tid = get_attr(attrs, "transcript_id")
        if tid is None:
            bad.add("missing_transcript_id", f"{gid} {feat}")
            continue
        t = tx.setdefault(
            tid, {"gene": gid, "chrom": chrom, "strand": strand, "exons": [], "span": None}
        )
        if t["gene"] != gid:
            bad.add("transcript_in_two_genes", f"{tid}: {t['gene']} / {gid}")
        if (t["chrom"], t["strand"]) != (chrom, strand):
            bad.add("transcript_on_two_loci", tid)
        gene_loci[gid].add((chrom, strand))
        if feat == "exon":
            t["exons"].append((s, e))
        elif feat == "transcript":
            if t["span"] is not None:
                bad.add("duplicate_transcript_line", tid)
            t["span"] = (s, e)

    gene_tx_span = {}
    for tid, t in tx.items():
        if not t["exons"]:
            bad.add("transcript_without_exons", tid)
            continue
        ex = sorted(t["exons"])
        if any(ex[i][1] >= ex[i + 1][0] for i in range(len(ex) - 1)):
            bad.add("overlapping_exons", tid)
        lo, hi = ex[0][0], max(x[1] for x in ex)
        if t["span"] is not None and t["span"] != (lo, hi):
            bad.add("transcript_span_mismatch", f"{tid}: line {t['span']} exons {(lo, hi)}")
        g = gene_tx_span.setdefault(t["gene"], [lo, hi])
        g[0], g[1] = min(g[0], lo), max(g[1], hi)

    for gid, loci in gene_loci.items():
        if len(loci) > 1:
            bad.add("gene_on_two_loci", f"{gid}: {sorted(loci)}")
    for gid, lines in gene_lines.items():
        if len(lines) > 1:
            bad.add("duplicate_gene_line", gid)
        span = gene_tx_span.get(gid)
        if span is None:
            bad.add("gene_line_without_transcripts", gid)
        elif (lines[0][2], lines[0][3]) != tuple(span):
            bad.add("gene_span_mismatch", f"{gid}: line {lines[0][2:]} transcripts {tuple(span)}")

    stats = {
        "transcripts": len(tx),
        "genes": len(gene_loci),
        "gene_lines": len(gene_lines),
    }
    return stats, bad


def count_fasta(path):
    with open_read(path) as fh:
        return sum(1 for line in fh if line.startswith(">"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--transcripts-fa", default=None)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    stats, bad = validate(a.gtf)
    if a.transcripts_fa:
        n = count_fasta(a.transcripts_fa)
        stats["transcripts_extracted"] = n
        if n != stats["transcripts"]:
            bad.add(
                "gffread_extraction",
                f"{n} sequences extracted for {stats['transcripts']} transcripts",
            )
    report = {"gtf": a.gtf, "ok": not bad.items, "stats": stats, "problems": bad.summary()}
    with open(a.out, "w") as fh:
        json.dump(report, fh, indent=2)
    if bad.items:
        print(f"GTF validation FAILED for {a.gtf}:", file=sys.stderr)
        for check, info in report["problems"].items():
            print(f"  {check}: {info['count']}  e.g. {info['examples'][:3]}", file=sys.stderr)
        sys.exit(1)
    print(f"GTF validation passed: {stats}", file=sys.stderr)


if __name__ == "__main__":
    main()
