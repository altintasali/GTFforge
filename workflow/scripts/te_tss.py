#!/usr/bin/env python3
"""Annotate novel transcripts with the transposable elements at their start.

For every novel transcript in --novel-tsv, finds the TE insertions (from a
TEtranscripts-style rmsk GTF: gene_id = subfamily, family_id, class_id)
overlapping
  * its transcription start site (1 bp, strand-aware) -> te_at_tss
  * its first exon                                     -> te_first_exon
reported as subfamily:family:class(sense|anti), orientation relative to the
transcript. Adds those two columns to the table and writes a family summary
for the report.
"""
import argparse
import bisect
import json
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gtf_utils import get_attr, iter_gtf, open_read, open_write


class IntervalIndex:
    """Per-chromosome sorted intervals; overlap queries via bisect on starts.

    Every interval starting within max_len of the query's start is a
    candidate, so a query costs O(log n + candidates) without an interval tree.
    """

    def __init__(self):
        self._by_chrom = defaultdict(list)
        self._starts = {}
        self._max_len = defaultdict(int)

    def add(self, chrom, start, end, payload):
        self._by_chrom[chrom].append((start, end, payload))
        self._max_len[chrom] = max(self._max_len[chrom], end - start + 1)

    def build(self):
        for chrom, ivs in self._by_chrom.items():
            ivs.sort(key=lambda x: (x[0], x[1]))
            self._starts[chrom] = [x[0] for x in ivs]
        return self

    def overlapping(self, chrom, start, end):
        ivs = self._by_chrom.get(chrom)
        if not ivs:
            return []
        starts = self._starts[chrom]
        lo = bisect.bisect_left(starts, start - self._max_len[chrom] + 1)
        hi = bisect.bisect_right(starts, end)
        return [iv[2] for iv in ivs[lo:hi] if iv[1] >= start]


def load_te(path):
    idx = IntervalIndex()
    for p in iter_gtf(path):
        if p[2] != "exon":
            continue
        a = p[8]
        label = f"{get_attr(a, 'gene_id')}:{get_attr(a, 'family_id')}:{get_attr(a, 'class_id')}"
        idx.add(p[0], p[3], p[4], (label, p[6]))
    return idx.build()


def first_exons(gtf, novel_ids):
    """{tid: (chrom, strand, first_exon_start, first_exon_end)} for novel transcripts."""
    exons = defaultdict(list)
    for p in iter_gtf(gtf):
        if p[2] != "exon":
            continue
        tid = get_attr(p[8], "transcript_id")
        if tid in novel_ids:
            exons[tid].append((p[0], p[6], p[3], p[4]))
    out = {}
    for tid, ex in exons.items():
        chrom, strand = ex[0][0], ex[0][1]
        ex.sort(key=lambda x: x[2])
        s, e = (ex[0][2], ex[0][3]) if strand == "+" else (ex[-1][2], ex[-1][3])
        out[tid] = (chrom, strand, s, e)
    return out


def describe(hits, strand):
    labels = sorted({f"{lab}({'sense' if s == strand else 'anti'})" for lab, s in hits})
    return ",".join(labels) or "-"


def annotate(idx, fe):
    """{tid: (te_at_tss, te_first_exon)}."""
    out = {}
    for tid, (chrom, strand, s, e) in fe.items():
        tss = s if strand == "+" else e
        out[tid] = (
            describe(idx.overlapping(chrom, tss, tss), strand),
            describe(idx.overlapping(chrom, s, e), strand),
        )
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--te-gtf", required=True)
    ap.add_argument("--novel-tsv", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-stats", required=True)
    a = ap.parse_args(argv)

    with open_read(a.novel_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    ids = {r[0] for r in rows}
    ann = annotate(load_te(a.te_gtf), first_exons(a.gtf, ids))

    families = Counter()
    with open_write(a.out_tsv) as oh:
        oh.write("\t".join([*header, "te_at_tss", "te_first_exon"]) + "\n")
        for r in rows:
            tss, fe = ann.get(r[0], ("-", "-"))
            oh.write("\t".join([*r, tss, fe]) + "\n")
            if tss != "-":
                for fam in {x.split(":")[1] for x in tss.split(",")}:
                    families[fam] += 1
    stats = {
        "novel_transcripts": len(rows),
        "with_te_at_tss": sum(1 for r in rows if ann.get(r[0], ("-",))[0] != "-"),
        "te_families_at_tss": dict(families.most_common()),
    }
    with open(a.out_stats, "w") as fh:
        json.dump(stats, fh, indent=2)
    print(f"novel transcripts with a TE at the TSS: {stats['with_te_at_tss']} / {len(rows)}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
