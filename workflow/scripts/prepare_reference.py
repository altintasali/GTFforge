#!/usr/bin/env python3
"""Clean the reference annotation once, before anything is compared to it.

Drops transcripts that every later step would misread:
  * a transcript_id used for several copies (different chromosome/strand,
    repeated transcript or first-exon lines, overlapping exons) -- UCSC
    refGene GTFs do this for duplicated genes. Left in, one such id spans
    every copy and the region between, and gffcompare then calls genuinely
    intergenic transcripts "intronic" to it;
  * transcripts on contigs outside --keep-contigs-regex.
Gene lines are kept for genes with at least one surviving transcript; comment
lines and gzip are removed so gffcompare/StringTie/gffread can all read it.
"""
import argparse
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from finalize import Reference, plan_reference
from gtf_utils import get_attr, iter_gtf, open_write


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--keep-contigs-regex", default="")
    ap.add_argument("--out-gtf", required=True)
    ap.add_argument("--out-stats", required=True)
    a = ap.parse_args(argv)

    contig_re = re.compile(a.keep_contigs_regex) if a.keep_contigs_regex else None

    def keep_contig(c):
        return contig_re is None or contig_re.search(c) is not None

    ref = Reference()
    for p in iter_gtf(a.gtf):
        ref.add(p)
    kept, _, drops, _ = plan_reference(ref, keep_contig)
    kept_genes = {ref.tx_gene[t] for t in kept}

    with open_write(a.out_gtf) as oh:
        for p in iter_gtf(a.gtf):
            if p[2] == "gene":
                if get_attr(p[8], "gene_id") not in kept_genes:
                    continue
            elif get_attr(p[8], "transcript_id") not in kept:
                continue
            oh.write("\t".join(map(str, p)) + "\n")

    stats = {
        "reference_transcripts_in": len(ref.tx_gene),
        "reference_transcripts_kept": len(kept),
        "dropped": dict(drops),
    }
    with open(a.out_stats, "w") as fh:
        json.dump(stats, fh, indent=2)
    print(json.dumps(stats), file=sys.stderr)


if __name__ == "__main__":
    main()
