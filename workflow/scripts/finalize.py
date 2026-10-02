#!/usr/bin/env python3
"""Build the final GTF: the reference annotation plus filtered novel transcripts.

Inputs
  --ref         reference GTF (plain or .gz)
  --merged      stringtie --merge output (reference + novel transcripts)
  --classified  gffcompare -r <ref> <merged> annotated GTF (class codes)
  --support     gffcompare -r <merged> <group supported GTFs...> tracking file
  --groups      group names, in the order their GTFs were given to that run

Reference transcripts are copied verbatim, except:
  * transcripts on contigs outside --keep-contigs-regex are dropped;
  * a transcript_id that occurs at more than one locus (different
    chromosome/strand, or overlapping exons -- UCSC refGene GTFs do this for
    duplicated genes) is dropped entirely, as in Staubli et al.;
  * a gene_id spread over several chromosomes/strands is split into
    <gene_id>_<chrom><strand> so every gene has one locus;
  * gene lines are regenerated so they span every transcript of the gene,
    including novel isoforms (or dropped, with --gene-lines drop).

Novel transcripts (merged transcript_ids absent from the reference) are kept
when their class code is in --keep-classes, they carry a strand, sit on a
kept contig and (unless --no-require-support) have an exact intron-chain
match in at least one group. Classes in ATTACH_CLASSES share junctions with,
or contain, a reference transcript and join that reference gene; every other
class becomes, or joins, a novel gene named after its StringTie locus.
"""
import argparse
import json
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gtf_utils import format_attrs, get_attr, iter_gtf, open_read, open_write

# gffcompare classes that mean "a new isoform of this reference gene".
ATTACH_CLASSES = set("=jk")
GENE_TYPE_KEYS = ("gene_type", "gene_biotype")
TRANSCRIPT_TYPE_KEYS = ("transcript_type", "transcript_biotype")


# ----------------------------------------------------------------------------
# reference
# ----------------------------------------------------------------------------
class Reference:
    """Everything finalize needs to know about the reference, from one pass."""

    def __init__(self):
        self.tx_gene = {}  # transcript_id -> gene_id
        self.tx_loci = defaultdict(set)  # transcript_id -> {(chrom, strand)}
        self.tx_exons = defaultdict(list)  # transcript_id -> [(start, end)]
        self.tx_span = {}  # transcript_id -> [start, end]
        self.tx_repeats = Counter()  # transcript_id -> extra transcript/first-exon lines
        self.gene_line = {}  # gene_id -> original gene line columns
        self.gene_name = {}
        self.gene_type = {}
        self.keys_seen = Counter()

    def add(self, p):
        attrs = p[8]
        gid = get_attr(attrs, "gene_id")
        if p[2] == "gene":
            if gid is not None:
                self.gene_line[gid] = p
            self._gene_info(gid, attrs)
            return
        tid = get_attr(attrs, "transcript_id")
        if tid is None or gid is None:
            return
        self.tx_gene[tid] = gid
        self.tx_loci[tid].add((p[0], p[6]))
        span = self.tx_span.setdefault(tid, [p[3], p[4]])
        span[0], span[1] = min(span[0], p[3]), max(span[1], p[4])
        if p[2] == "exon":
            self.tx_exons[tid].append((p[3], p[4]))
            first = get_attr(attrs, "exon_number") == "1"
        else:
            first = p[2] == "transcript"
        if first:
            self.tx_repeats[(tid, p[2])] += 1
        self._gene_info(gid, attrs)

    def multi_locus(self, tid):
        """A transcript_id used for more than one copy of a transcript."""
        return (
            len(self.tx_loci[tid]) > 1
            or self.tx_repeats[(tid, "transcript")] > 1
            or self.tx_repeats[(tid, "exon")] > 1
            or has_overlapping_exons(self.tx_exons.get(tid, []))
        )

    def _gene_info(self, gid, attrs):
        if gid in self.gene_type:
            return
        for k in ("gene_name", *GENE_TYPE_KEYS, *TRANSCRIPT_TYPE_KEYS):
            if f"{k} " in attrs:
                self.keys_seen[k] += 1
        name = get_attr(attrs, "gene_name")
        if name is not None:
            self.gene_name[gid] = name
        gtype = next((v for k in GENE_TYPE_KEYS if (v := get_attr(attrs, k)) is not None), None)
        if gtype is not None:
            self.gene_type[gid] = gtype

    def key(self, candidates, default):
        """The attribute name this annotation uses (gene_type vs gene_biotype ...)."""
        seen = [k for k in candidates if self.keys_seen[k]]
        return max(seen, key=lambda k: self.keys_seen[k]) if seen else default


def has_overlapping_exons(exons):
    exons = sorted(exons)
    return any(exons[i][1] >= exons[i + 1][0] for i in range(len(exons) - 1))


def plan_reference(ref, keep_contig):
    """Decide which reference transcripts survive and which genes get split.

    Returns (kept_tx: set, gene_of: {tid: final gene_id}, drops: Counter,
    number of reference genes split).
    """
    drops = Counter()
    kept = set()
    for tid, loci in ref.tx_loci.items():
        if ref.multi_locus(tid):
            drops["ref_transcript_multi_locus"] += 1
            continue
        (chrom, _strand), = loci
        if not keep_contig(chrom):
            drops["ref_transcript_contig"] += 1
            continue
        kept.add(tid)

    gene_loci = defaultdict(set)
    for tid in kept:
        gene_loci[ref.tx_gene[tid]].add(next(iter(ref.tx_loci[tid])))
    gene_of = {}
    for tid in kept:
        gid = ref.tx_gene[tid]
        if len(gene_loci[gid]) > 1:
            chrom, strand = next(iter(ref.tx_loci[tid]))
            gene_of[tid] = f"{gid}_{chrom}{strand}"
        else:
            gene_of[tid] = gid
    n_split = sum(len(v) > 1 for v in gene_loci.values())
    return kept, gene_of, drops, n_split


# ----------------------------------------------------------------------------
# novel transcripts
# ----------------------------------------------------------------------------
def read_classes(path):
    """{transcript_id: (class_code, cmp_ref)} from a gffcompare annotated GTF."""
    out = {}
    for p in iter_gtf(path):
        if p[2] == "transcript":
            out[get_attr(p[8], "transcript_id")] = (
                get_attr(p[8], "class_code"),
                get_attr(p[8], "cmp_ref"),
            )
    return out


def read_support(path, groups):
    """{merged transcript_id: {groups with an exact (=) intron-chain match}}."""
    support = defaultdict(set)
    with open_read(path) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 5 or c[3] != "=" or "|" not in c[2]:
                continue
            mtid = c[2].split("|", 1)[1]
            for g, q in zip(groups, c[4:]):
                if q != "-":
                    support[mtid].add(g)
    return support


def read_novel(path, ref_tids):
    """Novel transcripts of the merged GTF: {tid: dict(chrom, strand, locus, exons)}."""
    tx = {}
    for p in iter_gtf(path):
        tid = get_attr(p[8], "transcript_id")
        if tid is None or tid in ref_tids:
            continue
        t = tx.setdefault(
            tid,
            {"chrom": p[0], "strand": p[6], "locus": get_attr(p[8], "gene_id"), "exons": []},
        )
        if p[2] == "exon":
            t["exons"].append((p[3], p[4]))
    for t in tx.values():
        t["exons"].sort()
    return {tid: t for tid, t in tx.items() if t["exons"]}


def select_novel(novel, classes, support, keep_classes, keep_contig, require_support):
    """Filter novel transcripts; return (kept {tid: t+code+cref}, drop Counter)."""
    kept, drops = {}, Counter()
    for tid, t in novel.items():
        code, cref = classes.get(tid, (None, None))
        if code not in keep_classes:
            drops[f"novel_class_{code}"] += 1
        elif t["strand"] not in ("+", "-"):
            drops["novel_no_strand"] += 1
        elif not keep_contig(t["chrom"]):
            drops["novel_contig"] += 1
        elif require_support and not support.get(tid):
            drops["novel_no_exact_group_support"] += 1
        else:
            kept[tid] = {**t, "code": code, "cref": cref or ""}
    return kept, drops


def assign_genes(kept, ref_gene_of, ref_tx_locus):
    """Give every kept novel transcript a gene_id; returns n novel genes split.

    ATTACH classes join the reference gene of cmp_ref when that transcript
    survived and sits on the same chromosome/strand. Everything else forms a
    novel gene from its StringTie locus id, split per chromosome/strand if
    StringTie's locus spans more than one (it never should, but a gene_id on
    two strands breaks every downstream counter).
    """
    for t in kept.values():
        cref = t["cref"]
        if (
            t["code"] in ATTACH_CLASSES
            and cref in ref_gene_of
            and ref_tx_locus[cref] == (t["chrom"], t["strand"])
        ):
            t["gid"], t["novel_gene"] = ref_gene_of[cref], False
        else:
            t["gid"], t["novel_gene"] = t["locus"], True

    loci = defaultdict(set)
    for t in kept.values():
        if t["novel_gene"]:
            loci[t["gid"]].add((t["chrom"], t["strand"]))
    split = {g for g, s in loci.items() if len(s) > 1}
    for t in kept.values():
        if t["novel_gene"] and t["gid"] in split:
            t["gid"] = f"{t['gid']}_{t['chrom']}{t['strand']}"
    return len(split)


# ----------------------------------------------------------------------------
# writing
# ----------------------------------------------------------------------------
def gene_spans(ref, ref_kept, ref_gene_of, novel_kept):
    spans = {}
    for tid in ref_kept:
        chrom, strand = next(iter(ref.tx_loci[tid]))
        s, e = ref.tx_span[tid]
        _extend(spans, ref_gene_of[tid], chrom, strand, s, e)
    for t in novel_kept.values():
        _extend(spans, t["gid"], t["chrom"], t["strand"], t["exons"][0][0], t["exons"][-1][1])
    return spans


def _extend(spans, gid, chrom, strand, s, e):
    if gid in spans:
        sp = spans[gid]
        sp[2], sp[3] = min(sp[2], s), max(sp[3], e)
    else:
        spans[gid] = [chrom, strand, s, e]


def gene_line(gid, span, ref, original_gid, keys, novel_attrs=None):
    chrom, strand, s, e = span
    orig = ref.gene_line.get(original_gid)
    if orig is not None:
        attrs = orig[8]
        if gid != original_gid:
            attrs = attrs.replace(f'gene_id "{original_gid}"', f'gene_id "{gid}"', 1)
        source = orig[1]
    elif novel_attrs is not None:
        attrs, source = novel_attrs, "GTFforge"
    else:
        attrs = format_attrs(
            [
                ("gene_id", gid),
                ("gene_name", ref.gene_name.get(original_gid, original_gid)),
                (keys["gene_type"], ref.gene_type.get(original_gid)),
            ]
        )
        source = "GTFforge"
    return "\t".join([chrom, source, "gene", str(s), str(e), ".", strand, ".", attrs]) + "\n"


def write_gtf(args, ref, ref_kept, ref_gene_of, novel_kept, groups_of, keys, out):
    spans = gene_spans(ref, ref_kept, ref_gene_of, novel_kept)
    regenerate = args.gene_lines == "regenerate"
    written_genes = set()
    gid_re = re.compile(r'gene_id "[^"]*"')

    with open_write(out) as oh:
        # reference, in its own order; each gene line just before its first feature
        for p in iter_gtf(args.ref):
            if p[2] == "gene":
                continue
            tid = get_attr(p[8], "transcript_id")
            if tid not in ref_kept:
                continue
            gid = ref_gene_of[tid]
            if regenerate and gid not in written_genes:
                oh.write(gene_line(gid, spans[gid], ref, ref.tx_gene[tid], keys))
                written_genes.add(gid)
            attrs = p[8]
            if gid != ref.tx_gene[tid]:
                attrs = gid_re.sub(f'gene_id "{gid}"', attrs, count=1)
            oh.write("\t".join([*map(str, p[:8]), attrs]) + "\n")

        # novel transcripts, grouped by gene, in genomic order
        by_gene = defaultdict(list)
        for tid, t in novel_kept.items():
            by_gene[t["gid"]].append(tid)
        order = sorted(by_gene, key=lambda g: (spans[g][0], spans[g][2], g))
        for gid in order:
            tids = sorted(by_gene[gid], key=lambda x: (novel_kept[x]["exons"][0][0], x))
            for tid in tids:
                t = novel_kept[tid]
                attrs = novel_attrs(tid, t, ref, groups_of.get(tid, []), keys)
                if regenerate and gid not in written_genes:
                    gene_attrs = format_attrs(
                        [
                            ("gene_id", gid),
                            ("gene_name", t["gname"]),
                            (keys["gene_type"], t["gtype"]),
                        ]
                    )
                    oh.write(gene_line(gid, spans[gid], ref, gid, keys, gene_attrs))
                    written_genes.add(gid)
                base = [t["chrom"], "GTFforge", None, None, None, ".", t["strand"], "."]
                s, e = t["exons"][0][0], t["exons"][-1][1]
                oh.write("\t".join([*_cols(base, "transcript", s, e), attrs]) + "\n")
                for i, (es, ee) in enumerate(t["exons"], 1):
                    oh.write(
                        "\t".join([*_cols(base, "exon", es, ee), attrs + f' exon_number "{i}";'])
                        + "\n"
                    )


def _cols(base, feature, s, e):
    cols = list(base)
    cols[2], cols[3], cols[4] = feature, str(s), str(e)
    return cols


def novel_attrs(tid, t, ref, groups, keys):
    return format_attrs(
        [
            ("gene_id", t["gid"]),
            ("transcript_id", tid),
            ("gene_name", t["gname"]),
            (keys["gene_type"], t["gtype"]),
            (keys["transcript_type"], f"novel_{t['code']}"),
            ("transcript_name", tid),
            ("class_code", t["code"]),
            ("cmp_ref", t["cref"]),
            ("support_groups", ",".join(groups)),
            ("novel", "1"),
        ]
    )


def nfcore_notes(keys):
    lines = [
        "Using this GTF with nf-core/rnaseq",
        "==================================",
        "--gtf <this file>",
        "Do NOT pass --gencode: novel transcripts do not follow GENCODE's ID scheme.",
        "Do not reuse STAR/Salmon indices built from the original annotation: STAR",
        "projects reads onto the transcripts stored in its index, so novel",
        "transcripts would get no reads. Let nf-core rebuild both (--save_reference).",
    ]
    if keys["gene_type"] != "gene_biotype":
        lines.append(
            f"--featurecounts_group_type {keys['gene_type']}   (nf-core's default, "
            "gene_biotype, is not in this annotation; without it the biotype QC fails)"
        )
    return "\n".join(lines) + "\n"


# ----------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ref", required=True)
    ap.add_argument("--merged", required=True)
    ap.add_argument("--classified", required=True)
    ap.add_argument("--support", required=True)
    ap.add_argument("--groups", nargs="+", required=True)
    ap.add_argument("--keep-classes", required=True)
    ap.add_argument("--keep-contigs-regex", default="")
    ap.add_argument("--ref-stats", default=None, help="prepare_reference.py stats to fold in")
    ap.add_argument("--no-require-support", action="store_true")
    ap.add_argument("--gene-lines", choices=["regenerate", "drop"], default="regenerate")
    ap.add_argument("--out-gtf", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-stats", required=True)
    ap.add_argument("--out-nfcore", required=True)
    a = ap.parse_args(argv)

    contig_re = re.compile(a.keep_contigs_regex) if a.keep_contigs_regex else None

    def keep_contig(c):
        return contig_re is None or contig_re.search(c) is not None

    ref = Reference()
    for p in iter_gtf(a.ref):
        ref.add(p)
    keys = {
        "gene_type": ref.key(GENE_TYPE_KEYS, "gene_type"),
        "transcript_type": ref.key(TRANSCRIPT_TYPE_KEYS, "transcript_type"),
    }
    ref_kept, ref_gene_of, drops, n_ref_split = plan_reference(ref, keep_contig)
    n_ref_in = len(ref.tx_gene)
    if a.ref_stats:
        with open(a.ref_stats) as fh:
            prep = json.load(fh)
        drops.update(prep["dropped"])
        n_ref_in = prep["reference_transcripts_in"]
    ref_tx_locus = {t: next(iter(ref.tx_loci[t])) for t in ref_kept}

    classes = read_classes(a.classified)
    support = read_support(a.support, a.groups)
    novel = read_novel(a.merged, set(ref.tx_gene))
    novel_kept, novel_drops = select_novel(
        novel, classes, support, set(a.keep_classes), keep_contig, not a.no_require_support
    )
    drops.update(novel_drops)
    n_split = assign_genes(novel_kept, ref_gene_of, ref_tx_locus)

    final_gene_of_original = {}
    for tid in ref_kept:
        final_gene_of_original[ref_gene_of[tid]] = ref.tx_gene[tid]
    for t in novel_kept.values():
        if t["novel_gene"]:
            t["gname"], t["gtype"] = t["gid"], "novel"
        else:
            orig = final_gene_of_original[t["gid"]]
            t["gname"] = ref.gene_name.get(orig, orig)
            t["gtype"] = ref.gene_type.get(orig, "")

    groups_of = {tid: [g for g in a.groups if g in support.get(tid, ())] for tid in novel_kept}
    write_gtf(a, ref, ref_kept, ref_gene_of, novel_kept, groups_of, keys, a.out_gtf)

    with open_write(a.out_tsv) as oh:
        oh.write(
            "transcript_id\tgene_id\tgene_name\tgene_type\tclass_code\tcmp_ref\t"
            "support_groups\tchrom\tstart\tend\tstrand\tn_exons\n"
        )
        for tid, t in sorted(
            novel_kept.items(), key=lambda kv: (kv[1]["chrom"], kv[1]["exons"][0][0], kv[0])
        ):
            oh.write(
                "\t".join(
                    map(
                        str,
                        [
                            tid, t["gid"], t["gname"], t["gtype"], t["code"], t["cref"],
                            ",".join(groups_of[tid]), t["chrom"], t["exons"][0][0],
                            t["exons"][-1][1], t["strand"], len(t["exons"]),
                        ],
                    )
                )
                + "\n"
            )

    with open(a.out_nfcore, "w") as fh:
        fh.write(nfcore_notes(keys))

    by_class = Counter(t["code"] for t in novel_kept.values())
    by_group = Counter(g for gs in groups_of.values() for g in gs)
    stats = {
        "attribute_keys": keys,
        "reference_transcripts_in": n_ref_in,
        "reference_transcripts_kept": len(ref_kept),
        "novel_transcripts_in_merge": len(novel),
        "novel_transcripts_kept": len(novel_kept),
        "novel_isoforms_of_reference_genes": sum(not t["novel_gene"] for t in novel_kept.values()),
        "transcripts_in_novel_genes": sum(t["novel_gene"] for t in novel_kept.values()),
        "novel_genes": len({t["gid"] for t in novel_kept.values() if t["novel_gene"]}),
        "novel_genes_split_multi_locus": n_split,
        "reference_genes_split_multi_locus": n_ref_split,
        "kept_by_class": dict(sorted(by_class.items())),
        "kept_by_group": dict(sorted(by_group.items())),
        "dropped": dict(sorted(drops.items())),
    }
    with open(a.out_stats, "w") as fh:
        json.dump(stats, fh, indent=2)
    print(json.dumps(stats, indent=2), file=sys.stderr)


if __name__ == "__main__":
    main()
