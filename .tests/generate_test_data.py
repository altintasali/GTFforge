#!/usr/bin/env python3
"""Deterministically generate GTFforge's tiny test dataset (config/test.yaml).

Not part of the workflow DAG -- rerun by hand if the design changes:

    python3 .tests/generate_test_data.py

Layout
    .tests/resources/genome.fa        chrT (40 kb) + unplacedT (5 kb)
    .tests/resources/reference.gtf    GENCODE-style: G1..G5 (G5 on unplacedT)
    .tests/resources/te.gtf           TEtranscripts-style rmsk GTF
    .tests/assemblies/{sample}.gtf    StringTie-style assemblies, groups ctrl/trt
    .tests/bams/{sample}.bam          spliced single-end reads, group bam
    .tests/samples.csv

Every novel structure has a known fate (support needed with 3 samples:
max(2, ceil(3 * 0.33)) = 2; with 2 samples: 2):

    N_j     G1 exon-skipping isoform, all groups      kept, class j, gene G1
    N_u     intergenic, trt 2/3                       kept, class u, novel gene
    N_x     antisense to G2, ctrl 2/3                 kept, class x, novel gene
    N_scaf  intergenic on unplacedT, ctrl 3/3         kept (dropped with a contig regex)
    N_rare  intergenic, ctrl 1/3                      dropped: below support
    N_se    single exon, no strand, all samples       dropped: single-exon
    N_ri    G3 with a retained intron, all samples    dropped: class m
    (N_x has its own intron: one shared with G2 on the other strand is class s)
    N_b     intergenic, only in the BAM group         kept, class u (BAM mode works)
"""
import os
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RES, ASM, BAMS = ROOT / "resources", ROOT / "assemblies", ROOT / "bams"

CONTIGS = {"chrT": 40000, "unplacedT": 5000}

# id: (contig, strand, exons, gene_id, gene_name)
REFERENCE = {
    "ENSTT0000001.1": ("chrT", "+", [(1001, 1200), (2001, 2200), (3001, 3200)], "ENSGT0000001.1", "G1"),
    "ENSTT0000002.1": ("chrT", "+", [(6001, 6300), (7001, 7300)], "ENSGT0000002.1", "G2"),
    "ENSTT0000003.1": ("chrT", "-", [(11001, 11200), (12001, 12300), (13001, 13200)], "ENSGT0000003.1", "G3"),
    "ENSTT0000004.1": ("chrT", "+", [(16001, 16500)], "ENSGT0000004.1", "G4"),
    "ENSTT0000005.1": ("unplacedT", "+", [(1001, 1200), (2001, 2200)], "ENSGT0000005.1", "G5"),
}

NOVEL = {
    "N_j": ("chrT", "+", [(1001, 1200), (3001, 3200)]),
    "N_u": ("chrT", "+", [(20001, 20300), (21001, 21300)]),
    "N_x": ("chrT", "-", [(6101, 6250), (6901, 7200)]),
    "N_scaf": ("unplacedT", "+", [(3001, 3200), (3501, 3700)]),
    "N_rare": ("chrT", "+", [(24001, 24200), (25001, 25200)]),
    "N_se": ("chrT", ".", [(27001, 27500)]),
    "N_ri": ("chrT", "-", [(11001, 12300), (13001, 13200)]),
    "N_b": ("chrT", "+", [(30001, 30300), (31001, 31300)]),
}

# which samples assembled each novel structure (GTF groups)
PRESENT = {
    "N_j": {"ctrl_1", "ctrl_2", "ctrl_3", "trt_1", "trt_2", "trt_3"},
    "N_u": {"trt_1", "trt_2"},
    "N_x": {"ctrl_1", "ctrl_3"},
    "N_scaf": {"ctrl_1", "ctrl_2", "ctrl_3"},
    "N_rare": {"ctrl_2"},
    "N_se": {"ctrl_1", "ctrl_2", "ctrl_3", "trt_1", "trt_2", "trt_3"},
    "N_ri": {"ctrl_1", "ctrl_2", "ctrl_3", "trt_1", "trt_2", "trt_3"},
}
GTF_SAMPLES = ["ctrl_1", "ctrl_2", "ctrl_3", "trt_1", "trt_2", "trt_3"]
BAM_SAMPLES = ["bam_1", "bam_2"]
# what the BAM group expresses (reference ids and novel ids)
BAM_EXPRESSED = ["ENSTT0000001.1", "ENSTT0000003.1", "N_j", "N_b"]

TES = [  # contig, start, end, strand, subfamily, family, class
    ("chrT", 19901, 20100, "+", "MT2_Mm", "ERVL", "LTR"),  # at N_u's TSS (sense)
    ("chrT", 7101, 7400, "+", "B2_Mm1a", "B2", "SINE"),  # at N_x's TSS (antisense)
    ("chrT", 35001, 36000, "-", "L1Md_A", "L1", "LINE"),  # nowhere near a transcript
]


def write_genome(rng):
    seqs = {c: "".join(rng.choice("ACGT") for _ in range(n)) for c, n in CONTIGS.items()}
    with open(RES / "genome.fa", "w") as fh:
        for c, s in seqs.items():
            fh.write(f">{c}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i : i + 60] + "\n")
    return seqs


def gtf_line(contig, source, feature, s, e, strand, attrs):
    a = " ".join(f'{k} "{v}";' for k, v in attrs)
    return f"{contig}\t{source}\t{feature}\t{s}\t{e}\t.\t{strand}\t.\t{a}\n"


def write_reference():
    with open(RES / "reference.gtf", "w") as fh:
        fh.write("##description: GTFforge test reference (GENCODE-style attributes)\n")
        for tid, (c, strand, exons, gid, name) in REFERENCE.items():
            g = [("gene_id", gid), ("gene_type", "protein_coding"), ("gene_name", name)]
            fh.write(gtf_line(c, "HAVANA", "gene", exons[0][0], exons[-1][1], strand, g))
            t = [("gene_id", gid), ("transcript_id", tid), ("gene_type", "protein_coding"),
                 ("gene_name", name), ("transcript_type", "protein_coding"),
                 ("transcript_name", f"{name}-201")]
            fh.write(gtf_line(c, "HAVANA", "transcript", exons[0][0], exons[-1][1], strand, t))
            for i, (s, e) in enumerate(exons, 1):
                fh.write(gtf_line(c, "HAVANA", "exon", s, e, strand, [*t, ("exon_number", i)]))


def write_te():
    with open(RES / "te.gtf", "w") as fh:
        for i, (c, s, e, strand, sub, fam, cls) in enumerate(TES, 1):
            fh.write(gtf_line(c, "rmsk", "exon", s, e, strand,
                              [("gene_id", sub), ("transcript_id", f"{sub}_dup{i}"),
                               ("family_id", fam), ("class_id", cls)]))


def write_assembly(sample, rng):
    """A StringTie-like GTF: every reference transcript plus this sample's novel ones."""
    lines = [f"# stringtie test assembly for {sample}\n"]
    items = [(tid, *REFERENCE[tid][:3], tid) for tid in REFERENCE]
    items += [(nid, *NOVEL[nid], None) for nid in NOVEL if sample in PRESENT.get(nid, ())]
    for n, (_key, c, strand, exons, ref_id) in enumerate(items, 1):
        cov = f"{rng.uniform(5, 50):.4f}"
        t = [("gene_id", f"STRG.{n}"), ("transcript_id", f"STRG.{n}.1")]
        if ref_id:
            t += [("reference_id", ref_id)]
        t += [("cov", cov), ("FPKM", cov), ("TPM", cov)]
        lines.append(gtf_line(c, "StringTie", "transcript", exons[0][0], exons[-1][1], strand, t))
        for i, (s, e) in enumerate(exons, 1):
            lines.append(gtf_line(c, "StringTie", "exon", s, e, strand,
                                  [t[0], t[1], ("exon_number", i), ("cov", cov)]))
    with open(ASM / f"{sample}.gtf", "w") as fh:
        fh.writelines(lines)


def structure(key):
    return REFERENCE[key][:3] if key in REFERENCE else NOVEL[key]


def write_bam(sample, rng, read_len=76, reads_per_tx=400):
    import pysam

    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": c, "LN": n} for c, n in CONTIGS.items()]}
    tid_of = {c: i for i, c in enumerate(CONTIGS)}
    reads = []
    for key in BAM_EXPRESSED:
        contig, strand, exons = structure(key)
        # transcript coordinates -> genome blocks
        tx_len = sum(e - s + 1 for s, e in exons)
        for r in range(reads_per_tx):
            off = rng.randrange(0, tx_len - read_len + 1)
            blocks, remaining, pos = [], read_len, off
            for s, e in exons:
                ln = e - s + 1
                if pos >= ln:
                    pos -= ln
                    continue
                take = min(ln - pos, remaining)
                blocks.append((s + pos, take))
                remaining -= take
                pos = 0
                if not remaining:
                    break
            cigar, prev_end = [], None
            for bs, bl in blocks:
                if prev_end is not None:
                    cigar.append((3, bs - prev_end - 1))  # N
                cigar.append((0, bl))  # M
                prev_end = bs + bl - 1
            reads.append((tid_of[contig], blocks[0][0] - 1, cigar, strand, f"{key}_{r}"))
    reads.sort(key=lambda x: (x[0], x[1]))
    out = BAMS / f"{sample}.bam"
    with pysam.AlignmentFile(str(out), "wb", header=header) as bam:
        for ref_id, start, cigar, strand, name in reads:
            a = pysam.AlignedSegment()
            a.query_name = f"{sample}_{name}"
            a.query_sequence = "A" * read_len
            a.flag = 0 if rng.random() < 0.5 else 16  # unstranded library
            a.reference_id = ref_id
            a.reference_start = start
            a.mapping_quality = 60
            a.cigartuples = cigar
            a.query_qualities = pysam.qualitystring_to_array("I" * read_len)
            a.set_tag("NH", 1)
            if any(op == 3 for op, _ in cigar):
                a.set_tag("XS", strand)  # what STAR/HISAT2 add to spliced reads
            bam.write(a)
    pysam.index(str(out))


def main():
    rng = random.Random(42)
    for d in (RES, ASM, BAMS):
        d.mkdir(parents=True, exist_ok=True)
    write_genome(rng)
    write_reference()
    write_te()
    for s in GTF_SAMPLES:
        write_assembly(s, rng)
    for s in BAM_SAMPLES:
        write_bam(s, rng)
    with open(ROOT / "samples.csv", "w") as fh:
        fh.write("sample,group,gtf,bam,strandedness\n")
        for s in GTF_SAMPLES:
            fh.write(f"{s},{s.split('_')[0]},.tests/assemblies/{s}.gtf,,unstranded\n")
        for s in BAM_SAMPLES:
            fh.write(f"{s},bam,,.tests/bams/{s}.bam,unstranded\n")
    print(f"test data written under {os.path.relpath(ROOT)}")


if __name__ == "__main__":
    main()
