import finalize as fz


def _ref(lines):
    ref = fz.Reference()
    for chrom, feat, s, e, strand, attrs in lines:
        ref.add([chrom, "src", feat, s, e, ".", strand, ".", attrs])
    return ref


def _exon(tid, gid, n, extra=""):
    return f'gene_id "{gid}"; transcript_id "{tid}"; exon_number "{n}";{extra}'


def test_transcript_reused_on_two_chromosomes_is_dropped():
    ref = _ref([
        ("chr1", "exon", 100, 200, "+", _exon("T1", "G1", 1)),
        ("chr2", "exon", 100, 200, "+", _exon("T1", "G1", 1)),
        ("chr1", "exon", 500, 600, "+", _exon("T2", "G2", 1)),
    ])
    kept, _gene_of, drops, _n_split = fz.plan_reference(ref, lambda c: True)
    assert kept == {"T2"} and drops["ref_transcript_multi_locus"] == 1


def test_transcript_reused_twice_on_one_strand_is_dropped():
    # UCSC refGene style: same id, two non-overlapping copies, exon_number restarts
    ref = _ref([
        ("chr1", "exon", 100, 200, "+", _exon("T1", "G1", 1)),
        ("chr1", "exon", 300, 400, "+", _exon("T1", "G1", 2)),
        ("chr1", "exon", 9100, 9200, "+", _exon("T1", "G1", 1)),
        ("chr1", "exon", 9300, 9400, "+", _exon("T1", "G1", 2)),
    ])
    kept, _, _drops, _ = fz.plan_reference(ref, lambda c: True)
    assert kept == set()


def test_gene_on_two_chromosomes_is_split_per_locus():
    ref = _ref([
        ("chr1", "exon", 100, 200, "+", _exon("T1", "G1", 1)),
        ("chrY", "exon", 100, 200, "+", _exon("T2", "G1", 1)),
    ])
    _kept, gene_of, _, n_split = fz.plan_reference(ref, lambda c: True)
    assert gene_of == {"T1": "G1_chr1+", "T2": "G1_chrY+"} and n_split == 1


def test_contig_filter_drops_reference_transcripts():
    ref = _ref([("GL1.1", "exon", 1, 50, "+", _exon("T1", "G1", 1))])
    kept, _, drops, _ = fz.plan_reference(ref, lambda c: c.startswith("chr"))
    assert not kept and drops["ref_transcript_contig"] == 1


def test_attribute_keys_are_detected():
    gencode = _ref([("chr1", "exon", 1, 50, "+", _exon("T1", "G1", 1, ' gene_type "lncRNA";'))])
    ensembl = _ref([("chr1", "exon", 1, 50, "+", _exon("T1", "G1", 1, ' gene_biotype "lncRNA";'))])
    assert gencode.key(fz.GENE_TYPE_KEYS, "x") == "gene_type"
    assert ensembl.key(fz.GENE_TYPE_KEYS, "x") == "gene_biotype"


def _novel(chrom="chr1", strand="+", locus="MSTRG.1"):
    return {"chrom": chrom, "strand": strand, "locus": locus, "exons": [(10, 20), (30, 40)]}


def test_select_novel_applies_every_filter():
    novel = {"keep": _novel(), "cls": _novel(), "nostrand": _novel(strand="."),
             "contig": _novel(chrom="GL1"), "unsupported": _novel()}
    classes = {"keep": ("u", None), "cls": ("e", None), "nostrand": ("u", None),
               "contig": ("u", None), "unsupported": ("u", None)}
    support = {"keep": {"A"}, "cls": {"A"}, "nostrand": {"A"}, "contig": {"A"}}
    kept, drops = fz.select_novel(novel, classes, support, set("ju"),
                                  lambda c: c.startswith("chr"), True)
    assert set(kept) == {"keep"}
    assert drops == {"novel_class_e": 1, "novel_no_strand": 1, "novel_contig": 1,
                     "novel_no_exact_group_support": 1}


def test_isoforms_join_their_reference_gene_and_others_form_novel_genes():
    kept = {
        "iso": {**_novel(), "code": "j", "cref": "T1"},
        "anti": {**_novel(strand="-"), "code": "x", "cref": "T1"},
        "orphan": {**_novel(), "code": "j", "cref": "T_dropped"},
    }
    fz.assign_genes(kept, {"T1": "G1"}, {"T1": ("chr1", "+")})
    assert kept["iso"]["gid"] == "G1" and not kept["iso"]["novel_gene"]
    assert kept["anti"]["novel_gene"] and kept["orphan"]["novel_gene"]


def test_novel_locus_on_two_strands_is_split():
    kept = {
        "a": {**_novel(strand="+"), "code": "u", "cref": ""},
        "b": {**_novel(strand="-"), "code": "u", "cref": ""},
    }
    assert fz.assign_genes(kept, {}, {}) == 1
    assert {kept["a"]["gid"], kept["b"]["gid"]} == {"MSTRG.1_chr1+", "MSTRG.1_chr1-"}


def test_read_support_uses_group_order_and_exact_matches_only(tmp_path):
    p = tmp_path / "support_map.tracking"
    p.write_text(
        "TCONS_1\tXLOC_1\tMSTRG.1|MSTRG.1.1\t=\tq1:x\t-\tq3:y\n"
        "TCONS_2\tXLOC_1\tMSTRG.1|MSTRG.1.2\tj\tq1:x\t-\t-\n"
    )
    assert fz.read_support(p, ["A", "B", "C"]) == {"MSTRG.1.1": {"A", "C"}}


def test_nfcore_notes_only_ask_for_group_type_when_needed():
    assert "--featurecounts_group_type gene_type" in fz.nfcore_notes(
        {"gene_type": "gene_type", "transcript_type": "transcript_type"})
    assert "featurecounts_group_type" not in fz.nfcore_notes(
        {"gene_type": "gene_biotype", "transcript_type": "transcript_biotype"})
