import validate_gtf as vg


def _write(tmp_path, rows):
    p = tmp_path / "x.gtf"
    p.write_text("".join("\t".join(map(str, r)) + "\n" for r in rows))
    return p


def _row(chrom, feat, s, e, strand, gid, tid=None):
    attrs = f'gene_id "{gid}";' + (f' transcript_id "{tid}";' if tid else "")
    return [chrom, "src", feat, s, e, ".", strand, ".", attrs]


def test_clean_gtf_passes(tmp_path):
    p = _write(tmp_path, [
        _row("chr1", "gene", 10, 40, "+", "G1"),
        _row("chr1", "transcript", 10, 40, "+", "G1", "T1"),
        _row("chr1", "exon", 10, 20, "+", "G1", "T1"),
        _row("chr1", "exon", 30, 40, "+", "G1", "T1"),
    ])
    stats, bad = vg.validate(p)
    assert not bad.items and stats["transcripts"] == 1


def test_every_problem_class_is_caught(tmp_path):
    p = _write(tmp_path, [
        _row("chr1", "exon", 10, 20, "+", "G1", "T1"),
        _row("chr2", "exon", 10, 20, "+", "G1", "T1"),  # transcript + gene on two loci
        _row("chr1", "exon", 10, 30, "+", "G2", "T2"),
        _row("chr1", "exon", 25, 40, "+", "G2", "T2"),  # overlapping exons
        _row("chr1", "exon", 50, 60, ".", "G3", "T3"),  # no strand
        _row("chr1", "transcript", 70, 90, "+", "G4", "T4"),  # no exons
        _row("chr1", "gene", 1, 5, "+", "G2"),  # gene line does not span T2
    ])
    _, bad = vg.validate(p)
    assert {"transcript_on_two_loci", "gene_on_two_loci", "overlapping_exons", "no_strand",
            "transcript_without_exons", "gene_span_mismatch"} <= set(bad.items)


def test_main_exits_nonzero_on_failure(tmp_path):
    p = _write(tmp_path, [_row("chr1", "exon", 50, 60, ".", "G3", "T3")])
    out = tmp_path / "v.json"
    try:
        vg.main(["--gtf", str(p), "--out", str(out)])
    except SystemExit as e:
        assert e.code == 1
    else:
        raise AssertionError("validation should have failed")
    assert '"ok": false' in out.read_text()
