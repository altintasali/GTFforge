import te_tss


def test_interval_index_finds_long_and_short_overlaps():
    idx = te_tss.IntervalIndex()
    idx.add("chr1", 100, 10000, "long")
    idx.add("chr1", 5000, 5010, "short")
    idx.add("chr2", 5000, 5010, "other_chrom")
    idx.build()
    assert sorted(idx.overlapping("chr1", 5005, 5005)) == ["long", "short"]
    assert idx.overlapping("chr1", 10001, 10001) == []
    assert idx.overlapping("chrX", 1, 10) == []


def test_tss_and_first_exon_are_strand_aware():
    idx = te_tss.IntervalIndex()
    idx.add("chr1", 90, 110, ("TE_at_plus_start:F:C", "+"))
    idx.add("chr1", 590, 610, ("TE_at_minus_start:F:C", "+"))
    idx.build()
    fe = {
        "plus": ("chr1", "+", 100, 200),  # TSS 100
        "minus": ("chr1", "-", 500, 600),  # TSS 600
    }
    ann = te_tss.annotate(idx, fe)
    assert ann["plus"][0] == "TE_at_plus_start:F:C(sense)"
    assert ann["minus"][0] == "TE_at_minus_start:F:C(anti)"


def test_first_exon_of_minus_strand_transcript_is_the_last_in_genome_order(tmp_path):
    p = tmp_path / "x.gtf"
    p.write_text(
        'chr1\tGTFforge\texon\t100\t200\t.\t-\t.\ttranscript_id "N1"; novel "1";\n'
        'chr1\tGTFforge\texon\t500\t600\t.\t-\t.\ttranscript_id "N1"; novel "1";\n'
    )
    assert te_tss.first_exons(p, {"N1"}) == {"N1": ("chr1", "-", 500, 600)}
