import gtf_utils


def test_parse_attrs_quoted_unquoted_and_repeated_keys():
    a = gtf_utils.parse_attrs('gene_id "G1"; level 2; tag "basic"; tag "CCDS";')
    assert a == {"gene_id": "G1", "level": "2", "tag": "basic"}


def test_get_attr_does_not_match_key_suffixes():
    text = 'ref_gene_id "R1"; gene_id "G1"; transcript_id "T1";'
    assert gtf_utils.get_attr(text, "gene_id") == "G1"
    assert gtf_utils.get_attr(text, "ref_gene_id") == "R1"
    assert gtf_utils.get_attr(text, "missing") is None


def test_iter_gtf_skips_comments_and_casts_coordinates(tmp_path):
    p = tmp_path / "x.gtf"
    p.write_text('# comment\nchr1\tsrc\texon\t10\t20\t.\t+\t.\tgene_id "G";\n\n')
    rows = list(gtf_utils.iter_gtf(p))
    assert len(rows) == 1 and rows[0][3] == 10 and rows[0][4] == 20


def test_gz_roundtrip(tmp_path):
    p = tmp_path / "x.gtf.gz"
    with gtf_utils.open_write(p) as fh:
        fh.write("a\n")
    assert p.read_bytes()[:2] == b"\x1f\x8b"
    with gtf_utils.open_read(p) as fh:
        assert fh.read() == "a\n"
