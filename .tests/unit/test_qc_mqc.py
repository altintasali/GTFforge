import config_used_mqc
import qc_mqc

STATS = """#= Summary for dataset: merged.gtf
#-----------------| Sensitivity | Precision  |
        Base level:   100.0     |    65.7    |
  Transcript level:    80.9     |    28.4    |
"""


def test_gffcompare_stats_are_parsed(tmp_path):
    p = tmp_path / "classify.stats"
    p.write_text(STATS)
    acc = qc_mqc.parse_gffcompare_stats(p)
    assert acc["Transcript level"] == {"Sensitivity": 80.9, "Precision": 28.4}


def test_config_row_labels_have_no_dots_or_underscores():
    assert config_used_mqc.row_label("filter.keep_classes") == "filter > keep classes"


def test_config_flatten_marks_empty_values():
    flat = config_used_mqc.flatten({"ref": {"gtf": "a.gtf", "fasta": ""}})
    assert flat == {"ref.gtf": "a.gtf", "ref.fasta": "(none)"}
