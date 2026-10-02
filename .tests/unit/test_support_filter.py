import support_filter as sf


def test_min_support_matches_the_zpf68_thresholds():
    # max(2, ceil(n * 0.33)) for the group sizes of the ZPF68 cohort
    assert [sf.min_support(n, 2, 0.33) for n in (6, 8, 10, 12)] == [2, 3, 4, 4]


def test_min_support_floating_point_does_not_round_up():
    assert sf.min_support(12, 1, 0.25) == 3


def test_min_samples_floor_wins_for_small_groups():
    assert sf.min_support(3, 2, 0.33) == 2


def test_select_drops_low_support_and_unstranded_single_exon():
    support = {"multi": 3, "rare": 1, "se_unstranded": 3, "se_stranded": 3}
    structure = {
        "multi": ("+", 2),
        "rare": ("+", 2),
        "se_unstranded": (".", 1),
        "se_stranded": ("-", 1),
    }
    kept, dropped = sf.select(support, structure, needed=2)
    assert kept == {"multi", "se_stranded"}
    assert dropped == {"below_min_support": 1, "single_exon_unstranded": 1}


def test_read_tracking_counts_present_samples(tmp_path):
    p = tmp_path / "cmp.tracking"
    p.write_text(
        "TCONS_1\tXLOC_1\tG|T\t=\tq1:a|b|2\t-\tq3:c|d|2\n"
        "TCONS_2\tXLOC_1\t-\tu\t-\t-\tq3:e|f|2\n"
    )
    support, n = sf.read_tracking(p)
    assert n == 3 and support == {"TCONS_1": 2, "TCONS_2": 1}
