from bio_algos.gc_skew import calculate_gc_skew, gc_skew_plot


def test_calculate_gc_skew():
    vals = calculate_gc_skew("GCGC", window=2)
    assert vals[0] == 0.0
    assert len(vals) == 3


def test_gc_skew_plot_tmp(tmp_path):
    out = tmp_path / "skew.png"
    gc_skew_plot("ATGCGCGTAT", output=str(out))
    assert out.exists()
