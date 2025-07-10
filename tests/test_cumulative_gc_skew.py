from bio_algos.cumulative_gc_skew import (
    calculate_cumulative_gc_skew,
    cumulative_gc_skew_plot,
)


def test_calculate_cumulative_gc_skew():
    vals = calculate_cumulative_gc_skew("GCGC")
    assert vals == [1, 0, 1, 0]


def test_cumulative_gc_skew_plot(tmp_path):
    out = tmp_path / "cumul.png"
    cumulative_gc_skew_plot("ATGCGC", output=str(out))
    assert out.exists()
