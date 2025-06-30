import os
from bio_algos.gc_content_heatmap import gc_content_heatmap


def test_gc_content_heatmap_tmp(tmp_path):
    out = tmp_path / "heat.png"
    sequences = ["ATGCATGC", "GGGCCC", "ATATATAT"]
    labels = ["seq1", "seq2", "seq3"]
    gc_content_heatmap(sequences, labels, window=2, output=str(out))
    assert out.exists()
