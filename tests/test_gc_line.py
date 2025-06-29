import os
from bio_algos.gc_content_line import calculate_gc_content, gc_line_graph

def test_calculate_gc_content():
    vals = calculate_gc_content("GCGCGC", window=2)
    assert vals[0] == 1.0
    assert len(vals) == 5

def test_gc_line_graph_tmp(tmp_path):
    out = tmp_path / "out.png"
    gc_line_graph("ATGCATGC", output=str(out))
    assert out.exists()
