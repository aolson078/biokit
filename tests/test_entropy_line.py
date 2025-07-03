from bio_algos.entropy_line import shannon_entropy, calculate_entropy, entropy_line_graph


def test_shannon_entropy_values():
    assert abs(shannon_entropy('AAAA') - 0.0) < 1e-6
    assert abs(shannon_entropy('ATGC') - 2.0) < 1e-6


def test_calculate_entropy():
    vals = calculate_entropy('ATGCATGC', window=4)
    assert len(vals) == 5
    assert abs(vals[0] - 2.0) < 1e-6


def test_entropy_line_graph(tmp_path):
    out = tmp_path / 'entropy.png'
    entropy_line_graph('ATGCATGC', output=str(out))
    assert out.exists()
