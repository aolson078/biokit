from bio_algos.nucleotide_pie import nucleotide_pie_chart, nucleotide_frequencies


def test_nucleotide_frequencies():
    freqs = nucleotide_frequencies('AATTGGCC')
    assert abs(freqs['A'] - 0.25) < 1e-6
    assert sum(freqs.values()) == 1.0


def test_nucleotide_pie_chart(tmp_path):
    out = tmp_path / 'pie.png'
    nucleotide_pie_chart('ATGCATGC', output=str(out))
    assert out.exists()
