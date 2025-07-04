from bio_algos.kmer_histogram import kmer_frequencies, kmer_histogram


def test_kmer_frequencies_basic():
    freqs = kmer_frequencies("ATAT", k=2)
    assert freqs["AT"] == 2
    assert sum(freqs.values()) == 3


def test_kmer_histogram(tmp_path):
    out = tmp_path / "kmer.png"
    kmer_histogram("ATGCATGC", k=2, output=str(out))
    assert out.exists()
