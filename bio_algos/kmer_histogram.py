import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from typing import Dict, Optional


def kmer_frequencies(seq: str, k: int = 2) -> Dict[str, int]:
    """Count occurrences of each k-mer in a sequence."""
    if k < 1:
        raise ValueError("k must be >= 1")
    seq = seq.upper()
    freqs: Dict[str, int] = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        freqs[kmer] = freqs.get(kmer, 0) + 1
    return freqs


def kmer_histogram(seq: str, k: int = 2, output: Optional[str] = None, show: bool = False) -> None:
    """Plot a histogram of k-mer frequencies for a sequence."""
    freqs = kmer_frequencies(seq, k)
    labels = list(freqs.keys())
    counts = list(freqs.values())
    fig_width = max(6, len(labels) * 0.4)
    fig, ax = plt.subplots(figsize=(fig_width, 4))
    ax.bar(labels, counts, color="teal")
    ax.set_xlabel(f"{k}-mer")
    ax.set_ylabel("Count")
    ax.set_title(f"{k}-mer Frequency")
    plt.xticks(rotation=90)
    plt.tight_layout()
    if output:
        plt.savefig(output, format="png", dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)
