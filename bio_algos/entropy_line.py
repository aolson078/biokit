import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from typing import List, Optional


def shannon_entropy(window_seq: str) -> float:
    """Calculate Shannon entropy of a nucleotide window."""
    window_seq = window_seq.upper()
    length = len(window_seq)
    if length == 0:
        return 0.0
    bases = ['A', 'T', 'G', 'C']
    freqs = [window_seq.count(b) / length for b in bases]
    return -sum(p * np.log2(p) for p in freqs if p > 0)


def calculate_entropy(sequence: str, window: int = 50) -> List[float]:
    """Return entropy values for each sliding window in the sequence."""
    sequence = sequence.upper()
    if window < 1 or window > len(sequence):
        window = len(sequence)
    values = []
    for i in range(len(sequence) - window + 1):
        window_seq = sequence[i:i + window]
        values.append(shannon_entropy(window_seq))
    return values


def entropy_line_graph(
    nuc_string: str,
    output: Optional[str] = None,
    window: int = 50,
    show: bool = False,
) -> None:
    """Plot sequence entropy along a nucleotide sequence."""
    entropies = calculate_entropy(nuc_string, window)
    positions = list(range(1, len(entropies) + 1))
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(positions, entropies, color="orange", linewidth=1)
    ax.set_xlabel("Position")
    ax.set_ylabel("Shannon entropy")
    ax.set_title("Sequence Entropy Across Window")
    plt.tight_layout()
    if output:
        plt.savefig(output, format="png", dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)
