import matplotlib.pyplot as plt
from typing import List, Optional


def calculate_gc_skew(sequence: str, window: int = 50) -> List[float]:
    """Return GC skew for each sliding window in the sequence."""
    sequence = sequence.upper()
    if window < 1 or window > len(sequence):
        window = len(sequence)
    skew_values: List[float] = []
    for i in range(len(sequence) - window + 1):
        window_seq = sequence[i:i + window]
        g = window_seq.count('G')
        c = window_seq.count('C')
        total = g + c
        skew = (g - c) / total if total else 0.0
        skew_values.append(skew)
    return skew_values


def gc_skew_plot(nuc_string: str, output: Optional[str] = None,
                  window: int = 50, show: bool = False) -> None:
    """Plot GC skew across a nucleotide sequence."""
    skew = calculate_gc_skew(nuc_string, window)
    positions = list(range(1, len(skew) + 1))
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(positions, skew, color='purple', linewidth=1)
    ax.axhline(0, color='black', linestyle='--', linewidth=0.8)
    ax.fill_between(positions, 0, skew, where=[s >= 0 for s in skew],
                    color='violet', alpha=0.3)
    ax.fill_between(positions, 0, skew, where=[s < 0 for s in skew],
                    color='lightblue', alpha=0.3)
    ax.set_xlabel("Position")
    ax.set_ylabel("GC Skew")
    ax.set_title("GC Skew Across Sequence")
    plt.tight_layout()
    if output:
        plt.savefig(output, format='png', dpi=300, bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)
