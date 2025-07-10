import matplotlib.pyplot as plt
from typing import List, Optional


def calculate_cumulative_gc_skew(sequence: str) -> List[int]:
    """Return the cumulative GC skew (G count minus C count) at each position."""
    sequence = sequence.upper()
    g_count = 0
    c_count = 0
    skew: List[int] = []
    for base in sequence:
        if base == 'G':
            g_count += 1
        elif base == 'C':
            c_count += 1
        skew.append(g_count - c_count)
    return skew


def cumulative_gc_skew_plot(
    nuc_string: str,
    output: Optional[str] = None,
    show: bool = False,
) -> None:
    """Plot cumulative GC skew across a nucleotide sequence."""
    skew = calculate_cumulative_gc_skew(nuc_string)
    positions = list(range(1, len(skew) + 1))
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(positions, skew, color='green', linewidth=1)
    ax.set_xlabel("Position")
    ax.set_ylabel("Cumulative GC Skew (G - C)")
    ax.set_title("Cumulative GC Skew")
    plt.tight_layout()
    if output:
        plt.savefig(output, format='png', dpi=300, bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)
