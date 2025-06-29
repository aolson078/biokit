import matplotlib.pyplot as plt
from typing import List, Optional


def calculate_gc_content(sequence: str, window: int = 50) -> List[float]:
    """Return GC fraction for each sliding window in the sequence."""
    sequence = sequence.upper()
    if window < 1 or window > len(sequence):
        window = len(sequence)
    values = []
    for i in range(len(sequence) - window + 1):
        window_seq = sequence[i:i + window]
        gc = (window_seq.count('G') + window_seq.count('C')) / window
        values.append(gc)
    return values


def gc_line_graph(nuc_string: str, output: Optional[str] = None, window: int = 50, show: bool = False) -> None:
    """Plot GC content along a sequence."""
    gc_content = calculate_gc_content(nuc_string, window)
    positions = list(range(1, len(gc_content) + 1))
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(positions, gc_content, linewidth=1)
    ax.set_xlabel("Position")
    ax.set_ylabel("GC content")
    ax.set_title("GC Content Across Sequence")
    plt.tight_layout()
    if output:
        plt.savefig(output, format='png', dpi=300, bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)
