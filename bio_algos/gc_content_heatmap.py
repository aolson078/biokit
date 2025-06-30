import numpy as np
import matplotlib.pyplot as plt
from typing import List, Optional

from .gc_content_line import calculate_gc_content


def gc_content_heatmap(
    sequences: List[str],
    labels: List[str],
    window: int = 50,
    output: Optional[str] = None,
    show: bool = False,
) -> None:
    """Generate a heatmap of GC content across multiple sequences."""
    if len(sequences) != len(labels):
        raise ValueError("Number of sequences and labels must match.")

    gc_matrix: List[List[float]] = []
    max_len = 0
    for seq in sequences:
        gc_vals = calculate_gc_content(seq, window)
        gc_matrix.append(gc_vals)
        max_len = max(max_len, len(gc_vals))

    # Pad shorter sequences with NaNs for uniform heatmap size
    padded = [vals + [np.nan] * (max_len - len(vals)) for vals in gc_matrix]
    data = np.array(padded)

    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(data, aspect='auto', cmap='viridis', interpolation='none')

    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_xlabel("Window position")
    ax.set_title("GC Content Heatmap")
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("GC fraction", rotation=-90, va="bottom")

    plt.tight_layout()

    if output:
        plt.savefig(output, format='png', dpi=300, bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)
