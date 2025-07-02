import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from typing import Optional, Dict


def nucleotide_frequencies(seq: str) -> Dict[str, float]:
    """Calculate normalized nucleotide frequencies for A, T, G, C."""
    seq = seq.upper()
    counts = {base: seq.count(base) for base in 'ATGC'}
    total = sum(counts.values())
    if total == 0:
        raise ValueError("Sequence contains no nucleotides")
    return {base: counts[base] / total for base in 'ATGC'}


def nucleotide_pie_chart(
    seq: str,
    output: Optional[str] = None,
    show: bool = False
) -> None:
    """Generate a pie chart of nucleotide composition."""
    freqs = nucleotide_frequencies(seq)
    labels = [f"{b} ({freqs[b]*100:.1f}%)" for b in 'ATGC']
    sizes = [freqs[b] for b in 'ATGC']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

    fig, ax = plt.subplots(figsize=(4, 4))
    wedges, texts = ax.pie(sizes, labels=labels, colors=colors, startangle=90)
    ax.axis('equal')
    ax.set_title('Nucleotide Composition')

    if output:
        plt.savefig(output, format='png', dpi=300, bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)
