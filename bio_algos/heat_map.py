import numpy as np
import matplotlib.pyplot as plt
from typing import List, Optional

def lcs_length(seq1: str, seq2: str) -> int:
    """Compute the length of the longest common subsequence (LCS) between two sequences."""
    m, n = len(seq1), len(seq2)
    dp = np.zeros((m + 1, n + 1), dtype=int)
    for i in range(m):
        for j in range(n):
            if seq1[i] == seq2[j]:
                dp[i + 1, j + 1] = dp[i, j] + 1
            else:
                dp[i + 1, j + 1] = max(dp[i, j + 1], dp[i + 1, j])
    return dp[m, n]

def compute_similarity_matrix(nuc_list: List[str]) -> np.ndarray:
    """Compute a normalized pairwise similarity matrix for nucleotide sequences."""
    n = len(nuc_list)
    matrix = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            if i == j:
                matrix[i, j] = 1.0  # Max similarity with self
            else:
                lcs = lcs_length(nuc_list[i], nuc_list[j])
                norm = min(len(nuc_list[i]), len(nuc_list[j]))
                matrix[i, j] = lcs / norm if norm else 0.0
    return matrix

def plot_heatmap(
    matrix: np.ndarray,
    labels: List[str],
    output_file: Optional[str] = None,
    show: bool = False,
    title: str = "Nucleotide Sequence Similarity Heatmap"
):
    """Plot and optionally save/display a similarity heatmap."""
    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(matrix, cmap='viridis')

    # Set ticks and labels
    ax.set_xticks(np.arange(len(labels)))
    ax.set_yticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, fontsize=10, rotation=45, ha='right')
    ax.set_yticklabels(labels, fontsize=10)

    # Annotate cells (optional, can comment if crowded)
    for i in range(len(labels)):
        for j in range(len(labels)):
            ax.text(j, i, f"{matrix[i, j]:.2f}", ha='center', va='center',
                    color='w' if matrix[i, j] < 0.5 else 'black', fontsize=8)

    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("Normalized Similarity", rotation=-90, va="bottom")
    ax.set_title(title)
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, format='png', dpi=500, bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)

def heat_maps(
    nuc_list: List[str],
    organisms: List[str],
    output_file: Optional[str] = None,
    show: bool = False
):
    """
    Generates and saves/displays a heatmap of pairwise similarity between nucleotide sequences.
    """
    if len(nuc_list) != len(organisms):
        raise ValueError("Number of sequences and organism labels must match.")
    if len(nuc_list) < 2:
        raise ValueError("Provide at least two sequences.")

    matrix = compute_similarity_matrix(nuc_list)
    plot_heatmap(matrix, organisms, output_file, show)

# Example usage:
nuc_list = ["ATGCATGC", "ATGCGGGC", "TTTTCATC"]
organisms = ["Organism A", "Organism B", "Organism C"]
heat_maps(nuc_list, organisms, "similarity_heatmap.png", show=True)
