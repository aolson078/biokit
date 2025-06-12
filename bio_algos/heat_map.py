import numpy as np
from matplotlib import pyplot as plt


# Heatmap - visualize pairwise similarity between multiple nuc strings using color intensity
def heat_maps(nuc_list, organisms, output_file):
    seq1, seq2 = nuc_list[0], nuc_list[1]
    list_dimensions = ((len(seq1) + 1), (len(seq2) + 1))
    matrix = np.zeros((list_dimensions[0], list_dimensions[1]))

    for i in range(len(seq1)):
        for j in range(len(seq2)):
            if seq1[i] == seq2[j]:
                matrix[i + 1][j + 1] = matrix[i][j] + 1
            else:
                matrix[i + 1][j + 1] = max(matrix[i][j + 1], matrix[i + 1][j]) - 1

    maximum = np.max(matrix)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            matrix[i][j] /= maximum

    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(matrix, cmap='viridis')

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(organisms[1], fontsize=12)
    ax.set_ylabel(organisms[0], fontsize=12)

    colbar = ax.figure.colorbar(im, ax=ax)
    colbar.ax.set_ylabel("Similarity", rotation=-90, va='bottom')

    plt.setp(ax.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
    ax.set_title("Nucleotide Sequence Heat Map")
    plt.tight_layout()

    plt.savefig(output_file, format='png', dpi=500, bbox_inches='tight')
    plt.close(fig)
