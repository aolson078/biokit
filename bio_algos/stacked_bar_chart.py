import numpy as np
import matplotlib.pyplot as plt


def count_nucs(nuc_string):
    nuc_counts = {'G': 0, 'A': 0, 'T': 0, 'C': 0}
    for nuc in nuc_string:
        if nuc in nuc_counts:
            nuc_counts[nuc] += 1
    return nuc_counts


def stacked_bar_chart(nuc_list, organisms, output):
    x = ['G', 'A', 'T', 'C']
    y_data = [count_nucs(nuc_string) for nuc_string in nuc_list]
    bottom = np.zeros(4)
    colors = plt.cm.viridis(np.linspace(0, 1, len(nuc_list)))

    fig, ax = plt.subplots(figsize=(10, 8))

    for i, y in enumerate(y_data):
        plt.bar(x, list(y.values()), bottom=bottom, color=colors[i])
        bottom += np.array(list(y.values()))

    ax.set_xlabel('Nucleotides')
    ax.set_ylabel('Frequency')
    ax.set_title('Stacked Bar Graph - Nucleotide Frequencies')
    plt.legend(organisms, loc='upper right')

    plt.tight_layout()

    # save the figure as a PNG image file
    plt.savefig(output, format='png', dpi=300, bbox_inches='tight')
    plt.close(fig)

# stacked_bar_graph(["GATTATACA", "ATTTAGACA", "GGAATTC"])
