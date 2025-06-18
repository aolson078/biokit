import numpy as np
import matplotlib.pyplot as plt

def count_nucs(nuc_string):
    """
    Counts the number of nucleotides in a nucleotide string

    :param nuc_string: string to count nucleotides
    :return: dictionary object with nucleotides as keys and their counts as values.
    """
    nucs = ['G', 'A', 'T', 'C']
    nuc_counts = {nuc: 0 for nuc in nucs}
    for nuc in nuc_string.upper():
        if nuc in nuc_counts:
            nuc_counts[nuc] += 1
    return nuc_counts

def grouped_bar_chart(nuc_list, organisms, output=None, show=False):
    """
    Creates a grouped bar chart of nucleotide counts for different organisms.

    :param nuc_list: List of nucleotide strings for each organism.
    :param organisms: List of organism names, same order as nuc_list.
    :param output: Optional path to save the plot.
    :param show: If True, displays the plot window.
    """
    nucs = ['G', 'A', 'T', 'C']
    data = np.array([[count_nucs(nuc_string)[nuc] for nuc in nucs] for nuc_string in nuc_list])

    x = np.arange(len(nucs))
    width = 0.8 / len(organisms)  # width of each bar

    fig, ax = plt.subplots(figsize=(10, 7))
    colors = plt.cm.viridis(np.linspace(0, 1, len(organisms)))

    for i, (org, counts) in enumerate(zip(organisms, data)):
        ax.bar(x + i*width, counts, width, label=org, color=colors[i])

    ax.set_xticks(x + width * (len(organisms)-1) / 2)
    ax.set_xticklabels(nucs)
    ax.set_xlabel('Nucleotide')
    ax.set_ylabel('Count')
    ax.set_title('Nucleotide Frequency per Organism')
    ax.legend(title='Organism')

    plt.tight_layout()

    if output:
        plt.savefig(output, format='png', dpi=300, bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)

# Example usage:
nuc_list = ['GATTTGCA', 'GATTACA', 'GCGCGCAT']
organisms = ['Human', 'Mouse', 'Yeast']
grouped_bar_chart(nuc_list, organisms, output='nuc_freq.png', show=True)
