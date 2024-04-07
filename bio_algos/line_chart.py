# GC content line graph - displays gc content along len of nuc string. High GC regions correspond to coding regions
from matplotlib import pyplot as plt

from bio_algos.extra.graph import calculate_gc_content


def gc_line_graph(nuc_string, output):
    gc_content = calculate_gc_content(nuc_string, 5)
    # create list for positions and gc content values
    positions = list(range(1, len(gc_content) + 1))

    # create line graph
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(positions, gc_content, linewidth=1)
    ax.set_xlabel("Position")
    ax.set_ylabel("GC content")
    ax.set_title("GC content line")

    plt.tight_layout()

    # save the figure as a PNG image file
    plt.savefig(output, format='png', dpi=300, bbox_inches='tight')
    plt.close(fig)