from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import matplotlib.pyplot as plt

def align_sequences(sequences, ids):
    """
    Aligns a list of nucleotide sequences by adding gaps to make them the same length
    :param sequences: List of nucleotide sequences
    :param ids: List of specimen IDs corresponding to the sequences
    :return: Aligned sequences
    """
    # Find the maximum length among the sequences
    max_length = max(len(seq) for seq in sequences)

    # Create SeqRecord objects for each sequence, adding gaps to make them the same length
    aligned_seqs = []
    for i, (seq, specimen_id) in enumerate(zip(sequences, ids)):
        gaps = '-' * (max_length - len(seq))
        aligned_seq = SeqRecord(Seq(seq + gaps), id=specimen_id)
        aligned_seqs.append(aligned_seq)

    # Create a multiple sequence alignment from the aligned sequences
    alignment = MultipleSeqAlignment(aligned_seqs)

    return alignment


def generate_tree(sequences, output_file, ids):
    """
    Generates and displays phylogenetic tree from a list of nucleotide sequences
    :param sequences: List of nucleotide sequences
    :param output_file: Filepath to save tree
    :param ids: List of specimen IDs corresponding to the sequences
    """
    # Align the sequences
    aligned_sequences = align_sequences(sequences, ids)

    # Create distance calculator object (Identity model calculates proportion of mismatches in sequence)
    calculator = DistanceCalculator('identity')

    # Calculate distance matrix (Represents pairwise distance between sequences)
    distance_matrix = calculator.get_distance(aligned_sequences)
    print("Calculated distance matrix:")
    print(distance_matrix)
    print()
    print("*****************************************")
    print()

    # Construct phylogenetic tree using UPGMA algorithm (unweighted pair group method w/ arithmetic mean)
    # Clade is the linear descendants on a phylo tree: https://en.wikipedia.org/wiki/Cladistics
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(distance_matrix)

    # Remove inner node labels
    for clade in tree.find_clades():
        if clade.confidence is not None:
            clade.confidence = None

    # Draw phylogenetic tree
    fig = plt.figure(figsize=(8, 6))
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=axes, do_show=False)

    # Customize axis labels
    plt.xlabel("Evolutionary Distance")
    plt.ylabel("Sequences")

    print("Calculated tree from distance matrix:")
    print(tree)
    plt.show()
    # Draw phylogenetic tree (Read documentation for customization)
    Phylo.draw(tree)
    Phylo.write(tree, output_file, "newick")




