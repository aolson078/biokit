import random
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt


def align_sequences(sequences, ids):
	"""
        Aligns a list of nucleotide sequences by adding gaps to make them the same length
        :param sequences: List of nucleotide sequences
        :param ids: List of specimen IDs corresponding to the sequences
        :return: Aligned sequences
    """
	# find the maximum length among the sequences
	max_length = max(len(seq) for seq in sequences)

	# create SeqRecord objects for each sequence, adding gaps to make them the same length
	aligned_seqs = []
	for i, (seq, specimen_id) in enumerate(zip(sequences, ids)):
		gaps = '-' * (max_length - len(seq))
		aligned_seq = SeqRecord(Seq(seq + gaps), id=specimen_id)
		aligned_seqs.append(aligned_seq)

	# create a multiple sequence alignment from the aligned sequences
	alignment = MultipleSeqAlignment(aligned_seqs)

	return alignment


def generate_tree(sequences, output_file, ids):
	"""
        Generates and displays phylogenetic tree from a list of nucleotide sequences
        :param sequences: List of nucleotide sequences
        :param output_file: Filepath to save tree
        :param ids: List of specimen IDs corresponding to the sequences
    """
	# align the sequences

	aligned_sequences = align_sequences(sequences, ids)

	# create distance calculator object (Identity model calculates proportion of mismatches in sequence)
	calculator = DistanceCalculator('identity')

	# calculate distance matrix (Represents pairwise distance between sequences)
	distance_matrix = calculator.get_distance(aligned_sequences)

	# construct phylogenetic tree using UPGMA algorithm (unweighted pair group method w/ arithmetic mean)
	# clade is the linear descendants on a phylo tree: https://en.wikipedia.org/wiki/Cladistics
	constructor = DistanceTreeConstructor()
	tree = constructor.upgma(distance_matrix)

	# remove inner node labels
	for clade in tree.find_clades():
		if clade.confidence is not None:
			clade.confidence = None

	colors = ['#' + ''.join([random.choice('0123456789ABCDEF') for _ in range(6)]) for _ in
	          range(len(tree.get_terminals()))]

	# assign colors to each clade
	for clade, color in zip(tree.get_terminals(), colors):
		clade.color = color

	# create a matplotlib figure
	fig = plt.figure(figsize=(10, 8))

	# Uncomment to draw tree
	axes = fig.add_subplot(1, 1, 1)
	Phylo.draw(tree, axes=axes)

	# save the figure as a PNG image file
	fig.savefig(output_file, format='png', dpi=300, bbox_inches='tight')
	plt.close(fig)
