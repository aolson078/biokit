import matplotlib.pyplot as plt
import numpy as np

# https://en.wikipedia.org/wiki/Dot_plot_(bioinformatics)
from bio_algos import utilities

"""
	Creates a graphical representation of two nucleotide sequences and identifies closely related regions.
	After sequence alignment, the sequences are compared base by base and if the bases match, a dot is put on the graph.
	Identical sequences will have a diagonal line running from left to right. 
	
	:param alignment: list of nucleotide sequences to graph
	:param organisms: list of organisms that correspond to the nucleotide sequences
	:param output_file: path of folder to save graph
	:return: dot plot object created with matplotlib and pyplot
"""


def dot_plot(alignment, organisms, output_file):
	# extract the aligned sequences
	seq1 = alignment[0]
	seq2 = alignment[1]

	# create a matrix of zeros len(seq1) x len(seq2)
	matrix = np.zeros((len(seq1), len(seq2)))

	# fill in 1s where bases match
	for i in range(len(seq1)):
		for j in range(len(seq2)):
			if seq1[i] == seq2[j]:
				matrix[i, j] = 1

	# Set figure size in inches
	plt.figure(figsize=(10, 10))
	# create a dot plot from the matrix
	plt.imshow(matrix, cmap='Greys', interpolation='None')

	plt.xticks([0, len(seq2) - 1], ['0', str(len(seq2) - 1)], fontsize=8, rotation=90)  # Display start and end indices
	plt.yticks([0, len(seq1) - 1], ['0', str(len(seq1) - 1)], fontsize=8)  # Display start and end indices

	plt.xlabel(organisms[0])
	plt.ylabel(organisms[1])
	plt.title('Dot Plot Comparison of Sequenced Nucleotides')

	# save plot
	plt.savefig(output_file)

	return output_file

