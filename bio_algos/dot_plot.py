import matplotlib.pyplot as plt
import numpy as np

# https://en.wikipedia.org/wiki/Dot_plot_(bioinformatics)
from bio_algos import utilities

"""
	Creates a graphical representation of two nucleotide sequences and identifies closely related regions.
	After sequence alignment, the sequences are compared base by base and if the bases match, a dot is put on the graph.
	Identical sequences will have a diagonal line running from left to right. 
	
	:param seq1: String representing first nucleotide sequence to compare (Y axis on graph)
	:param seq2: String representing second nucleotide sequence to compare (X axis on graph)
	:return: dot plot object created with matplotlib and pyplot
"""


def dot_plot(alignment):
	# Extract the aligned sequences
	seq1 = alignment[0]
	seq2 = alignment[1]

	# Create a matrix of zeros len(seq1) x len(seq2)
	matrix = np.zeros((len(seq1), len(seq2)))

	# Fill in 1s where bases match
	for i in range(len(seq1)):
		for j in range(len(seq2)):
			if seq1[i] == seq2[j]:
				matrix[i, j] = 1

	# Set figure size in inches
	plt.figure(figsize=(10, 10))

	# Create a dot plot from the matrix
	plt.imshow(matrix, cmap='Greys', interpolation='None')

	plt.xticks([0, len(seq2) - 1], ['0', str(len(seq2) - 1)], fontsize=8, rotation=90)  # Display start and end indices
	plt.yticks([0, len(seq1) - 1], ['0', str(len(seq1) - 1)], fontsize=8)  # Display start and end indices

	plt.xlabel('Sequence 2')
	plt.ylabel('Sequence 1')
	plt.title('Dot Plot Comparison of Sequenced Nucleotides')

	# Show the plot
	plt.show()


# Align the sequences
aligned, score = utilities.align_sequences([
	"GAAGCTGGACAGGTGTCTGGATGAGGAAGCTGCCCTGTGCAACTGTGCTGGCTGCCTCCTAACACTTTCTGAATTGA",
	"GCCATGGGGCACACTGCTTCAATGGCCGCGGAGGACACCGAGGCTGTGAGCGCTGTGCA"])
# Print the alignment and its score
print(aligned)
print("Alignment score: ", score)

# Create a dot plot of the aligned sequences
dot_plot(aligned)
