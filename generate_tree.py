from Bio import Phylo, AlignIO, SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq


def generate_tree(filepath):
	# Normalize sequences from FASTA file
	normalize_fasta(filepath, "./uploads/normal.fasta")

	# Concat sequences from FASTA file
	sequences = AlignIO.read("./uploads/normal.fasta", "fasta")

	# Create distance calculator object (Identity model calculates proportion of mismatches in sequence)
	calculator = DistanceCalculator('identity')

	# Calculate distance matrix (Represents pairwise distance between sequences)
	distance_matrix = calculator.get_distance(sequences)
	print("Calculated distance matrix:")
	print(distance_matrix)

	print()
	print("****************************************************")
	print()

	#  Construct phylogenetic tree using UPGMA algorithm (unweighted pair group method w/ arithmetic mean)
	# Clade is the linear descendants on a phylo tree: https://en.wikipedia.org/wiki/Cladistics
	constructor = DistanceTreeConstructor()
	tree = constructor.upgma(distance_matrix)
	print("Calculated tree from distance matrix:")
	print(tree)

	# Draw phylogenetic tree
	Phylo.draw(tree)


# Takes FASTA file as input and formats each sequence, so they are the same length, then outputs to a new FASTA
def normalize_fasta(input, output):
	# Read sequences from input FASTA file
	sequences = list(SeqIO.parse(input, "fasta"))

	# Find max length of sequences
	max_length = max(len(record.seq) for record in sequences)

	# Normalize the length of each sequence
	for record in sequences:
		if len(record.seq) < max_length:
			# Pad sequence
			record.seq = Seq(str(record.seq).ljust(max_length, '-'))

	# Write padded sequences to new FASTA file
	SeqIO.write(sequences, output, "fasta")
