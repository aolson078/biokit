from Bio import Phylo, AlignIO, SeqIO, Align
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq

# DNA codon table for nucleotide transcription
DNA_codons = {
	"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
	"TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
	"TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
	"TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
	"CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
	"CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
	"CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
	"CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
	"ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
	"ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
	"AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
	"AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
	"GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
	"GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
	"GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
	"GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}

# RNA codon table for nucleotide transcription
RNA_codons = {
	"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
	"UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
	"UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
	"UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
	"CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
	"CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
	"CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
	"CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
	"AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
	"ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
	"AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
	"AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
	"GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
	"GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
	"GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
	"GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}

# Start codon(s): AUG (can also be GUG (14%) and UUG (3%) in prokaryotes) (AUG can also M if not coding for start)
# Stop codons: UAG UAA UGA (Represented by * in table)
pro_start_codons = ["AUG", "GUG", "UUG"]
eu_start_codons = ["AUG"]
stop_codons = ["UAG", "UAA", "UGA"]


def dna_to_rna(dna_sequence):
	"""
	Converts a DNA sequence to an RNA sequence.
	:param: dna_sequence (str): The input DNA sequence.
	:return: str: The corresponding RNA sequence.
	"""
	# Define the mapping for complementary base pairs
	complement = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}

	# Initialize an empty string to store the RNA sequence
	rna_sequence = ''

	# Iterate over each nucleotide in the DNA sequence
	for nucleotide in dna_sequence:
		# Check if the nucleotide is valid (A, T, C, or G)
		if nucleotide in complement:
			# If valid, replace T with U for RNA
			if nucleotide == 'T':
				rna_sequence += complement['T']
			else:
				rna_sequence += complement[nucleotide]
		else:
			# If invalid, return an error message
			return "Invalid DNA sequence. Only A, T, C, and G are allowed."

	return rna_sequence

def generate_tree(file):
	"""
    Generates and displays phylogenetic tree from aligned sequences in FASTA file
    :param file: Filepath of FASTA file to generate tree from
	"""
	# Normalize sequences from FASTA file
	normalize_fasta(file, "./uploads/normal.fasta")

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

	# Draw phylogenetic tree (Read documentation for customization)
	Phylo.draw(tree)


# Takes FASTA file as input and formats each sequence, so they are the same length, then outputs to a new FASTA
def normalize_fasta(file, output):
	"""
	Formats sequences in FASTA file so they are the same length
	:param file: Filepath to FASTA file
	:param output: Desired filepath to output FASTA file
	"""
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


def GC_content(sequence):
	"""
    Calculates GC content, higher GC content implies higher thermal stability
    due to GC pairs having 3 hydrogen bonds instead of AT's 2
    :param sequence: String sequence of nucleotide bases
    :return: GC percentage
	"""
	gc_count = 0
	for nuc in sequence:
		if nuc in ["G", "C"]:
			gc_count += 1
	return gc_count / len(sequence)


def transcribe_dna(sequence):
	"""
    Takes codon(set of 3 nucleotides), transcribes to Amino Acid.
    If less than 3 nucleo. left, won't add incomplete codon
    :param sequence: String sequence of nucleotide bases
    :return: trasncribed DNA sequence
	"""
	transcribed = ""
	removed = ""
	while len(sequence) % 3 != 0:
		removed += sequence[-1]
		sequence = sequence[:-1]
	for i in range(0, len(sequence) - 1, 3):
		transcribed += DNA_codons[sequence[i: i + 3]]
	return transcribed


def align_sequences(sequences):
	"""
    Aligns two or more sequences using pairwise alignment.
    :param sequences: List of sequences
    :return: Aligned sequences
    """
	if len(sequences) < 2:
		raise ValueError("There must be at least two sequences to align")

	aligner = Align.PairwiseAligner()
	alignments = aligner.align(sequences[0], sequences[1])
	alignment = ""
	for a in alignments:
		alignment = a
	return alignment


def read_sequences_from_file(file="../uploads/normal.fasta", file_type="fasta"):
	"""
    Reads sequences from a file and returns a list of sequence records.
    :param file: Path to the input file
    :param file_type: File format ex: 'fasta', 'genbank'
    :return: List of sequence records
    """
	return [record.seq for record in SeqIO.parse(file, file_type)]



