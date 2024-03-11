from bio_algos.utilities import dna_to_rna
from Bio.Seq import Seq

def amino_acid_composition(seq):
	"""
		Translates RNA to amino acids and calculates the relative amounts of different amino acids in the sequence

		:param seq: String representing the RNA sequence.
		:return: Dictionary containing tuple with count of amino acids and their percentage of the total
	"""

	amino_acid_seq = Seq(seq).translate()

	amino_acids = {acid: (amino_acid_seq.count(acid), amino_acid_seq.count(acid)/len(amino_acid_seq)) for acid in set(amino_acid_seq)}
	return amino_acids

# How water-repellent or attracting an amino acid is using the Eisenberg Hydrophobicity scale
def hydrophobicity(seq):
	pass


# A sequence profile represents the distribution of specific properties (e.g., hydrophobicity, secondary structure propensity) along a protein sequence.
def sequence_profiles(seq):
	pass


rna_sequence = dna_to_rna(
	"GAAGCTGGACAGAGCCGGTTCCTGGAAAGAGCTGGTTCCCTGGCAGGCTGGAGGGCAGGAGCTGGGGCCACGCTGGTCTGGGATAGTTGGGCAGGGAGACGGAGTCTCGAT"
	"CTGTCACCCAGGCTGGAGTGCAGTGGCACAACCTTGGCTCACTGCAACCTCCGCCTCCCAGGTTCAAGTGATTCTCCTGCCTCAGCCTTCTGAGTAGCTGGAATTACAAGC"
	"TGTCTACCTGGTCTCCAGAATGGACGGCCCTGTGGCAGAGCATGCCAAGCAGGAGCCCTTTCACGTGGTCACACCTCTGTTGGAGAGCTGGGCGCTGTCCCAGGTGGCGGG"
	"CATGCCTGTCTTCCTCAAGTGTGAGAATGTGCAGCCCAGCGGCTCCTTCAAGATTCGGGGCATTGGGCATTTCTGCCAGGAGATGGCCAAGAAGGGATGCAGACACCTGGT"
	"GTGCTCCTCAGGGGGTAATGCGGGCATCGCTGCTGCCTATGCTGCTAGGAAGCTGGGCATTCCTGCCACCATCGTGCTCCCCGAGAGCACCTCCCTGCAGGTGGTGCAGAG"
	"GCTGCAGGGGGAGGGGGCCGAGGTTCAGCTGACTGGAAAGGTCTGGGACGAGGCCAATCTGAGGGCGCAAGAGTTGGCCAAGAGGGACGGCTGGGAGAATGTCCCCCCGTT"
	"TGACCACCCCCTAATATGGAAAGGCCACGCCAGCCTGGTGCAGGAGCTGAAAGCAGTGCTGAGGACCCCACCAGGTGCCCTGGTGCTGGCAGTTGGGGGTGGGGGTCTCCT"
	"GGCCGGGGTGGTGGCTGGCCTGCTGGAGGTGGGCTGGCAGCATGTACCCATCATTGCCATGGAGACCCATGGGGCACACTGCTTCAATGCGGCCATCACAGCCGGCAAGCT"
	"GGTCACACTTCCAGACATCACCAGTGTGGCCAAGAGCCTGGGTGCCAAGACGGTGGCCGCTCGGGCCCTGGAGTGCATGCAGGTGTGCAAGATTCACTCTGAAGTGGTGGA"
	"GGACACCGAGGCTGTGAGCGCTGTGCAGCAGCTCCTGGATGATGAGCGTATGCTGGTGGAGCCTGCCTGTGGGGCAGCCTTAGCAGCCATCTACTCAGGCCTCCTGCGGAG"
	"GCTCCAGGCCGAGGGCTGCCTGCCCCCTTCCCTGACTTCAGTTGTGGTAATCGTGTGTGGAGGCAACAACATCAACAGCCGAGAGCTGCAGGCTTTGAAAACCCACCTGGG"
	"CCAGGTCTGAGGGGTCCCATCCTGGCCCCAAAGACCCCTGAGAGGCCCATGGACAGTCCTGTGTCTGGATGAGGAGGACTCAGTGCTGGCAGATGGCAGTGGAAGCTGCCC"
	"TGTGCAACTGTGCTGGCTGCCTCCTGAAGGAAGCCCTCCTGGACTGCTTCTTTTGGCTCTCCGACAACTCCGGCCAATAAACACTTTCTGAATTGA")

print(amino_acid_composition(rna_sequence))