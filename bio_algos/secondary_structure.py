from bio_algos.utilities import dna_to_rna
from Bio.Seq import Seq


def amino_acid_composition(seq):
	"""
		Translates RNA to amino acids and calculates the relative amounts of different amino acids in the sequence

		:param seq: String representing the RNA sequence.
		:return: Dictionary containing tuple with count of amino acids and their percentage of the total
	"""

	amino_acid_seq = Seq(seq).translate()

	amino_acids = {acid: (amino_acid_seq.count(acid), amino_acid_seq.count(acid) / len(amino_acid_seq)) for acid in
	               set(amino_acid_seq)}
	return amino_acids


def hydrophobicity(amino_acid_dict):
	"""
		Calculates hydrophobicity of amino acid sequence using Eisenberg Hydrophobicity scale
		(Higher value == more hydrophobic)

		:param amino_acid_dict: Dictionary containing count of amino acid
		:return: Float score representing the hydrophobicity of the amino acid sequence
	"""
	Eisenberg_scale = {'A': 0.62, "C": 0.29, "D": -0.9, "E": -0.74, "F": 1.19, "G": 0.48, "H": -0.4,
	                   "I": 1.38, "K": -1.50, "L": 1.06, "M": 0.64, "N": -0.78, "P": 0.12, "Q": -0.85,
	                   "R": -2.53, "S": -.018, "T": -0.05, "V": 1.08, "W": 0.81, "Y": 0.26, "*": 0}

	# Sum up hydrophobicity value for each amino acid
	score = sum(Eisenberg_scale[acid] for acid in amino_acid_dict)
	return score


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

amino_acid_dict = amino_acid_composition(rna_sequence)
hydrophobicity_score = hydrophobicity(amino_acid_dict)
