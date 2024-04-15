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


def hydrophobicity(aminos):
	"""
        Calculates hydrophobicity of amino acid sequence using Eisenberg Hydrophobicity scale
        (Higher value == more hydrophobic)

        :param aminos: Dictionary containing count and relative percentage of amino acid
        :return: Float score representing the hydrophobicity of the amino acid sequence
    """
	Eisenberg_scale = {'A': 0.62, "C": 0.29, "D": -0.9, "E": -0.74, "F": 1.19, "G": 0.48, "H": -0.4,
	                   "I": 1.38, "K": -1.50, "L": 1.06, "M": 0.64, "N": -0.78, "P": 0.12, "Q": -0.85,
	                   "R": -2.53, "S": -.018, "T": -0.05, "V": 1.08, "W": 0.81, "Y": 0.26, "*": 0}

	# Sum up hydrophobicity value for each amino acid
	score = sum(Eisenberg_scale[acid] for acid in aminos)
	return score


def ss_propensity(amino_dicts):
	"""
        Calculates score based on each amino acids propensity to form each secondary structure
        (This function only considers proteins and their α-Helix ands β-Sheet propensities)

        :param amino_dicts: Dictionary containing count and relative percentage of amino acid
        :return: String value of most likely secondary structures formed by the amino acid sequence
    """
	propensities = {'A': ('α-Helix', 5), 'E': ('α-Helix', 5), 'K': ('α-Helix', 5), 'L': ('α-Helix', 5),
	                'M': ('α-Helix', 5),
	                'H': ('α-Helix', 3), 'I': ('α-Helix', 3), 'P': ('α-Helix', 1), 'V': ('β-Sheet', 5),
	                'W': ('β-Sheet', 5),
	                'Y': ('β-Sheet', 5), 'C': ('β-Sheet', 3), 'D': ('β-Sheet', 3), 'Q': ('β-Sheet', 3),
	                'G': ('β-Sheet', 1),
	                'R': ('β-Sheet', 1), 'S': ('β-Sheet', 1), 'T': ('β-Sheet', 1)}
	scores = {'α-Helix': 0, 'β-Sheet': 0}

	for acid, count in amino_dicts.items():
		structure, value = propensities.get(acid, (None, 0))
		if structure:
			scores[structure] += value * count[0]  # Multiply by count of amino acid
	prediction = 'α-Helix' if scores['α-Helix'] > scores['β-Sheet'] else 'β-Sheet'
	return prediction