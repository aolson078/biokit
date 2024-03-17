from Bio.Blast import NCBIWWW, NCBIXML
from Bio.SeqUtils import MeltingTemp, molecular_weight
from bio_algos.utilities import GC_content, dna_to_rna
from sklearn.metrics import jaccard_score


def select_target_sequence(seq, target_length=21):
	"""
	    Selects a target sequence within the given gene for siRNA design. The selected sequence
	    should have a GC content within a specified range to balance off-target risk, stability,
	    and specificity.

	    :param seq: String representing the gene sequence.
	    :param target_length: Length of the target sequence to be selected. Default is 21.
	    :return: Tuple containing the selected target sequence and its GC content.
	    :raises ValueError: If no valid target sequence is found within the specified GC content range.
	"""
	lower_bound = .1  # lowering lower bound reduces off-target risk, but reduces stability and specificity
	upper_bound = .9  # raising upper bound increases stability, but reduces stability and increases off-target risk

	gc_bounds = (lower_bound, upper_bound)
	# number of positions to move target sequence so that it falls in the desired gc_content bounds
	offset = 0
	# prevents infinite loop if no suitable target found
	max_offset = len(seq) - target_length
	print(seq[offset:target_length + offset])
	print(GC_content(seq[offset:target_length + offset]))
	# if gc_content of target ! in bounds, and there offset < max offset, increment offset and recalculate gc_content
	while (gc_content := GC_content(seq[offset:target_length + offset])) < gc_bounds[0] \
			or gc_content > gc_bounds[1]:
		offset += 1
		if offset > max_offset:
			raise ValueError("No valid target sequence found within the specified GC content range")

	target_sequence = seq[offset:target_length + offset]

	return target_sequence, GC_content(target_sequence)


def design_siRNA(target_sequence):
	"""
	    Generates potential 21-nt siRNA sequences from the given target sequence. The generated sequences
	    follow common design rules for effective siRNA, including starting with 'A' or 'U' and avoiding long GC stretches.

	    :param target_sequence: String representing the target sequence from which to design siRNAs.
	    :return: List of potential siRNA sequences.
	"""
	siRNA_candidates = []
	for i in range(len(target_sequence) - 20):
		candidate = target_sequence[i:i + 21]
		if candidate[0] in 'AU' and candidate[18] in 'GC':
			# Check AU richness in positions 1-7
			if candidate[1:8].count('A') + candidate[1:8].count('U') >= 4:
				# Avoid long GC stretches
				if 'GGG' not in candidate and 'CCC' not in candidate:
					siRNA_candidates.append(candidate)

	return siRNA_candidates


def create_rna_strands(target_sequence, overhang="UU"):
	"""
	    Creates the sense and antisense RNA strands. The sense strand is the same as the target mRNA sequence,
	    the antisense strand is the complement of the
	    sense strand with a 2-nucleotide overhang at the 3' end.

	    :param target_sequence: String representing the target mRNA sequence.
	    :param overhang: String representing the 2-nucleotide overhang to add to the 3' end of the antisense strand. Default is "UU".
	    :return: Tuple containing the sense and antisense RNA strands.
	"""
	complements = {"A": "U", "U": "A", "C": "G", "G": "C"}

	# Create complement RNA (sense strand)
	sense_strand = "".join(complements.get(nuc, nuc) for nuc in target_sequence)

	# Create reverse complement RNA (antisense strand)
	antisense_strand = "".join(complements.get(nuc, nuc) for nuc in reversed(target_sequence))
	antisense_strand += overhang

	return sense_strand, antisense_strand


def calculate_similarity(sense_strand, antisense_strand):
	"""
		Calculates the similarity between the sense and antisense RNA strands using the Jaccard similarity
		coefficient.

		:param sense_strand: String representing the sense RNA strand.
		:param antisense_strand: String representing the antisense RNA strand.
		:return: Float representing the Jaccard similarity coefficient between the sense and antisense strands.
	"""
	sense_list = list(sense_strand)
	antisense_list = list(antisense_strand)
	return jaccard_score(sense_list, antisense_list[:-2], average='micro')


def seed_region_analysis(sense_strand, utr_sequences, seed_region_lengths=None):
	"""
	    Analyzes the seed region of the sense RNA strand for matches in the provided 3' UTR
	    sequences. The seed region is a contiguous subsequence of the sense strand. The function checks for matches of
	    each possible seed region within the 3' UTR sequences.

	    :param sense_strand: String representing the sense RNA strand.
	    :param utr_sequences: List of strings representing the 3' UTR sequences to check for matches.
	    :param seed_region_lengths: List of integers representing the lengths of the seed regions to check.
	                                Default is [6, 7, 8].
	    :return: List of dictionaries with info about the seed regions that have at least one
	             match in the 3' UTR sequences.
	             - 'seed_region': The seed region sequence.
	             - 'seed_length': The length of the seed region.
	             - 'start_pos': The starting position of the seed region within the sense strand.
	             - 'matches': The number of matches found in the 3' UTR sequences.
	"""
	if seed_region_lengths is None:
		seed_region_lengths = [6, 7, 8]
	matches = []

	for seed_length in seed_region_lengths:
		for start_pos in range(len(sense_strand) - seed_length + 1):
			seed_region = sense_strand[start_pos:start_pos + seed_length]

			# Check for matches in the 3' UTR sequences
			seed_matches = 0
			for seq in utr_sequences:
				if seed_region in str(seq):
					seed_matches += 1

			if seed_matches > 0:
				matches.append({
					'seed_region': seed_region,
					'seed_length': seed_length,
					'start_pos': start_pos,
					'matches': seed_matches
				})

	if matches:
		print("Alternative seed regions found:")
		for match in matches:
			print(
				f"Seed region: {match['seed_region']}, Length: {match['seed_length']}, Start position: {match['start_pos']}, Matches: {match['matches']}")
	else:
		print("No alternative seed regions found.")


def calculate_melting_temp(seq):
	"""
		Calculates the melting temp of the RNA sequence using the nearest-neighbor (NN) method using BioPython's
		MeltingTemp class.
		Na= concentration of Na ions in solution (100 millimolar).
		Dnac= dna concentration in nanomolars. Higher values are more stable

		:param seq: String representing the RNA sequence.
		:return: Float representing the melting temperature (Tm) of the given RNA sequence in Celsius.
	"""

	tm = MeltingTemp.Tm_NN(seq, dnac1=100, Na=100, nn_table=MeltingTemp.RNA_NN1, saltcorr=7)
	return tm


def calculate_molecular_weight(seq):
	"""
	    Calculates the molecular weight of the RNA seq with BioPythons molecular_weight function.
	    The molecular weight is calculated based on the average weights of the RNA bases,
	     considering the sequence's length.

	    :param seq: String representing the RNA sequence.
	    :return: Float representing the molecular weight of the given RNA sequence in grams/mole.
	"""
	mw = molecular_weight(seq, seq_type="RNA")
	return mw




def search_genome(sequence):
	"""
	    This function searches the given sequence against the genome and returns a list of similar sequences.

	    :param sequence: The input sequence to search for.
	    :return list: A list of similar sequences found in the genome.
    """
	# Perform a BLAST search against the appropriate genome database
	result_handle = NCBIWWW.qblast("blastn", "nt", sequence)  # Replace "nt" with the desired database

	# Parse the BLAST results
	blast_records = NCBIXML.parse(result_handle)

	# Extract the similar sequences
	similar_sequences = []
	for blast_record in blast_records:
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				# Retrieve the aligned sequence from the genome
				aligned_sequence = hsp.sbjct
				similar_sequences.append(aligned_sequence)

	return similar_sequences


# 3. Evaluate siRNA Efficacy
def predict_efficiency(siRNA_sequence):
	"""
	    Predicts the efficiency of the given siRNA sequence based on GC content, molecular weight,
	    melting temperature, and specific nucleotides at certain positions.

	    :param siRNA_sequence: String representing the siRNA sequence.
	    :return: Integer representing the predicted efficiency score of the given siRNA sequence.
    """
	gc_content = GC_content(siRNA_sequence)
	mole_weight = calculate_molecular_weight(siRNA_sequence)
	tm = calculate_melting_temp(siRNA_sequence)

	# Initialize score
	score = 0

	# GC content
	if 30 <= gc_content < 40:
		score += 1
	elif 40 <= gc_content < 50:
		score += 2
	elif 50 <= gc_content < 60:
		score += 3
	elif gc_content >= 60:
		score += 4

	# molecular weight
	if 20000 <= mole_weight < 26000:
		score += 1
	elif 13000 <= mole_weight < 20000:
		score += 2
	elif mole_weight < 13000:
		score += 3

	# melting temperature
	if 50 <= tm < 60:
		score += 1
	elif 60 <= tm < 120:
		score += 2
	elif tm >= 120:
		score += 3

	# sequence motifs correlated with potency levels (Replace this with motif finding function)
	for i, nucleotide in enumerate(siRNA_sequence):
		if nucleotide == 'A' and i == 5:  # Check for A at position 6 (A6)
			score += 1
		if nucleotide == 'U' and i == 0:  # Check for U at position 1 (U1)
			score -= 1
		if nucleotide == 'G' and i == 18:  # Check for G at position 19 (G19)
			score -= 1
	return score


# -----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
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

	result = select_target_sequence(rna_sequence)
	selected_target = result[0]
	gc = result[1]

	print(f"Selected target sequence: {selected_target}")
	print(f"GC content of target: {gc:.4f}")

	cRNA = create_rna_strands(selected_target)
	sense = cRNA[0]
	antisense = cRNA[1]

	similarity = calculate_similarity(sense, antisense)

	print("Similarity between strands: ", similarity)

	print("Melting temperature: ", calculate_melting_temp(rna_sequence))

	print("Molecular weight: ", calculate_molecular_weight(rna_sequence))

	candidates1 = design_siRNA(rna_sequence)

	print("possible siRNA candidates: ", candidates1)

	efficiency_score = {}

	for candidate1 in candidates1:
		efficiency_score[candidate1] = predict_efficiency(candidate1)

	final_candidates = []

	max_score = max(efficiency_score.values())

	for candidate1, score1 in efficiency_score.items():
		if score1 == max_score:
			final_candidates.append(candidate1)


	print("Most efficient siRNA candidates: ", final_candidates)
	print("These are the optimal targets for siRNA design, they have high target specificity, with the lowest "
	      "off-target effects (unwanted effects on other genes, and longest half-life and stability")
