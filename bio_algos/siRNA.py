from Bio.SeqUtils import MeltingTemp, molecular_weight
from bio_algos.utilities import gc_content as calc_gc_content
from sklearn.metrics import jaccard_score
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')


def select_nontarget_sequence(seq, target):
	print("wrong sequence")

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
    lower_bound = .3
    upper_bound = .7
    gc_bounds = (lower_bound, upper_bound)
    
    offset = 0
    max_offset = len(seq) - target_length

    while offset <= max_offset:
        target_candidate = seq[offset:target_length + offset]
        gc_val = calc_gc_content(target_candidate)

        logging.info(f'Checking target sequence at offset {offset}: {target_candidate}, GC_content={gc_val:.2f}')

        if gc_bounds[0] <= gc_val <= gc_bounds[1]:
            logging.info(f'Target sequence selected: {target_candidate}')
            return target_candidate, gc_val

        offset += 1

    logging.error('No valid target sequence found within the specified GC content range')
    raise ValueError("No valid target sequence found within the specified GC content range")


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


def predict_efficiency(siRNA_sequence):
    """Predict the efficiency of the given siRNA sequence."""

    gc_val = calc_gc_content(siRNA_sequence)
    mole_weight = calculate_molecular_weight(siRNA_sequence)
    tm = calculate_melting_temp(siRNA_sequence)

    score = 0

    if 30 <= gc_val < 40:
        score += 1
    elif 40 <= gc_val < 50:
        score += 2
    elif 50 <= gc_val < 60:
        score += 3
    elif gc_val >= 60:
        score += 4

    if 20000 <= mole_weight < 26000:
        score += 1
    elif 13000 <= mole_weight < 20000:
        score += 2
    elif mole_weight < 13000:
        score += 3

    if 50 <= tm < 60:
        score += 1
    elif 60 <= tm < 120:
        score += 2
    elif tm >= 120:
        score += 3

    for i, nucleotide in enumerate(siRNA_sequence):
        if nucleotide == 'A' and i == 5:
            score += 1
        if nucleotide == 'U' and i == 0:
            score -= 1
        if nucleotide == 'G' and i == 18:
            score -= 1

    return score
