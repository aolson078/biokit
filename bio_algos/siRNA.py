from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.SeqUtils import MeltingTemp, molecular_weight

from bio_algos.utilities import GC_content, dna_to_rna
from sklearn.metrics import jaccard_score


# small interfering RNA design. Silences chosen genes for therapeutic uses and research

# selects sequence within gene for siRNA design
def select_target_sequence(seq, target_length=21):
    # represents gc % range that selected sequence should be in for siRNA design.
    # increase range of bounds if need to get better results (Total range difference shouldn't go above .4)
    lower_bound = .3  # lowering lower bound reduces off-target risk, but reduces stability and specificity
    upper_bound = .5  # raising upper bound increases stability, but reduces stability and increases off-target risk

    gc_bounds = (lower_bound, upper_bound)
    # number of positions to move target sequence so that it falls in the desired gc_content bounds
    offset = 0
    # prevents infinite loop if no suitable target found
    max_offset = len(seq) - target_length

    # if gc_content of target ! in bounds, and there offset < max offset, increment offset and recalculate gc_content
    while (gc_content := GC_content(seq[offset:target_length + offset])) < gc_bounds[0] \
            or gc_content > gc_bounds[1]:
        offset += 1
        if offset > max_offset:
            raise ValueError("No valid target sequence found within the specified GC content range")

    target_sequence = seq[offset:target_length + offset]

    return target_sequence, GC_content(target_sequence)


def design_siRNA(target_sequence):
    # Generate all possible 21-nt siRNA sequences
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


# 1. Design siRNA duplex. Consists of 2 strands
# 1: Sense strand, corresponds to target mRNA, and 2. Antisense, guides RNA-induced silencing complex (RISC)
# Sense strand is same as the target mRNA, Antisense strand is the complement of the sense strand.
# Antisense has 2-nucleotide overhang at 3' end (3rd carbon w/ hydroxyl group) of the antisense (Usually UU or TT)
def create_rna_strands(target_sequence, overhang="UU"):
    complements = {"A": "U", "U": "A", "C": "G", "G": "C"}

    # Create complement RNA (sense strand)
    sense_strand = "".join(complements.get(nuc, nuc) for nuc in target_sequence)

    # Create reverse complement RNA (antisense strand)
    antisense_strand = "".join(complements.get(nuc, nuc) for nuc in reversed(target_sequence))
    antisense_strand += overhang

    return sense_strand, antisense_strand


# 2. Check for off-target effects (occurs when siRNA binds to unintended mRNA sequences
# Consider thermodynamic stability (melting temp) ((Tm)) and seed region complementarity


# 2.1. Sequence Similarity search
def calculate_similarity(sense_strand, antisense_strand):
    sense_list = list(sense_strand)
    antisense_list = list(antisense_strand)
    return jaccard_score(sense_list, antisense_list[:-2], average='micro')


#  2.2 Seed Region Analysis
def seed_region_analysis(sense_strand, utr_sequences, seed_region_lengths=None):
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


# 2.3 Thermodynamic stability
# Na= concentration of Na ions in solution (100 millimolar)  Dnac= dna concentration in nanomolars, higher more stable
def calculate_melting_temp(seq):
    tm = MeltingTemp.Tm_NN(seq, dnac1=100, Na=100, nn_table=MeltingTemp.RNA_NN1, saltcorr=7)
    return tm


# Weight in grams per mole
def calculate_molecular_weight(seq):
    mw = molecular_weight(seq, seq_type="RNA")
    return mw


# Untranslated Regions (region after stop codon)
def get_3utrs(sequence):
    start_codon = 'AUG'
    stop_codons = ['UAA', 'UAG', 'UGA']
    three_utrs = []

    start_index = 0
    while True:
        # Find the start codon index
        start_codon_index = sequence.find(start_codon, start_index)
        if start_codon_index == -1:
            break

        # Find the stop codon index
        stop_codon_index = float('inf')
        for stop_codon in stop_codons:
            temp_index = sequence.find(stop_codon, start_codon_index + 3)
            if temp_index != -1 and temp_index < stop_codon_index:
                stop_codon_index = temp_index

        if stop_codon_index == float('inf'):
            break

        # Extract the 3' UTR
        three_utr = sequence[stop_codon_index + 3:]
        three_utrs.append(three_utr)

        start_index = stop_codon_index + 3

    return three_utrs


# 2.4 Design algos (siDirect, siRNA Wizard, DSIR)  2.5 Experiment Validation (qPCR or microarrays)

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
# Possible scoring methods: MysiRNA Score, Biopredsi, DISR, Thermocomposition21, i-score
# Then filter siRNAs based on thresholds calculated from experimental data

# Perform SNP and Off-Target Filtration
# Check for single nucleotide polymorphisms (SNPs) in target region
# Filter out siRNAs that have significant off-target matches


if __name__ == "__main__":
    rna_sequence = dna_to_rna(
        "GAAGCTGGACAGAGCCGGTTCCTGGAAAGAGCTGGTTCCCTGGCAGGCTGGAGGGCAGGAGCTGGGGCCACGCTGGTCTGGGATAGTTGGGCAGGGAGACGGAGTCTCGATCTGTCACCCAGGCTGGAGTGCAGTGGCACAACCTTGGCTCACTGCAACCTCCGCCTCCCAGGTTCAAGTGATTCTCCTGCCTCAGCCTTCTGAGTAGCTGGAATTACAAGCTGTCTACCTGGTCTCCAGAATGGACGGCCCTGTGGCAGAGCATGCCAAGCAGGAGCCCTTTCACGTGGTCACACCTCTGTTGGAGAGCTGGGCGCTGTCCCAGGTGGCGGGCATGCCTGTCTTCCTCAAGTGTGAGAATGTGCAGCCCAGCGGCTCCTTCAAGATTCGGGGCATTGGGCATTTCTGCCAGGAGATGGCCAAGAAGGGATGCAGACACCTGGTGTGCTCCTCAGGGGGTAATGCGGGCATCGCTGCTGCCTATGCTGCTAGGAAGCTGGGCATTCCTGCCACCATCGTGCTCCCCGAGAGCACCTCCCTGCAGGTGGTGCAGAGGCTGCAGGGGGAGGGGGCCGAGGTTCAGCTGACTGGAAAGGTCTGGGACGAGGCCAATCTGAGGGCGCAAGAGTTGGCCAAGAGGGACGGCTGGGAGAATGTCCCCCCGTTTGACCACCCCCTAATATGGAAAGGCCACGCCAGCCTGGTGCAGGAGCTGAAAGCAGTGCTGAGGACCCCACCAGGTGCCCTGGTGCTGGCAGTTGGGGGTGGGGGTCTCCTGGCCGGGGTGGTGGCTGGCCTGCTGGAGGTGGGCTGGCAGCATGTACCCATCATTGCCATGGAGACCCATGGGGCACACTGCTTCAATGCGGCCATCACAGCCGGCAAGCTGGTCACACTTCCAGACATCACCAGTGTGGCCAAGAGCCTGGGTGCCAAGACGGTGGCCGCTCGGGCCCTGGAGTGCATGCAGGTGTGCAAGATTCACTCTGAAGTGGTGGAGGACACCGAGGCTGTGAGCGCTGTGCAGCAGCTCCTGGATGATGAGCGTATGCTGGTGGAGCCTGCCTGTGGGGCAGCCTTAGCAGCCATCTACTCAGGCCTCCTGCGGAGGCTCCAGGCCGAGGGCTGCCTGCCCCCTTCCCTGACTTCAGTTGTGGTAATCGTGTGTGGAGGCAACAACATCAACAGCCGAGAGCTGCAGGCTTTGAAAACCCACCTGGGCCAGGTCTGAGGGGTCCCATCCTGGCCCCAAAGACCCCTGAGAGGCCCATGGACAGTCCTGTGTCTGGATGAGGAGGACTCAGTGCTGGCAGATGGCAGTGGAAGCTGCCCTGTGCAACTGTGCTGGCTGCCTCCTGAAGGAAGCCCTCCTGGACTGCTTCTTTTGGCTCTCCGACAACTCCGGCCAATAAACACTTTCTGAATTGA")

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

    utrs = get_3utrs(rna_sequence)

    seed_region_analysis(sense, utrs)

    candidates = design_siRNA(rna_sequence)

    print("possible siRNA candidates: ", )

    for candidate in candidates:
        print(search_genome(candidate))

    #get_assemblies(candidate)