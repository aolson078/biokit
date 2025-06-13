from Bio import Align, Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from custom_exceptions import InvalidQueryException

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
	complement = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C', '\n': ""}

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


def align_sequences(sequences, ids, seq_type="nucleotide"):
	"""
	    Aligns two sequences using the Needleman-Wunsch algorithm.
	    https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
	    :param sequences: List of sequences
	    :param ids: List of sequence IDs
	    :param seq_type: Type of genetic sequence to run alignment on. "nucleotide", "protein", "genome", or None
	    :return: Multiple sequence alignment
    """
	if len(sequences) < 2:
		raise ValueError("There must be at least two sequences to align")

	# Create a new PairwiseAligner object
	aligner = Align.PairwiseAligner()

	# Adjust the scoring system based on the sequence type (If seq_type==None, the only scoring is +1 for match)
	if seq_type == 'nucleotide':
		aligner.match_score = 3
		aligner.mismatch_score = -3
		aligner.open_gap_score = -7
		aligner.extend_gap_score = -2
	elif seq_type == 'protein':
		aligner.open_gap_score = -10
		aligner.extend_gap_score = -1
	elif seq_type == 'genome':
		aligner.match_score = 1
		aligner.mismatch_score = -1
		aligner.open_gap_score = -5
		aligner.extend_gap_score = -.5

	# Use the aligner to align the sequences
	alignments = aligner.align(sequences[0], sequences[1])
	result = next(alignments)

	# Create SeqRecord objects with appropriate IDs
	seq_records = [
		SeqRecord(Seq(str(result[0])), id=ids[0]),
		SeqRecord(Seq(str(result[1])), id=ids[1])
	]

	# Create a MultipleSeqAlignment object using the SeqRecord objects
	msa = Align.MultipleSeqAlignment(seq_records)

	return msa

# queries the selected database for term and returns the record with the nucleotide string
def fetch_records(query):
	Entrez.email = "aolson078@gmail.com"
	Entrez.apikey = "4a1d5a80996f3691a67335ded0b79299c708"

	# nucleotide, gene, or protein for db
	database = "nucleotide"

	# number of desired responses per query
	return_max = 5

	# fetch nucleotide record related to term
	IDs = Entrez.read(Entrez.esearch(db=database, term=query, field="Organism", retmax=return_max))["IdList"]
	records = []

	for ID in IDs:
		handle = Entrez.efetch(db=database, id=ID, rettype="fasta", retmod="text")
		fasta_record = handle.read()

		# Parse the FASTA record to extract the title and description
		lines = fasta_record.split("\n")
		# Extract the title from the first line
		title = lines[0].strip(">")
		# Join the remaining lines as the description
		description = "\n".join(lines[1:])

		record = {
			"id": ID,
			"title": title,
			"description": description
		}
		records.append(record)

	if len(records) == 0:
		raise InvalidQueryException(query, database)

	if len(records) == -1:
		raise InvalidQueryException(query, database)


	return records
