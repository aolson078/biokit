from Bio import AlignIO, SeqIO, Align

import utils

# Calculates GC content, higher GC content implies higher thermal stability
# due to GC pairs having 3 hydrogen bonds instead of AT's 2
from generate_tree import normalize_fasta


def GC_content(nucleotide_string):
    gc_count = 0
    for nuc in nucleotide_string:
        if nuc in ["G", "C"]:
            gc_count += 1
    return gc_count / len(nucleotide_string)


# Takes codon(set of 3 nucleotides), transcribes to Amino Acid. If less than 3 nucleo. left, won't add new Amino Acid
def transcribe_dna(nucleotide_string):
    transcribed = ""
    for i in range(0, len(nucleotide_string) - 1, 3):
        transcribed += utils.DNA_codons[nucleotide_string[i: i + 3]]
    return transcribed


# Takes two or more sequences of nucleotide data, and aligns to identify regions of similarity. Similar regions can
# represent functional, structural, or evolutionary relationships between sequences.
def sequence_alignment(file="./uploads/normal.fasta", file_type="fasta"):
    # Read sequences from file
    sequences = [record.seq for record in SeqIO.parse(file, file_type)]

    # Check that there is multiple sequences
    if len(sequences) < 2:
        raise ValueError("There must be at least two sequences to align")

    # Create pairwise aligner object to score alignment of sequences
    aligner = Align.PairwiseAligner()

    # Align sequences (Only uses first two) !!!!!!!!!!! WILL TAKE FOREVER ON LONG STRINGS.

    alignments = aligner.align(sequences[0], sequences[1])

    for alignment in alignments:
        print(alignment)


# Takes return value from API and returns a tuple with [0] description and [1] nucleotide string
def parse_fasta(fasta):
    output = ([],[])
    for i in range(len(fasta)):
        lines = fasta[i].split('\n')
        description = lines[0]
        output[0].append(description)
        output[1].append("".join(lines[1:]))
    return output

