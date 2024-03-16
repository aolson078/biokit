from Bio import motifs
from Bio.Seq import Seq

# Motif with example instances (set create param alphabet to use for RNA)
instances = [Seq("AACGT"), Seq("AGGCT"), Seq("ATTGC"), Seq("GAGTC"), Seq("TTGCA")]
m = motifs.create(instances)
instances = [Seq("AAAGT"), Seq("AAGCT"), Seq("AGGGC"), Seq("GAGAC"), Seq("TTGAA")]
n = motifs.create(instances)

print(m)
print(m.alignment.sequences)

print("--------------------------------------------------")
# How many times each nucleotide shows up in each position
print(m.counts)

# How many times 'A' shows up in each pos
print(m.counts['A'])

# How many times each nucleotide shows up in pos 3
print(m.counts[:, 3])

# How many times "T" shows up in pos 1
print(m.counts["T", 1])

print("--------------------------------------------------")
# Produces seq with most used nucleotide for each position
print("Consensus:", m.consensus)

# Produces seq with least used nucleotide for each position
print("Anti-consensus:", m.anticonsensus)

# Can produce ambiguous nucleotides if tied for min or max count,
# if so, uses IUPAC nucleotide ambiguity codes: https://droog.gs.washington.edu/parc/images/iupac.html
print("Degen-consensus:", m.degenerate_consensus)

print("--------------------------------------------------")
print(m.consensus)
# Identity % required to define consensus. .5 = 50% values must be nuc or value replaced with N or X (nuc or amino)
# Setcase  % threshold for letter case.  .7 = 70% of values must be nuc for letter to be uppercase
print(m.counts.calculate_consensus(identity=0.5, setcase=0.7))

print("--------------------------------------------------")
# Reverses the letters, and interchanges A <-> T and C <-> G
# Only works for DNA
r = m.reverse_complement()
print("Reverse complement: ", r.consensus)
print("--------------------------------------------------")
# Shows relative entropy (Measure of how much probability distribution is different from reference distribution)
print("Relative entropy matrix: ", m.relative_entropy)

# Relative entropy of entire motif
print("Total relative entropy: ", sum(m.relative_entropy))

# Graphic representation of sequence alignment. Height of stack indicates sequence conservation at that pos.
# Height of symbols in each stack indicate relative frequency of each amino or nucleic acid at that pos.
# Provides description of binding site
m.weblogo("./motifs/motif.png")

print("--------------------------------------------------")
# Set motif ID for jaspar
m.matrix_id = "My first motif :)"

# Writes motif in noted formats
print(format(m, "pfm"))
print(format(m, "jaspar"))
print(format(m, "transfac"))
print("--------------------------------------------------")
# Write multiple motifs
two_motifs = [m, n]
print(motifs.write(two_motifs, "jaspar"))

print("--------------------------------------------------")
# Normalize position weight probabilities (psuedocounts param is partly so we don't divide by 0)
# Can set background weight if not standard, ex: human genome gc content is about 40%
# pwm = m.counts.normalize(pseudocounts={"A": 0.6, "C": 0.4, "G": 0.4, "T": 0.6})

pwm = m.counts.normalize(pseudocounts=0.1)
print(pwm)
# Can derive consensus from pos weight matrix, degenerate might be different from the original due to psuedocounts
print(pwm.consensus)
print(pwm.anticonsensus)
print(pwm.degenerate_consensus)

print("--------------------------------------------------")
# Position specific scoring matrices
# positive values = more frequent in sample distribution, negative values = more frequent in standard distribution
pssm = pwm.log_odds()
print("Position-Specific Scoring Matrix:")
print(pssm)

# Can set specific background as before (40% gc content for human genome:
print("Position-Specific Scoring Matrix (Weighted for human genome):")
background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
pssm1 = pwm.log_odds(background)
print(pssm1)

print("Maximum obtainable score: %4.2f" % pssm.max)
print("Minimum obtainable score: %4.2f" % pssm.min)

print("Mean: %0.2f" % pssm.mean())
print("Standard deviation: %0.2f" % pssm.std())

print("--------------------------------------------------")

# 17.6
# Searches search_seq for substring and returns tuple with position(s) and sequence(s) if found
def search_exact_motif(search_seq=Seq("TACACTGCATTACAACCCAAGCATTA")):
    matches = []
    for pos, seq in search_seq.search(m.alignment):
        matches.append((pos, seq))
    return matches


def search_reverse_complement(search_seq=Seq("TACACTGCATTACAACCCAAGCATTA")):
    matches = []
    for pos, seq in search_seq.search(r.alignment):
        matches.append((pos, seq))
    return matches


# Searches search_seq for substring with position specific score matrix threshold (in log2, so 3 == x8)
def search_by_pssm_score(search_seq=Seq("TACACTGCATTACAACCCAAGCATTA")):
    for position, score in pssm.search(search_seq, threshold=3.0):
        print("Position %d: score = %5.3f" % (position, score))


print(search_exact_motif())
print("--------------------------")
print(search_reverse_complement())
print("--------------------------")
print(search_by_pssm_score())
print("--------------------------")


print("--------------------------")
# Setting score thresholds (score space distribution grows expo. w/ seq length, so we need to limit it)

distribution = pssm.distribution(background=background, precision=10**4)

# Set false_pos rate (Prob of finding motif in background generated sequence)
threshold = distribution.threshold_fpr(0.01)
print("%5.3f" % threshold)

# Set false_neg (Prob of not finding instance generated from motif)
threshold = distribution.threshold_fnr(0.1)
print("%5.3f" % threshold)

# Set treshold to relation between neg and pos
threshold = distribution.threshold_balanced(1000)
print("%5.3f" % threshold)