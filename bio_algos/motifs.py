from Bio import motifs
from Bio.Seq import Seq

def create_motif(instances):
    return motifs.create(instances)

def print_motif(m):
    return m.alignment.sequences

# How many times each nucleotide shows up in each position
def print_counts(m, letter=None, row=None):
    if (letter and row) is not None:
        return m.counts[letter, row]
    elif (letter is None) and row is not None:
        return m.counts[row]
    elif (row is None) and letter is not None:
        return m.counts[letter]
    else:
        return m.counts

# consensus_type can be 'cons' for regular consensus (default), 'anti' or 'degen'
# identity is the % needed for base to be valid consensus (0-1.0).
# setcase is % needed for uppercase instead of lower (0-1.0)
def print_consensus(m, consensus_type='normal', weighted=False, identity=0, setcase=0):
    if weighted:
        m.counts.calculate_consensus(identity=identity, setcase=setcase)
    # produces seq with most used nucleotide for each position
    if consensus_type is None:
        return m.consensus
    # produces seq with least used nucleotide for each position
    elif consensus_type is 'anti':
        return m.anticonsensus
    # can produce ambiguous nucleotides if tied for min or max count,
    # if so, uses IUPAC nucleotide ambiguity codes: https://droog.gs.washington.edu/parc/images/iupac.html
    elif consensus_type is 'degen':
        return m.degenerate_consensus

# reverses the letters, and interchanges A <-> T and C <-> G. Only works for DNA
def reverse_complement(m):
    r = m.reverse_complement()
    return r
# Shows relative entropy (Measure of how much probability distribution is different from reference distribution)
def relative_entropy(m):
    rel_entropy = m.relative_entropy
    total_rel_entropy = sum(rel_entropy)
    return rel_entropy, total_rel_entropy

def create_weblogo(m, filename):
    m.weblogo(filename)

# set motif ID for jaspar (or other output format)
def set_motif_id(m, identity):
    m.matrix_id = identity

# writes motif in noted formats ('jaspar', 'psm', 'transfac')
def write_motif(m, output_format):
    return format(m, output_format)

def write_multiple_motifs(motifs_list, output_format):
    return motifs.write(motifs_list, output_format)

def normalize_counts(m, pseudocounts):
    return m.counts.normalize(pseudocounts=pseudocounts)

def get_pssm(pwm):
    pssm = pwm.log_odds()
    return pssm

def pssm_stats(pssm):
    return pssm.max, pssm.min, pssm.mean(), pssm.std()

def search_exact_motif(search_seq, m):
    matches = []
    for pos, seq in search_seq.search(m.alignment):
        matches.append((pos, seq))
    return matches

def search_reverse_complement(search_seq, r):
    matches = []
    for pos, seq in search_seq.search(r.alignment):
        matches.append((pos, seq))
    return matches

def search_by_pssm_score(search_seq, pssm):
    matches = []
    for position, score in pssm.search(search_seq, threshold=3.0):
        matches.append((position, score))
    return matches

# set false to 'pos' for false positive, 'neg' for false negative, or 'balance' for balanced threshold
def set_score_thresholds(pssm, background, false=None, threshold=0):
    distribution = pssm.distribution(background=background, precision=10**4)
    if false is 'pos':
        threshold = distribution.threshold_fpr(threshold)
    elif false is 'neg':
        threshold = distribution.threshold_fnr(threshold)
    elif false is 'balance':
        threshold = distribution.threshold_balanced(1000)
