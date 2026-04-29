import sys
sys.path.append('.')

from bio_algos.utilities import SequenceTools


def test_find_orfs_detects_simple_orf():
    sequence_tools = SequenceTools()
    sequence = "CCCATGAAATAAGGG"
    orfs = sequence_tools.find_orfs(sequence, min_length=6)
    assert len(orfs) >= 1
    starts = [orf[0] for orf in orfs]
    assert 3 in starts


def test_find_orfs_returns_empty_for_high_threshold():
    sequence_tools = SequenceTools()
    sequence = "ATGAAATAA"
    orfs = sequence_tools.find_orfs(sequence, min_length=300)
    assert orfs == []


def test_find_orfs_requires_terminal_stop_codon():
    sequence_tools = SequenceTools()
    sequence = "ATGAAAAAA"
    orfs = sequence_tools.find_orfs_detailed(sequence, min_length=6)
    assert orfs == []


def test_find_orfs_reports_strand_and_frame():
    sequence_tools = SequenceTools()
    sequence = "CCCATGAAATAAGGG"
    orfs = sequence_tools.find_orfs_detailed(sequence, min_length=6)
    assert len(orfs) >= 1
    first = orfs[0]
    assert first["strand"] in {"+", "-"}
    assert first["frame"] in {1, 2, 3}
