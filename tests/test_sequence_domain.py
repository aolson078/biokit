import time

import pytest

from sequence_workspace.domain import (
    MAX_IMPORT_BYTES,
    SequenceDomainError,
    normalize_sequence,
    parse_import,
    reverse_complement,
    run_analysis,
)


def candidate(sequence, molecule_type="dna"):
    return {
        "client_id": "fixture",
        "source": "manual",
        "sequence": sequence,
        "molecule_type": molecule_type,
    }


def test_normalization_supports_iupac_and_rejects_mixed_t_u():
    normalized = normalize_sequence(" acgt ryswkmbdhvn\n", molecule_hint="dna")
    assert normalized.sequence == "ACGTRYSWKMBDHVN"
    assert reverse_complement(normalized.sequence, normalized.sequence_alphabet) == "NBDHVKMWSRYACGT"
    with pytest.raises(SequenceDomainError, match="both T and U"):
        normalize_sequence("ACTU", molecule_hint="dna")


def test_neutral_manual_sequence_requires_explicit_molecule_type():
    with pytest.raises(SequenceDomainError) as error:
        normalize_sequence("ACGNN")
    assert error.value.code == "molecule_type_required"
    assert normalize_sequence("ACGNN", molecule_hint="rna").molecule_type == "rna"


@pytest.mark.parametrize(
    ("sequence", "molecule_type"),
    [
        ("ACGTRYSWKMBDHVN", "dna"),
        ("ACGURYSWKMBDHVN", "rna"),
    ],
)
def test_full_iupac_reverse_complement_is_an_involution(sequence, molecule_type):
    normalized = normalize_sequence(sequence, molecule_hint=molecule_type)
    reversed_once = reverse_complement(normalized.sequence, normalized.sequence_alphabet)
    assert reverse_complement(reversed_once, normalized.sequence_alphabet) == normalized.sequence


def test_fasta_import_keeps_valid_and_invalid_records_independent():
    records = parse_import(">valid\nATGAAATAA\n>invalid\nATUX", molecule_hint="dna")
    assert records[0]["status"] == "ready"
    assert records[1]["status"] == "error"
    assert records[1]["error"]["code"] == "mixed_alphabet"


def test_import_rejects_oversized_and_non_ascii_whitespace_inputs():
    with pytest.raises(SequenceDomainError) as oversized:
        parse_import("A" * (MAX_IMPORT_BYTES + 1), molecule_hint="dna")
    assert oversized.value.code == "import_too_large"
    result = parse_import("ATG\u00a0TAA", molecule_hint="dna")
    assert result[0]["status"] == "error"
    assert result[0]["error"]["code"] == "invalid_sequence_symbols"


def test_translation_reports_partial_codon_and_orf_uses_first_stop():
    translation = run_analysis("translation", [candidate("ATGAAATAAC")], {"frame": 1, "genetic_code": 1})
    assert translation["protein"] == "MK*"
    assert translation["ignored_terminal_nucleotides"] == 1

    result = run_analysis("orf", [candidate("CCCATGAAATAACCC")], {"min_length": 30})
    assert result["count"] == 0
    result = run_analysis("orf", [candidate("ATG" + "AAA" * 9 + "TAA")], {"min_length": 30})
    assert result["count"] == 1
    assert result["orfs"][0]["length_nt"] == 33


def test_linear_100k_orf_boundary_is_bounded():
    sequence = ("ATG" * 33_332) + "TAA"
    started = time.perf_counter()
    result = run_analysis("orf", [candidate(sequence)], {"min_length": 90})
    elapsed = time.perf_counter() - started
    assert result["count"] == 33_304
    assert len(result["orfs"]) == 100
    assert elapsed < 3.0


def test_pairwise_similarity_counts_alignment_columns():
    result = run_analysis("pairwise_alignment", [candidate("ACGT"), candidate("ACCT")])
    assert result["similarity_percent"] == 75.0


def test_quadratic_comparison_boundary_rejects_oversized_sequence():
    with pytest.raises(SequenceDomainError) as error:
        run_analysis("pairwise_alignment", [candidate("A" * 2001), candidate("A")])
    assert error.value.code == "comparison_too_large"
