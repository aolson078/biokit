"""Authoritative sequence-analysis catalog exposed to the client."""

from __future__ import annotations

from typing import Any


def _entry(
    analysis_id: str,
    name: str,
    description: str,
    *,
    status: str = "verified",
    execution: str = "sync",
    min_records: int = 1,
    max_records: int = 1,
    molecule_types: tuple[str, ...] = ("dna", "rna"),
    unavailable_reason: str = "This analysis has not completed scientific verification.",
    parameters: tuple[dict[str, Any], ...] = (),
    units: str = "nucleotides",
    interpretation: str | None = None,
    limitations: str | None = None,
) -> dict[str, Any]:
    return {
        "id": analysis_id,
        "name": name,
        "description": description,
        "status": status,
        "execution": execution,
        "algorithm_version": "1.0.0",
        "min_records": min_records,
        "max_records": max_records,
        "molecule_types": list(molecule_types),
        "max_sequence_length": 100_000,
        "compatible_inputs": {
            "molecule_types": list(molecule_types),
            "min_records": min_records,
            "max_records": max_records,
        },
        "parameters": list(parameters),
        "limits": {
            "max_sequence_length": 100_000,
            "max_records": max_records,
        },
        "units": units,
        "interpretation": interpretation or description,
        "limitations": limitations or (
            "Results must be interpreted in the context of sequence provenance and ambiguity."
            if status == "verified"
            else unavailable_reason
        ),
        "unavailable_reason": unavailable_reason,
    }


ANALYSIS_CATALOG: dict[str, dict[str, Any]] = {
    "summary": _entry("summary", "Sequence summary", "Length, IUPAC composition, ambiguity, and GC content."),
    "reverse_complement": _entry("reverse_complement", "Reverse complement", "Full-IUPAC reverse complement."),
    "transcription": _entry("transcription", "DNA to RNA", "Transcribe a DNA working copy without altering the source.", molecule_types=("dna",)),
    "translation": _entry(
        "translation", "Translation", "Translate a selected frame and genetic code.", execution="adaptive",
        parameters=(
            {"id": "frame", "type": "integer", "minimum": 1, "maximum": 3, "default": 1},
            {"id": "strand", "type": "enum", "values": ["forward", "reverse"], "default": "forward"},
            {"id": "genetic_code", "type": "enum", "values": [1, 2, 3, 11], "default": 1},
        ),
        units="amino_acids",
    ),
    "orf": _entry(
        "orf", "ORF discovery", "Find complete open reading frames across both DNA strands.",
        execution="adaptive", molecule_types=("dna",),
        parameters=(
            {"id": "min_length", "type": "integer", "minimum": 30, "maximum": 5000, "default": 90, "units": "nt"},
            {"id": "genetic_code", "type": "enum", "values": [1, 2, 3, 11], "default": 1},
        ),
    ),
    "nucleotide_composition": _entry("nucleotide_composition", "Nucleotide composition", "Canonical counts plus an ambiguity bucket."),
    "gc_line": _entry("gc_line", "Windowed GC content", "GC percentage across configurable windows."),
    "gc_skew": _entry("gc_skew", "GC skew", "Windowed (G-C)/(G+C)."),
    "cumulative_gc_skew": _entry("cumulative_gc_skew", "Cumulative GC skew", "Running G-C across the sequence."),
    "entropy": _entry("entropy", "Sequence entropy", "Windowed Shannon entropy over observed IUPAC symbols."),
    "kmer_histogram": _entry("kmer_histogram", "3-mer histogram", "Canonical 3-mer counts with skipped ambiguity reporting."),
    "stacked_composition": _entry("stacked_composition", "Comparative composition", "Compare nucleotide composition across saved records.", execution="async", min_records=2, max_records=25),
    "pairwise_alignment": _entry("pairwise_alignment", "Pairwise alignment", "Bounded canonical global alignment and identity.", execution="async", min_records=2, max_records=2),
    "dot_plot": _entry("dot_plot", "Dot plot", "Bounded pairwise dot plot.", status="unavailable", execution="async", min_records=2, max_records=2, unavailable_reason="Private artifact rendering is not yet scientifically verified."),
    "heat_map": _entry("heat_map", "Similarity heat map", "Bounded pairwise similarity heat map.", status="unavailable", execution="async", min_records=2, max_records=2, unavailable_reason="Private artifact rendering is not yet scientifically verified."),
    "melting_temperature": _entry("melting_temperature", "Melting temperature", "Oligonucleotide melting estimate.", status="exploratory"),
    "molecular_weight": _entry("molecular_weight", "Molecular weight", "Estimated molecular weight.", status="exploratory"),
    "sirna": _entry("sirna", "siRNA design", "Candidate design and efficiency scoring.", status="exploratory", molecule_types=("rna",)),
    "sense_similarity": _entry("sense_similarity", "Sense/antisense similarity", "Exploratory strand similarity.", status="exploratory"),
    "phylogenetic_tree": _entry("phylogenetic_tree", "Phylogenetic tree", "Requires a validated multiple-sequence alignment workflow.", status="unavailable", execution="async", min_records=2, max_records=25, unavailable_reason="The legacy padding implementation is not scientifically valid."),
    "hydrophobicity": _entry("hydrophobicity", "Hydrophobicity", "Protein-level hydrophobicity prediction.", status="unavailable", unavailable_reason="The current nucleotide workflow does not provide a validated protein contract."),
    "secondary_structure": _entry("secondary_structure", "Secondary structure", "Protein secondary-structure prediction.", status="unavailable", unavailable_reason="The existing heuristic has not been scientifically validated."),
}


def public_catalog() -> list[dict[str, Any]]:
    return [dict(entry) for entry in ANALYSIS_CATALOG.values()]
