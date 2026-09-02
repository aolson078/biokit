"""Pure sequence normalization and verified analysis contracts."""

from __future__ import annotations

from collections import Counter
from bisect import bisect_right
from dataclasses import dataclass
import heapq
from math import log2
from typing import Any, Iterable
from uuid import uuid4

from Bio import Align
from Bio.Data import CodonTable
from Bio.Seq import Seq

MAX_SEQUENCE_LENGTH = 100_000
MAX_IMPORT_BYTES = 2 * 1024 * 1024
MAX_IMPORT_RECORDS = 25

DNA_IUPAC = frozenset("ACGTRYSWKMBDHVN")
RNA_IUPAC = frozenset("ACGURYSWKMBDHVN")
CANONICAL_DNA = frozenset("ACGT")
CANONICAL_RNA = frozenset("ACGU")
ASCII_WHITESPACE = frozenset(" \t\r\n\f\v")

IUPAC_DNA_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVN",
    "TGCAYRSWMKVHDBN",
)
IUPAC_RNA_COMPLEMENT = str.maketrans(
    "ACGURYSWKMBDHVN",
    "UGCAYRSWMKVHDBN",
)


class SequenceDomainError(ValueError):
    """A typed, user-safe sequence-domain error."""

    def __init__(self, code: str, message: str, *, retryable: bool = False):
        super().__init__(message)
        self.code = code
        self.message = message
        self.retryable = retryable


@dataclass(frozen=True)
class NormalizedSequence:
    sequence: str
    molecule_type: str
    sequence_alphabet: str


def _strip_ascii_whitespace(value: str) -> str:
    return "".join(character for character in value.upper() if character not in ASCII_WHITESPACE)


def normalize_sequence(
    value: str,
    *,
    molecule_hint: str | None = None,
    authoritative_molecule: bool = False,
) -> NormalizedSequence:
    """Normalize one IUPAC nucleotide sequence without changing its source alphabet."""
    if not isinstance(value, str):
        raise SequenceDomainError("invalid_sequence", "Sequence must be text.")

    sequence = _strip_ascii_whitespace(value)
    if not sequence:
        raise SequenceDomainError("empty_sequence", "Sequence cannot be empty.")
    if len(sequence) > MAX_SEQUENCE_LENGTH:
        raise SequenceDomainError(
            "sequence_too_long",
            f"Sequence exceeds the {MAX_SEQUENCE_LENGTH:,} nucleotide limit.",
        )
    if "T" in sequence and "U" in sequence:
        raise SequenceDomainError("mixed_alphabet", "A sequence cannot contain both T and U.")

    if "T" in sequence:
        alphabet = "dna"
        allowed = DNA_IUPAC
    elif "U" in sequence:
        alphabet = "rna"
        allowed = RNA_IUPAC
    else:
        alphabet = "neutral"
        allowed = DNA_IUPAC & RNA_IUPAC

    invalid = sorted(set(sequence) - allowed)
    if invalid:
        raise SequenceDomainError(
            "invalid_sequence_symbols",
            f"Unsupported nucleotide symbols: {', '.join(invalid)}.",
        )

    hint = molecule_hint.lower() if isinstance(molecule_hint, str) else None
    if hint not in {None, "dna", "rna", "unknown"}:
        raise SequenceDomainError("invalid_molecule_type", "Molecule type must be DNA or RNA.")

    if authoritative_molecule and hint in {"dna", "rna", "unknown"}:
        molecule_type = hint
    elif alphabet in {"dna", "rna"}:
        molecule_type = alphabet
    elif hint in {"dna", "rna"}:
        molecule_type = hint
    else:
        raise SequenceDomainError(
            "molecule_type_required",
            "Choose DNA or RNA for a sequence that contains neither T nor U.",
        )

    return NormalizedSequence(sequence, molecule_type, alphabet)


def parse_import(text: str, *, molecule_hint: str | None = None, source_name: str | None = None) -> list[dict[str, Any]]:
    """Parse plain text or FASTA and retain independent record errors."""
    if not isinstance(text, str):
        raise SequenceDomainError("invalid_import", "Import content must be text.")
    if len(text.encode("utf-8")) > MAX_IMPORT_BYTES:
        raise SequenceDomainError("import_too_large", "Import content exceeds 2 MiB.")

    lines = text.splitlines()
    records: list[tuple[str, str]] = []
    if any(line.lstrip().startswith(">") for line in lines):
        header: str | None = None
        sequence_lines: list[str] = []
        for line in lines:
            stripped = line.strip()
            if stripped.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence_lines)))
                header = stripped[1:].strip() or "Untitled FASTA record"
                sequence_lines = []
            elif stripped:
                if header is None:
                    records.append(("Text before first FASTA header", stripped))
                else:
                    sequence_lines.append(stripped)
        if header is not None:
            records.append((header, "".join(sequence_lines)))
    else:
        records.append((source_name or "Pasted sequence", text))

    if len(records) > MAX_IMPORT_RECORDS:
        raise SequenceDomainError(
            "too_many_records",
            f"Import contains more than {MAX_IMPORT_RECORDS} records.",
        )

    results: list[dict[str, Any]] = []
    for index, (title, raw_sequence) in enumerate(records, start=1):
        client_id = f"local:{uuid4()}"
        try:
            normalized = normalize_sequence(raw_sequence, molecule_hint=molecule_hint)
            results.append({
                "status": "ready",
                "client_id": client_id,
                "source": "manual",
                "source_title": title[:500],
                "source_name": source_name,
                "molecule_type": normalized.molecule_type,
                "sequence_alphabet": normalized.sequence_alphabet,
                "sequence": normalized.sequence,
                "length": len(normalized.sequence),
            })
        except SequenceDomainError as error:
            results.append({
                "status": "error",
                "client_id": client_id,
                "record_number": index,
                "source_title": title[:500],
                "error": {"code": error.code, "message": error.message, "retryable": error.retryable},
            })
    return results


def sequence_summary(record: dict[str, Any]) -> dict[str, Any]:
    normalized = normalize_record(record)
    sequence = normalized["sequence"]
    canonical = CANONICAL_DNA if normalized["sequence_alphabet"] != "rna" else CANONICAL_RNA
    counts = Counter(sequence)
    canonical_total = sum(counts[symbol] for symbol in canonical)
    ambiguous_total = len(sequence) - canonical_total
    gc_count = counts["G"] + counts["C"]
    return {
        "length": len(sequence),
        "molecule_type": normalized["molecule_type"],
        "sequence_alphabet": normalized["sequence_alphabet"],
        "symbol_counts": dict(sorted(counts.items())),
        "canonical_count": canonical_total,
        "ambiguous_count": ambiguous_total,
        "gc_percent": round((gc_count / canonical_total) * 100, 4) if canonical_total else None,
    }


def normalize_record(record: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(record, dict):
        raise SequenceDomainError("invalid_record", "Analysis record must be an object.")
    authoritative = record.get("source") == "ncbi"
    normalized = normalize_sequence(
        record.get("sequence", ""),
        molecule_hint=record.get("molecule_type"),
        authoritative_molecule=authoritative,
    )
    result = dict(record)
    result.update({
        "sequence": normalized.sequence,
        "molecule_type": normalized.molecule_type,
        "sequence_alphabet": normalized.sequence_alphabet,
        "length": len(normalized.sequence),
    })
    result.setdefault("client_id", f"local:{uuid4()}")
    return result


def reverse_complement(sequence: str, alphabet: str) -> str:
    table = IUPAC_RNA_COMPLEMENT if alphabet == "rna" else IUPAC_DNA_COMPLEMENT
    return sequence.translate(table)[::-1]


def translate_sequence(record: dict[str, Any], parameters: dict[str, Any]) -> dict[str, Any]:
    normalized = normalize_record(record)
    frame = _bounded_int(parameters.get("frame", 1), "frame", 1, 3)
    genetic_code = _bounded_int(parameters.get("genetic_code", 1), "genetic_code", 1, 11)
    if genetic_code not in {1, 2, 3, 11}:
        raise SequenceDomainError("unsupported_genetic_code", "Genetic code must be 1, 2, 3, or 11.")
    strand = parameters.get("strand", "forward")
    if strand not in {"forward", "reverse"}:
        raise SequenceDomainError("invalid_strand", "Strand must be forward or reverse.")
    if strand == "reverse" and normalized["molecule_type"] != "dna":
        raise SequenceDomainError("unsupported_strand", "Reverse-strand translation is available for DNA only.")

    working = normalized["sequence"]
    if strand == "reverse":
        working = reverse_complement(working, normalized["sequence_alphabet"])
    if "U" in working:
        working = working.replace("U", "T")
    offset = frame - 1
    coding = working[offset:]
    trailing = len(coding) % 3
    complete = coding[:-trailing] if trailing else coding
    protein = str(Seq(complete).translate(table=genetic_code, to_stop=False)) if complete else ""
    return {
        "protein": protein,
        "frame": frame,
        "strand": strand,
        "genetic_code": genetic_code,
        "ignored_terminal_nucleotides": trailing,
    }


def find_orfs(record: dict[str, Any], parameters: dict[str, Any]) -> dict[str, Any]:
    normalized = normalize_record(record)
    if normalized["molecule_type"] != "dna":
        raise SequenceDomainError("unsupported_molecule", "ORF discovery is available for DNA only.")
    genetic_code = _bounded_int(parameters.get("genetic_code", 1), "genetic_code", 1, 11)
    if genetic_code not in {1, 2, 3, 11}:
        raise SequenceDomainError("unsupported_genetic_code", "Genetic code must be 1, 2, 3, or 11.")
    min_length = _bounded_int(parameters.get("min_length", 90), "min_length", 30, 5000)
    table = CodonTable.unambiguous_dna_by_id[genetic_code]
    start_codons = set(table.start_codons)
    stop_codons = set(table.stop_codons)
    original = normalized["sequence"].replace("U", "T")
    candidate_count = 0
    longest: list[tuple[int, int, str, int, int, int]] = []

    for strand, working in (("forward", original), ("reverse", reverse_complement(original, "dna"))):
        for frame in range(3):
            starts: list[int] = []
            for position in range(frame, len(working) - 2, 3):
                codon = working[position:position + 3]
                if codon in start_codons:
                    starts.append(position)
                if codon not in stop_codons:
                    continue
                end = position + 3
                last_eligible = bisect_right(starts, end - min_length)
                candidate_count += last_eligible
                for start in starts[:last_eligible]:
                    length_nt = end - start
                    item = (length_nt, -start, strand, frame + 1, start, end)
                    if len(longest) < 100:
                        heapq.heappush(longest, item)
                    elif item > longest[0]:
                        heapq.heapreplace(longest, item)
                starts = []

    candidates: list[dict[str, Any]] = []
    for length_nt, _, strand, frame, start, end in longest:
        working = original if strand == "forward" else reverse_complement(original, "dna")
        if strand == "forward":
            display_start, display_end = start + 1, end
        else:
            display_start, display_end = len(original) - end + 1, len(original) - start
        protein = str(Seq(working[start:end]).translate(table=genetic_code, to_stop=False))
        candidates.append({
            "start": display_start,
            "end": display_end,
            "strand": strand,
            "frame": frame,
            "length_nt": length_nt,
            "length_aa": max(0, len(protein) - 1),
            "protein": protein[:-1] if protein.endswith("*") else protein,
        })
    candidates.sort(key=lambda item: (-item["length_nt"], item["start"], item["strand"], item["frame"]))
    return {"count": candidate_count, "orfs": candidates, "truncated": candidate_count > 100}


def windowed_metrics(record: dict[str, Any], parameters: dict[str, Any], kind: str) -> dict[str, Any]:
    normalized = normalize_record(record)
    sequence = normalized["sequence"]
    if len(sequence) < 10:
        raise SequenceDomainError("sequence_too_short", "Windowed analysis requires at least 10 nucleotides.")
    default_window = min(100, len(sequence))
    window = _bounded_int(parameters.get("window", default_window), "window", 10, min(10_000, len(sequence)))
    step = _bounded_int(parameters.get("step", max(1, window // 2)), "step", 1, window)
    points: list[dict[str, Any]] = []
    for start in range(0, len(sequence) - window + 1, step):
        segment = sequence[start:start + window]
        counts = Counter(segment)
        if kind == "gc_line":
            denominator = sum(counts[base] for base in ("A", "C", "G", "T", "U"))
            value = ((counts["G"] + counts["C"]) / denominator) * 100 if denominator else None
        elif kind == "gc_skew":
            denominator = counts["G"] + counts["C"]
            value = (counts["G"] - counts["C"]) / denominator if denominator else None
        else:
            probabilities = [count / len(segment) for count in counts.values()]
            value = -sum(probability * log2(probability) for probability in probabilities)
        points.append({"start": start + 1, "end": start + window, "value": value})
    return {"window": window, "step": step, "points": points}


def run_analysis(analysis_id: str, records: Iterable[dict[str, Any]], parameters: dict[str, Any] | None = None) -> dict[str, Any]:
    from .catalog import ANALYSIS_CATALOG

    parameters = parameters or {}
    normalized_records = [normalize_record(record) for record in records]
    entry = ANALYSIS_CATALOG.get(analysis_id)
    if entry is None:
        raise SequenceDomainError("unknown_analysis", "Unknown analysis.")
    if entry["status"] != "verified":
        raise SequenceDomainError("analysis_unavailable", entry["unavailable_reason"])
    if len(normalized_records) < entry["min_records"] or len(normalized_records) > entry["max_records"]:
        raise SequenceDomainError("invalid_record_count", "This analysis does not support the selected record count.")
    if any(record["molecule_type"] not in entry["molecule_types"] for record in normalized_records):
        raise SequenceDomainError("unsupported_molecule", "This analysis is not compatible with the selected molecule type.")

    record = normalized_records[0]
    if analysis_id == "summary":
        return sequence_summary(record)
    if analysis_id == "reverse_complement":
        return {"sequence": reverse_complement(record["sequence"], record["sequence_alphabet"])}
    if analysis_id == "transcription":
        if record["molecule_type"] != "dna":
            raise SequenceDomainError("unsupported_molecule", "Transcription requires DNA.")
        return {"sequence": record["sequence"].replace("T", "U"), "source_unchanged": True}
    if analysis_id == "translation":
        return translate_sequence(record, parameters)
    if analysis_id == "orf":
        return find_orfs(record, parameters)
    if analysis_id in {"gc_line", "gc_skew", "entropy"}:
        return windowed_metrics(record, parameters, analysis_id)
    if analysis_id == "cumulative_gc_skew":
        running = 0
        points = []
        for index, symbol in enumerate(record["sequence"], start=1):
            running += 1 if symbol == "G" else -1 if symbol == "C" else 0
            points.append({"position": index, "value": running})
        return {"points": points}
    if analysis_id == "nucleotide_composition":
        counts = Counter(record["sequence"])
        canonical = {base: counts[base] for base in ("A", "C", "G", "T", "U") if counts[base]}
        canonical["Ambiguous"] = len(record["sequence"]) - sum(canonical.values())
        return {"counts": canonical}
    if analysis_id == "kmer_histogram":
        canonical = CANONICAL_RNA if record["sequence_alphabet"] == "rna" else CANONICAL_DNA
        counts: Counter[str] = Counter()
        skipped = 0
        for index in range(len(record["sequence"]) - 2):
            kmer = record["sequence"][index:index + 3]
            if set(kmer) <= canonical:
                counts[kmer] += 1
            else:
                skipped += 1
        return {"k": 3, "counts": dict(counts.most_common()), "skipped": skipped}
    if analysis_id == "stacked_composition":
        return {"records": [
            {"client_id": item["client_id"], **run_analysis("nucleotide_composition", [item])}
            for item in normalized_records
        ]}
    if analysis_id == "pairwise_alignment":
        first, second = normalized_records
        _require_canonical_pair(first, second)
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 1
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.5
        alignment = aligner.align(first["sequence"], second["sequence"])[0]
        coordinates = alignment.coordinates
        matches = 0
        aligned_length = 0
        for index in range(coordinates.shape[1] - 1):
            first_start, first_end = int(coordinates[0, index]), int(coordinates[0, index + 1])
            second_start, second_end = int(coordinates[1, index]), int(coordinates[1, index + 1])
            first_length = first_end - first_start
            second_length = second_end - second_start
            aligned_length += max(first_length, second_length)
            if first_length and second_length:
                matches += sum(
                    left == right
                    for left, right in zip(
                        first["sequence"][first_start:first_end],
                        second["sequence"][second_start:second_end],
                    )
                )
        return {
            "score": alignment.score,
            "similarity_percent": (matches / max(aligned_length, 1)) * 100,
            "alignment": str(alignment),
        }
    raise SequenceDomainError("analysis_not_implemented", "This verified analysis is not implemented.")


def _require_canonical_pair(first: dict[str, Any], second: dict[str, Any]) -> None:
    if len(first["sequence"]) > 2000 or len(second["sequence"]) > 2000:
        raise SequenceDomainError("comparison_too_large", "Pairwise analysis is limited to 2,000 nt per sequence.")
    if len(first["sequence"]) * len(second["sequence"]) > 4_000_000:
        raise SequenceDomainError("comparison_too_large", "Pairwise analysis exceeds 4,000,000 comparison cells.")
    for record in (first, second):
        canonical = CANONICAL_RNA if record["sequence_alphabet"] == "rna" else CANONICAL_DNA
        if not set(record["sequence"]) <= canonical:
            raise SequenceDomainError("unsupported_ambiguity", "Pairwise analysis requires canonical bases only.")


def _bounded_int(value: Any, name: str, minimum: int, maximum: int) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as error:
        raise SequenceDomainError("invalid_parameter", f"{name} must be an integer.") from error
    if parsed < minimum or parsed > maximum:
        raise SequenceDomainError("invalid_parameter", f"{name} must be between {minimum} and {maximum}.")
    return parsed
