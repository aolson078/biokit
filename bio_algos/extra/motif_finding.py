"""motif_analysis.py
====================
A concise, modular demonstration of common Biopython *motifs* operations.

This script rewrites the original single-block example into self-contained, re-usable
functions and a simple CLI so you can:

* build one or more motifs from example instances;
* inspect counts, consensus, reverse complements and information content;
* derive PWM/PSSM objects with custom pseudocounts or background models;
* search a query sequence on either strand using an exact match **or** a score
  threshold; and
* dump the motifs to popular formats (PFM, JASPAR, TRANSFAC) in a single call.

Run `python motif_analysis.py --help` to see all options.
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Iterable, List, Tuple

from Bio import motifs
from Bio.Seq import Seq

# -----------------------------------------------------------------------------
# Configuration & logging
# -----------------------------------------------------------------------------

LOG_FMT = "%(asctime)s | %(levelname)-8s | %(message)s"
logging.basicConfig(level=logging.INFO, format=LOG_FMT)
logger = logging.getLogger(__name__)

# -----------------------------------------------------------------------------
# Example data (replace with your own as required)
# -----------------------------------------------------------------------------

EXAMPLE_INSTANCES_A: List[Seq] = [
    Seq("AACGT"),
    Seq("AGGCT"),
    Seq("ATTGC"),
    Seq("GAGTC"),
    Seq("TTGCA"),
]

EXAMPLE_INSTANCES_B: List[Seq] = [
    Seq("AAAGT"),
    Seq("AAGCT"),
    Seq("AGGGC"),
    Seq("GAGAC"),
    Seq("TTGAA"),
]

# -----------------------------------------------------------------------------
# Motif utilities
# -----------------------------------------------------------------------------

def build_motif(instances: Iterable[Seq], motif_id: str | None = None) -> motifs.Motif:
    """Create a *motifs.Motif* object from an iterable of sequences."""
    m = motifs.create(list(instances))
    if motif_id:
        m.matrix_id = motif_id
    return m


def motif_summary(m: motifs.Motif) -> str:
    """Return a human-readable multi-line summary of *m*."""
    lines = [
        f"Motif length      : {m.length}",
        f"Alphabet          : {m.alphabet}",
        f"Instances         : {len(m.instances)}",
        "--- consensus ---",
        f"  consensus       : {m.consensus}",
        f"  anti-consensus  : {m.anticonsensus}",
        f"  degenerate      : {m.degenerate_consensus}",
        "--- counts matrix ---",
        str(m.counts),
        "--- relative entropy / bit score ---",
        f"  per-pos bits    : {[round(x, 3) for x in m.relative_entropy]}",
        f"  total bits      : {round(sum(m.relative_entropy), 3)}",
    ]
    return "\n".join(lines)


def build_pwm_pssm(
    m: motifs.Motif,
    *,
    pseudocount: float | dict[str, float] = 0.1,
    background: dict[str, float] | None = None,
):
    """Return (PWM, PSSM) pair with given pseudocounts and background model."""
    pwm = m.counts.normalize(pseudocounts=pseudocount)
    pssm = pwm.log_odds(background)
    return pwm, pssm


# -----------------------------------------------------------------------------
# Search helpers
# -----------------------------------------------------------------------------

def search_exact(m: motifs.Motif, query: Seq) -> List[Tuple[int, Seq]]:
    """Find exact occurrences of *m* in *query* (forward strand only)."""
    return [(pos, seq) for pos, seq in query.search(m.alignment)]


def search_pssm(
    pssm: motifs.matrix.PositionSpecificScoringMatrix,
    query: Seq,
    threshold: float = 3.0,
    both_strands: bool = True,
) -> List[Tuple[int, float]]:
    """Return [(position, score)] hits with *score >= threshold*.

    The *threshold* is in log₂ space; 3.0 means 2³ = 8-fold enrichment over
    background.
    """
    hits = [(pos, score) for pos, score in pssm.search(query, threshold=threshold)]
    if both_strands:
        rc_hits = [
            (pos, score)
            for pos, score in pssm.reverse_complement().search(query, threshold=threshold)
        ]
        hits.extend(rc_hits)
    return sorted(hits, key=lambda t: t[0])


# -----------------------------------------------------------------------------
# WebLogo util (optional, requires *weblogo* binary)
# -----------------------------------------------------------------------------

def save_weblogo(m: motifs.Motif, outfile: Path, fmt: str = "png") -> None:
    """Generate a WebLogo for *m* if *weblogo* is installed."""
    try:
        m.weblogo(outfile.as_posix(), format=fmt)
        logger.info("Saved WebLogo → %s", outfile)
    except (RuntimeError, OSError) as exc:
        logger.warning("Could not create WebLogo: %s", exc)


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Biopython motif utility demo")
    p.add_argument("--search-seq", help="Sequence to scan", type=str)
    p.add_argument("--threshold", help="PSSM score threshold (log2)", type=float, default=3.0)
    p.add_argument("--logo", help="Path to save WebLogo (optional)")
    return p.parse_args()


def main() -> None:
    args = _parse_args()

    # Build example motifs (replace with your own data)
    motif_a = build_motif(EXAMPLE_INSTANCES_A, "Motif-A")
    motif_b = build_motif(EXAMPLE_INSTANCES_B, "Motif-B")

    logger.info("\n%s", motif_summary(motif_a))

    # Build PWM / PSSM using default settings
    pwm, pssm = build_pwm_pssm(motif_a)
    logger.info("Max score %.2f | Min score %.2f", pssm.max, pssm.min)

    # Search the provided sequence (both strands)
    if args.search_seq:
        query_seq = Seq(args.search_seq.upper())
        hits = search_pssm(pssm, query_seq, threshold=args.threshold)
        for pos, score in hits:
            logger.info("Hit at %d (score %.2f)", pos, score)
        if not hits:
            logger.info("No hits ≥ %.2f found", args.threshold)

    # Optional WebLogo
    if args.logo:
        save_weblogo(motif_a, Path(args.logo))

    # Dump combined motifs in JASPAR format
    jaspar_out = motifs.write([motif_a, motif_b], "jaspar")
    logger.info("\nJASPAR output (2 motifs) →\n%s", jaspar_out)


if __name__ == "__main__":
    main()
