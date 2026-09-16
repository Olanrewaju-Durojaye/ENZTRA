"""Reference-relative kinetic selection for ENZTRA.

The strict gate intentionally uses logical AND. A high catalytic efficiency
cannot rescue a candidate that fails either the kcat or Km threshold.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Iterable

from .models import Candidate, ReferenceThresholds, SelectionResult


def evaluate_candidate(
    candidate: Candidate,
    reference: ReferenceThresholds,
) -> SelectionResult:
    """Evaluate one candidate against strict reference thresholds.

    Equality does not pass: improvement must be strictly higher for kcat and
    strictly lower for Km.
    """

    passes_kcat = candidate.kcat_s > reference.kcat_s
    passes_km = candidate.km_mm < reference.km_mm
    return SelectionResult(
        candidate=candidate,
        passes_kcat=passes_kcat,
        passes_km=passes_km,
        passes_strict_gate=passes_kcat and passes_km,
        kcat_fold_change=candidate.kcat_s / reference.kcat_s,
        km_fold_change=candidate.km_mm / reference.km_mm,
        efficiency_fold_change=(
            candidate.catalytic_efficiency / reference.catalytic_efficiency
        ),
    )


def evaluate_candidates(
    candidates: Iterable[Candidate],
    reference: ReferenceThresholds,
) -> list[SelectionResult]:
    """Evaluate all candidates and rank only strict survivors by efficiency."""

    results = [evaluate_candidate(candidate, reference) for candidate in candidates]
    survivors = sorted(
        (result for result in results if result.passes_strict_gate),
        key=lambda result: (
            result.efficiency_fold_change,
            result.kcat_fold_change,
            -result.km_fold_change,
        ),
        reverse=True,
    )
    ranks = {result.candidate.design_id: rank for rank, result in enumerate(survivors, 1)}
    return [
        replace(result, rank=ranks.get(result.candidate.design_id))
        for result in results
    ]

