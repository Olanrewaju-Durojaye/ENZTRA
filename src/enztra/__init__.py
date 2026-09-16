"""ENZTRA core package."""

from .models import Candidate, ReferenceThresholds, SelectionResult
from .selection import evaluate_candidate, evaluate_candidates

__all__ = [
    "Candidate",
    "ReferenceThresholds",
    "SelectionResult",
    "evaluate_candidate",
    "evaluate_candidates",
]

__version__ = "1.0.1"
