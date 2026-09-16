"""Typed records shared by model adapters, the API, and the user interface."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite


def _require_positive_finite(name: str, value: float) -> None:
    if not isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be a positive finite number; received {value!r}")


@dataclass(frozen=True, slots=True)
class ReferenceThresholds:
    """Predicted kinetic baseline for one reference enzyme–substrate pair."""

    enzyme_id: str
    substrate_id: str
    kcat_s: float
    km_mm: float

    def __post_init__(self) -> None:
        if not self.enzyme_id.strip():
            raise ValueError("enzyme_id must not be empty")
        if not self.substrate_id.strip():
            raise ValueError("substrate_id must not be empty")
        _require_positive_finite("kcat_s", self.kcat_s)
        _require_positive_finite("km_mm", self.km_mm)

    @property
    def catalytic_efficiency(self) -> float:
        """Return kcat/Km in s^-1 mM^-1."""

        return self.kcat_s / self.km_mm


@dataclass(frozen=True, slots=True)
class Candidate:
    """Kinetic predictions for one designed enzyme sequence."""

    design_id: str
    sequence: str
    kcat_s: float
    km_mm: float
    kcat_uncertainty: float | None = None
    km_uncertainty: float | None = None

    def __post_init__(self) -> None:
        if not self.design_id.strip():
            raise ValueError("design_id must not be empty")
        if not self.sequence.strip():
            raise ValueError("sequence must not be empty")
        _require_positive_finite("kcat_s", self.kcat_s)
        _require_positive_finite("km_mm", self.km_mm)
        for name, value in (
            ("kcat_uncertainty", self.kcat_uncertainty),
            ("km_uncertainty", self.km_uncertainty),
        ):
            if value is not None and (not isfinite(value) or value < 0):
                raise ValueError(f"{name} must be non-negative and finite")

    @property
    def catalytic_efficiency(self) -> float:
        return self.kcat_s / self.km_mm


@dataclass(frozen=True, slots=True)
class SelectionResult:
    """Decision and derived metrics for one candidate."""

    candidate: Candidate
    passes_kcat: bool
    passes_km: bool
    passes_strict_gate: bool
    kcat_fold_change: float
    km_fold_change: float
    efficiency_fold_change: float
    rank: int | None = None

    @property
    def decision(self) -> str:
        return "SURVIVOR" if self.passes_strict_gate else "REJECTED"

