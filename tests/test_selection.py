import pytest

from enztra.models import Candidate, ReferenceThresholds
from enztra.selection import evaluate_candidate, evaluate_candidates


REFERENCE = ReferenceThresholds(
    enzyme_id="1DKP_A",
    substrate_id="IHP",
    kcat_s=314.0739,
    km_mm=0.181112,
)


def candidate(design_id: str, kcat: float, km: float) -> Candidate:
    return Candidate(design_id=design_id, sequence="ACDEFG", kcat_s=kcat, km_mm=km)


@pytest.mark.parametrize(
    ("kcat", "km", "passes"),
    [
        (400.0, 0.10, True),
        (400.0, 0.30, False),
        (200.0, 0.10, False),
        (200.0, 0.30, False),
        (REFERENCE.kcat_s, 0.10, False),
        (400.0, REFERENCE.km_mm, False),
    ],
)
def test_strict_and_gate(kcat: float, km: float, passes: bool) -> None:
    result = evaluate_candidate(candidate("d", kcat, km), REFERENCE)
    assert result.passes_strict_gate is passes


def test_efficiency_does_not_rescue_failed_kcat() -> None:
    design = candidate("high_efficiency_but_low_kcat", 300.0, 0.01)
    result = evaluate_candidate(design, REFERENCE)
    assert result.efficiency_fold_change > 1
    assert not result.passes_kcat
    assert not result.passes_strict_gate


def test_only_survivors_receive_ranks() -> None:
    results = evaluate_candidates(
        [
            candidate("survivor_2", 400.0, 0.10),
            candidate("rejected", 500.0, 0.20),
            candidate("survivor_1", 500.0, 0.05),
        ],
        REFERENCE,
    )
    by_id = {result.candidate.design_id: result for result in results}
    assert by_id["survivor_1"].rank == 1
    assert by_id["survivor_2"].rank == 2
    assert by_id["rejected"].rank is None


def test_reference_efficiency() -> None:
    assert REFERENCE.catalytic_efficiency == pytest.approx(1734.141856972481)

