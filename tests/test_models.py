import pytest

from enztra.models import Candidate, ReferenceThresholds


@pytest.mark.parametrize("value", [0.0, -1.0, float("inf"), float("nan")])
def test_reference_rejects_invalid_kinetics(value: float) -> None:
    with pytest.raises(ValueError):
        ReferenceThresholds("enzyme", "substrate", kcat_s=value, km_mm=1.0)


def test_candidate_rejects_empty_sequence() -> None:
    with pytest.raises(ValueError):
        Candidate("design", "", kcat_s=1.0, km_mm=1.0)

