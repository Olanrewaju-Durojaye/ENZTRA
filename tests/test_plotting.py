from xml.etree import ElementTree

from enztra.models import Candidate, ReferenceThresholds
from enztra.plotting import display_design_label, write_kinetic_selection_plot
from enztra.selection import evaluate_candidates


def test_writes_standalone_kinetic_selection_svg(tmp_path):
    reference = ReferenceThresholds("reference", "IHP", 100.0, 0.2)
    results = evaluate_candidates(
        [
            Candidate("survivor", "ACDE", 120.0, 0.1),
            Candidate("rejected", "FGHI", 90.0, 0.3),
        ],
        reference,
    )
    path = write_kinetic_selection_plot(
        tmp_path / "results" / "kinetic_selection.svg", reference, results
    )

    assert path.is_file()
    ElementTree.parse(path)
    content = path.read_text()
    assert "reference" in content
    assert "survivor" in content
    assert "rejected" in content  # retained in the point's hover title
    assert ">rejected</text>" not in content
    assert ">survivor</text>" in content
    assert "#13795b" in content
    assert "#c55a4b" in content


def test_plot_handles_reference_only(tmp_path):
    reference = ReferenceThresholds("reference", "IHP", 100.0, 0.2)
    path = write_kinetic_selection_plot(
        tmp_path / "kinetic_selection.svg", reference, []
    )

    ElementTree.parse(path)
    assert path.stat().st_size > 1000


def test_long_survivor_id_is_shortened_only_in_visible_label(tmp_path):
    reference = ReferenceThresholds("reference", "IHP", 100.0, 0.2)
    long_id = "run_project_cond0_0-atomized-bb-False_sequence_001"
    results = evaluate_candidates(
        [Candidate(long_id, "ACDE", 120.0, 0.1)], reference
    )
    path = write_kinetic_selection_plot(tmp_path / "plot.svg", reference, results)
    content = path.read_text()

    assert f"<title>{long_id}:" in content
    assert f">{long_id}</text>" not in content
    assert ">bb000 · sequence_001</text>" in content


def test_generated_label_preserves_backbone_and_sequence_numbers():
    first = "run_release-500_cond0_0-atomized-bb-False_seq0007"
    second = "run_release-500_cond0_11-atomized-bb-False_seq0007"

    assert display_design_label(first) == "bb000 · seq0007"
    assert display_design_label(second) == "bb011 · seq0007"
    assert display_design_label(first) != display_design_label(second)


def test_axes_never_display_negative_kinetic_values(tmp_path):
    reference = ReferenceThresholds("reference", "IHP", 10.0, 0.01)
    results = evaluate_candidates(
        [Candidate("design", "ACDE", 20.0, 0.001)], reference
    )
    path = write_kinetic_selection_plot(tmp_path / "plot.svg", reference, results)

    assert ">-" not in path.read_text()
