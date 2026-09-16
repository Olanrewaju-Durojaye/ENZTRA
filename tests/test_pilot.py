import csv
import json

from enztra.pilot import generate_pilot_report, prepare_pilot


def _project(tmp_path, backbones=10, sequences_per_backbone=5):
    job = tmp_path / "jobs/pilot-50"
    job.mkdir(parents=True)
    total = backbones * sequences_per_backbone
    (job / "design_request.json").write_text(json.dumps({
        "job_id": "pilot-50", "backbone_count": backbones,
        "sequences_per_backbone": sequences_per_backbone,
        "total_designs": total, "design_length": "150-180",
    }))
    (job / "status.json").write_text(json.dumps({"execution_ready": True}))
    config = tmp_path / "config.json"
    config.write_text(json.dumps({
        "rfdiffusion2_root": str(tmp_path / "rf"),
        "dlkcat_root": str(tmp_path / "dlkcat"),
        "catpred_root": str(tmp_path / "catpred"),
        "catpred_data_root": str(tmp_path / "catpred-data"),
    }))
    return job, config


def test_recommended_pilot_preflight(tmp_path, monkeypatch):
    job, config = _project(tmp_path)
    monkeypatch.setattr("enztra.pilot.run_doctor", lambda config: {
        "status": "ready", "checks": [{"name": "GPU", "ok": True}],
    })
    result = prepare_pilot(job, config)
    assert result["ready"] is True
    assert result["planned_sequences"] == 50
    assert not any("pilot_size" in warning for warning in result["warnings"])
    assert (job / "pilot/preflight.json").is_file()


def test_small_pilot_receives_scale_warnings(tmp_path, monkeypatch):
    job, config = _project(tmp_path, backbones=5, sequences_per_backbone=5)
    monkeypatch.setattr("enztra.pilot.run_doctor", lambda config: {
        "status": "ready", "checks": [],
    })
    result = prepare_pilot(job, config)
    assert any("pilot_size_outside" in warning for warning in result["warnings"])
    assert any("limited_backbone_diversity" in warning for warning in result["warnings"])


def test_pilot_report_verifies_expected_counts(tmp_path):
    job, _ = _project(tmp_path, backbones=2, sequences_per_backbone=2)
    output = job / "outputs/rfdiffusion2"
    output.mkdir(parents=True)
    for name in ("a", "b"):
        (output / f"{name}.pdb").write_text("MODEL\n")
        (output / f"{name}.trb").write_bytes(b"trb")
    fasta = job / "outputs/ligandmpnn/designs.fasta"
    fasta.parent.mkdir(parents=True)
    fasta.write_text(">a1\nACDE\n>a2\nACDF\n>b1\nACDG\n>b2\nACDH\n")
    with (job / "selection.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["design_id"])
        writer.writeheader()
        for index in range(4):
            writer.writerow({"design_id": f"d{index}"})
    (job / "summary.json").write_text(json.dumps({
        "survivor_count": 0, "boltz2_queue": [],
    }))
    publication = job / "publication"
    publication.mkdir()
    (publication / "publication_summary.json").write_text(json.dumps({
        "candidate_count": 4,
    }))
    (job / "workflow_state.json").write_text(json.dumps({
        "workflow_status": "completed", "stages": {},
    }))
    report = generate_pilot_report(job)
    assert report["pilot_status"] == "passed"
    assert not report["failed_checks"]
    assert (job / "pilot/integration_report.md").is_file()


def test_pilot_report_exposes_partial_outputs(tmp_path):
    job, _ = _project(tmp_path, backbones=2, sequences_per_backbone=2)
    (job / "workflow_state.json").write_text(json.dumps({
        "workflow_status": "failed", "stages": {},
    }))
    report = generate_pilot_report(job)
    assert report["pilot_status"] == "incomplete_or_failed"
    assert "backbone_count" in report["failed_checks"]
    assert "sequence_count" in report["failed_checks"]
