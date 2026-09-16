import csv
import json
from collections import namedtuple

from enztra.demonstration import (
    generate_release_report,
    prepare_release_demonstration,
)


def _project(tmp_path, backbones=50, sequences_per_backbone=10):
    job = tmp_path / "jobs/release-500"
    job.mkdir(parents=True)
    total = backbones * sequences_per_backbone
    (job / "design_request.json").write_text(json.dumps({
        "job_id": "release-500",
        "backbone_count": backbones,
        "sequences_per_backbone": sequences_per_backbone,
        "total_designs": total,
        "design_length": "150-180",
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


def _ready(monkeypatch, free_gib=100):
    monkeypatch.setattr("enztra.demonstration.run_doctor", lambda config: {
        "status": "ready", "checks": [{"name": "GPU", "ok": True}],
    })
    usage = namedtuple("usage", "total used free")
    monkeypatch.setattr(
        "enztra.demonstration.shutil.disk_usage",
        lambda path: usage(200 * 1024**3, 0, free_gib * 1024**3),
    )


def test_exact_release_plan_passes_preflight(tmp_path, monkeypatch):
    job, config = _project(tmp_path)
    _ready(monkeypatch)
    report = prepare_release_demonstration(job, config)
    assert report["ready"] is True
    assert report["plan"]["total_sequences"] == 500
    assert (job / "release/preflight.json").is_file()


def test_nonrelease_plan_is_blocked(tmp_path, monkeypatch):
    job, config = _project(tmp_path, backbones=10, sequences_per_backbone=5)
    _ready(monkeypatch)
    report = prepare_release_demonstration(job, config)
    assert report["ready"] is False
    assert any("exactly 50 backbones" in item for item in report["blockers"])


def test_insufficient_disk_is_blocked(tmp_path, monkeypatch):
    job, config = _project(tmp_path)
    _ready(monkeypatch, free_gib=10)
    report = prepare_release_demonstration(job, config)
    assert report["ready"] is False
    assert any("insufficient disk space" in item for item in report["blockers"])


def test_complete_500_sequence_run_qualifies(tmp_path):
    job, _ = _project(tmp_path)
    output = job / "outputs/rfdiffusion2"
    output.mkdir(parents=True)
    for index in range(50):
        (output / f"backbone_{index}.pdb").write_text("MODEL\n")
        (output / f"backbone_{index}.trb").write_bytes(b"trb")
    fasta = job / "outputs/ligandmpnn/designs.fasta"
    fasta.parent.mkdir(parents=True)
    fasta.write_text("".join(f">design_{index}\nACDE\n" for index in range(500)))
    with (job / "selection.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["design_id"])
        writer.writeheader()
        for index in range(500):
            writer.writerow({"design_id": f"design_{index}"})
    (job / "summary.json").write_text(json.dumps({
        "survivor_count": 0, "boltz2_queue": [],
    }))
    publication = job / "publication"
    publication.mkdir()
    (publication / "publication_summary.json").write_text(json.dumps({
        "candidate_count": 500,
    }))
    (job / "workflow_state.json").write_text(json.dumps({
        "workflow_status": "completed", "stages": {},
    }))
    report = generate_release_report(job)
    assert report["qualification_status"] == "passed"
    assert report["github_release_ready"] is True
    assert (job / "release/qualification_report.md").is_file()
