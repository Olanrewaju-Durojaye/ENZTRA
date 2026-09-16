import csv
import json
from pathlib import Path
from types import SimpleNamespace

from enztra.boltz2 import (
    find_boltz2_ready_jobs,
    prepare_boltz2_execution,
    refresh_boltz2_reports,
    run_boltz2,
)
from enztra.config import EnztraConfig


def _kinetic_job(tmp_path):
    job = tmp_path / "jobs/kinetic-test"
    job.mkdir(parents=True)
    (job / "manifest.json").write_text(json.dumps({
        "substrate_name": "IHP", "substrate_smiles": "C[C@H](O)P(=O)(O)O"
    }))
    (job / "summary.json").write_text(json.dumps({
        "survivor_count": 2,
        "survivor_ids": ["design/1", "design-2"],
        "boltz2_queue": ["design/1", "design-2"],
    }))
    with (job / "selection.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=(
            "design_id", "sequence", "passes_strict_gate"
        ))
        writer.writeheader()
        writer.writerows([
            {"design_id": "design/1", "sequence": "ACDE", "passes_strict_gate": "True"},
            {"design_id": "design-2", "sequence": "FGHI", "passes_strict_gate": "True"},
            {"design_id": "rejected", "sequence": "KLMN", "passes_strict_gate": "False"},
        ])
    config = EnztraConfig(
        "conda", "apptainer", tmp_path, tmp_path, tmp_path, tmp_path,
        tmp_path, tmp_path / "cache"
    )
    return job, config


def test_plan_writes_only_strict_survivor_yaml(tmp_path):
    job, config = _kinetic_job(tmp_path)
    plan = prepare_boltz2_execution(job, config)
    assert plan["survivor_ids"] == ["design/1", "design-2"]
    assert plan["affinity_requested"] is False
    assert "--use_msa_server" not in plan["command"]
    files = sorted((job / "boltz2/inputs").glob("*/*.yaml"))
    assert [path.name for path in files] == ["design-2.yaml", "design_1.yaml"]
    content = files[0].read_text()
    assert "msa: empty" in content
    assert "properties:" not in content
    assert "smiles:" in content


def test_server_mode_omits_empty_msa_and_adds_flag(tmp_path):
    job, config = _kinetic_job(tmp_path)
    plan = prepare_boltz2_execution(job, config, "server", 3)
    assert "--use_msa_server" in plan["command"]
    assert plan["diffusion_samples"] == 3
    assert all("msa: empty" not in path.read_text() for path in (
        job / "boltz2/inputs"
    ).glob("*/*.yaml"))


def test_ready_jobs_require_nonempty_strict_queue(tmp_path):
    job, _ = _kinetic_job(tmp_path)
    ready = find_boltz2_ready_jobs(tmp_path / "jobs")
    assert [item["job_dir"] for item in ready] == [job]
    summary = json.loads((job / "summary.json").read_text())
    summary["boltz2_queue"] = []
    (job / "summary.json").write_text(json.dumps(summary))
    assert find_boltz2_ready_jobs(tmp_path / "jobs") == []


def test_successful_run_collects_every_expected_model(tmp_path):
    job, config = _kinetic_job(tmp_path)

    def runner(command, **kwargs):
        results_dir = Path(command[command.index("--out_dir") + 1])
        for stem in ("design_1", "design-2"):
            target = results_dir / "run/predictions" / stem
            target.mkdir(parents=True, exist_ok=True)
            (target / f"{stem}_model_0.cif").write_text("data_model\n")
            (target / f"confidence_{stem}_model_0.json").write_text(json.dumps({
                "confidence_score": 0.8, "ptm": 0.7, "iptm": 0.6,
                "ligand_iptm": 0.65, "protein_iptm": 0.0,
                "complex_plddt": 0.75, "complex_iplddt": 0.7,
                "complex_pde": 2.0, "complex_ipde": 3.0,
            }))
        return SimpleNamespace(returncode=0)

    assert run_boltz2(job, config, runner=runner) == 0
    status = json.loads((job / "boltz2/status.json").read_text())
    assert status["stage"] == "boltz2_complete"
    with (job / "boltz2/confidence.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert [row["design_id"] for row in rows] == ["design/1", "design-2"]
    summary = json.loads((job / "boltz2/summary.json").read_text())
    assert summary["affinity_requested"] is False
    assert len(summary["best_models"]) == 2
    assert (job / "boltz2/results/reports/best_models.csv").is_file()
    assert (job / "boltz2/results/reports/boltz2_confidence.svg").is_file()


def test_multiple_samples_are_ranked_and_reports_can_be_refreshed(tmp_path):
    job, config = _kinetic_job(tmp_path)

    def runner(command, **kwargs):
        results_dir = Path(command[command.index("--out_dir") + 1])
        for stem_index, stem in enumerate(("design_1", "design-2")):
            target = results_dir / "run/predictions" / stem
            target.mkdir(parents=True, exist_ok=True)
            for model_index, score in enumerate((0.55 + stem_index * 0.1, 0.85 - stem_index * 0.1)):
                (target / f"{stem}_model_{model_index}.cif").write_text("data_model\n")
                (target / f"confidence_{stem}_model_{model_index}.json").write_text(json.dumps({
                    "confidence_score": score, "ptm": 0.7, "iptm": 0.6,
                    "ligand_iptm": score - 0.05, "protein_iptm": 0.0,
                    "complex_plddt": score - 0.1, "complex_iplddt": 0.7,
                    "complex_pde": 2.0, "complex_ipde": 3.0,
                }))
        return SimpleNamespace(returncode=0)

    assert run_boltz2(job, config, diffusion_samples=2, runner=runner) == 0
    with (job / "boltz2/results/reports/best_models.csv").open(newline="") as handle:
        best = list(csv.DictReader(handle))
    assert [row["design_id"] for row in best] == ["design/1", "design-2"]
    assert [row["model_index"] for row in best] == ["1", "1"]
    assert [row["structural_rank"] for row in best] == ["1", "2"]

    best_csv, plot_path = refresh_boltz2_reports(job)
    assert best_csv.is_file()
    assert plot_path.is_file()


def test_missing_confidence_output_marks_run_failed(tmp_path):
    job, config = _kinetic_job(tmp_path)

    def runner(command, **kwargs):
        return SimpleNamespace(returncode=0)

    assert run_boltz2(job, config, runner=runner) == 2
    status = json.loads((job / "boltz2/status.json").read_text())
    assert status["stage"] == "boltz2_failed"
    assert "expected one confidence_" in status["validation_error"]
