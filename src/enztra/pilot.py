"""Controlled pilot preflight, execution, and integration reporting."""

from __future__ import annotations

import csv
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from . import __version__
from .config import EnztraConfig
from .doctor import run_doctor
from .fasta import read_fasta
from .workflow import STAGES, run_workflow


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _read_json(path: Path, required: bool = True) -> dict[str, Any]:
    if not path.is_file():
        if required:
            raise FileNotFoundError(f"required pilot file is missing: {path}")
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def _csv_count(path: Path) -> int:
    if not path.is_file():
        return 0
    with path.open(newline="", encoding="utf-8") as handle:
        return sum(1 for _ in csv.DictReader(handle))


def prepare_pilot(job_dir: Path, config_path: Path) -> dict[str, Any]:
    """Validate a prepared job and record a controlled-pilot preflight."""

    job_dir = job_dir.expanduser().resolve()
    config_path = config_path.expanduser().resolve()
    request = _read_json(job_dir / "design_request.json")
    status = _read_json(job_dir / "status.json")
    if not status.get("execution_ready"):
        raise ValueError("pilot requires a PDB project that passed RFdiffusion2 preflight")
    backbones = int(request["backbone_count"])
    per_backbone = int(request["sequences_per_backbone"])
    total = backbones * per_backbone
    warnings: list[str] = []
    if not 50 <= total <= 100:
        warnings.append(
            f"pilot_size_outside_recommended_range: planned {total}; recommended 50-100"
        )
    if backbones < 10:
        warnings.append(
            "limited_backbone_diversity: fewer than 10 RFdiffusion2 backbones"
        )
    usage = shutil.disk_usage(job_dir)
    free_gib = usage.free / 1024**3
    if free_gib < 20:
        warnings.append(
            f"low_free_disk_space: {free_gib:.1f} GiB available; at least 20 GiB recommended"
        )
    doctor = run_doctor(EnztraConfig.from_json(config_path))
    failed_checks = [check["name"] for check in doctor["checks"] if not check["ok"]]
    preflight = {
        "schema_version": 1,
        "enztra_version": __version__,
        "created_at": _now(),
        "job_id": request["job_id"],
        "profile": "controlled_native_pilot",
        "planned_backbones": backbones,
        "sequences_per_backbone": per_backbone,
        "planned_sequences": total,
        "designed_length": request["design_length"],
        "free_disk_gib": round(free_gib, 2),
        "installation_status": doctor["status"],
        "failed_installation_checks": failed_checks,
        "warnings": warnings,
        "ready": doctor["status"] == "ready" and not failed_checks,
    }
    _write_json(job_dir / "pilot/preflight.json", preflight)
    return preflight


def _stage_durations(state: dict[str, Any]) -> dict[str, float | None]:
    durations: dict[str, float | None] = {}
    for name in STAGES:
        stage = state.get("stages", {}).get(name, {})
        try:
            start = datetime.fromisoformat(stage["started_at"])
            finish = datetime.fromisoformat(stage["finished_at"])
            durations[name] = round((finish - start).total_seconds(), 3)
        except (KeyError, TypeError, ValueError):
            durations[name] = None
    return durations


def generate_pilot_report(job_dir: Path) -> dict[str, Any]:
    """Compare planned and observed outputs after a pilot attempt."""

    job_dir = job_dir.expanduser().resolve()
    request = _read_json(job_dir / "design_request.json")
    state = _read_json(job_dir / "workflow_state.json", required=False)
    planned_backbones = int(request["backbone_count"])
    planned_sequences = int(request["total_designs"])
    rfd_dir = job_dir / "outputs/rfdiffusion2"
    pdb_stems = {path.stem for path in rfd_dir.glob("*.pdb")}
    trb_stems = {path.stem for path in rfd_dir.glob("*.trb")}
    completed_backbones = len(pdb_stems & trb_stems)
    fasta_path = job_dir / "outputs/ligandmpnn/designs.fasta"
    sequence_count = len(read_fasta(fasta_path)) if fasta_path.is_file() else 0
    kinetic_count = _csv_count(job_dir / "selection.csv")
    summary = _read_json(job_dir / "summary.json", required=False)
    survivor_count = int(summary.get("survivor_count", 0))
    structural_count = _csv_count(
        job_dir / "boltz2/results/reports/best_models.csv"
    )
    publication = _read_json(
        job_dir / "publication/publication_summary.json", required=False
    )
    checks = [
        {"name": "backbone_count", "expected": planned_backbones,
         "observed": completed_backbones,
         "ok": completed_backbones == planned_backbones},
        {"name": "sequence_count", "expected": planned_sequences,
         "observed": sequence_count, "ok": sequence_count == planned_sequences},
        {"name": "kinetic_candidate_count", "expected": planned_sequences,
         "observed": kinetic_count, "ok": kinetic_count == planned_sequences},
        {"name": "structural_survivor_count", "expected": survivor_count,
         "observed": structural_count,
         "ok": survivor_count == 0 or structural_count == survivor_count},
        {"name": "publication_candidate_count", "expected": kinetic_count,
         "observed": int(publication.get("candidate_count", 0)),
         "ok": bool(publication) and int(publication.get("candidate_count", 0)) == kinetic_count},
    ]
    failures = [check["name"] for check in checks if not check["ok"]]
    report = {
        "schema_version": 1,
        "enztra_version": __version__,
        "generated_at": _now(),
        "job_id": request["job_id"],
        "workflow_status": state.get("workflow_status", "not_started"),
        "planned": {"backbones": planned_backbones, "sequences": planned_sequences},
        "observed": {
            "completed_backbones": completed_backbones,
            "sequences": sequence_count,
            "kinetic_candidates": kinetic_count,
            "strict_survivors": survivor_count,
            "structurally_ranked": structural_count,
        },
        "stage_duration_seconds": _stage_durations(state),
        "integration_checks": checks,
        "failed_checks": failures,
        "pilot_status": (
            "passed" if state.get("workflow_status") == "completed" and not failures
            else "incomplete_or_failed"
        ),
    }
    _write_json(job_dir / "pilot/integration_report.json", report)
    lines = [
        f"# ENZTRA controlled pilot: {request['job_id']}", "",
        f"- Pilot status: {report['pilot_status']}",
        f"- Workflow status: {report['workflow_status']}",
        f"- Backbones: {completed_backbones}/{planned_backbones}",
        f"- Protein sequences: {sequence_count}/{planned_sequences}",
        f"- Kinetic candidates: {kinetic_count}/{planned_sequences}",
        f"- Strict survivors: {survivor_count}",
        f"- Structurally ranked: {structural_count}", "",
        "## Integration checks", "",
    ]
    lines.extend(
        f"- [{'x' if check['ok'] else ' '}] {check['name']}: "
        f"expected {check['expected']}, observed {check['observed']}"
        for check in checks
    )
    (job_dir / "pilot/integration_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    return report


def run_controlled_pilot(
    job_dir: Path, config_path: Path, *, msa_mode: str = "single",
    diffusion_samples: int = 1, seed: int = 42,
) -> int:
    """Preflight, run/resume, and report a controlled native pilot."""

    preflight = prepare_pilot(job_dir, config_path)
    if not preflight["ready"]:
        generate_pilot_report(job_dir)
        return 2
    try:
        returncode = run_workflow(
            job_dir, config_path, msa_mode=msa_mode,
            diffusion_samples=diffusion_samples, seed=seed,
        )
    except BaseException:
        generate_pilot_report(job_dir)
        raise
    generate_pilot_report(job_dir)
    return returncode
