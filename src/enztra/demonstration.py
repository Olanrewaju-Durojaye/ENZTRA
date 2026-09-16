"""Release qualification for the native 500-sequence demonstration."""

from __future__ import annotations

import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from . import __version__
from .config import EnztraConfig
from .doctor import run_doctor
from .pilot import generate_pilot_report
from .workflow import run_workflow

EXPECTED_BACKBONES = 50
EXPECTED_SEQUENCES_PER_BACKBONE = 10
EXPECTED_SEQUENCES = 500
RESERVED_FREE_BYTES = 5 * 1024**3


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"required demonstration file is missing: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def _tree_size(path: Path) -> int:
    if not path.is_dir():
        return 0
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def _estimated_bytes(calibration_job: Path | None) -> tuple[int, str]:
    if calibration_job is None or not calibration_job.is_dir():
        return 50 * 1024**3, "conservative 50 GiB estimate; no calibration job supplied"
    request_path = calibration_job / "design_request.json"
    if not request_path.is_file():
        return 50 * 1024**3, "conservative 50 GiB estimate; calibration request missing"
    request = _read_json(request_path)
    backbones = max(1, int(request.get("backbone_count", 1)))
    sequences = max(1, int(request.get("total_designs", 1)))
    rfd_size = _tree_size(calibration_job / "outputs/rfdiffusion2")
    total_size = _tree_size(calibration_job)
    other_size = max(0, total_size - rfd_size)
    estimate = (
        rfd_size * EXPECTED_BACKBONES / backbones
        + other_size * EXPECTED_SEQUENCES / sequences
    )
    # Preserve a 25% allowance for variable trajectories, logs, and survivors.
    return max(int(estimate * 1.25), 5 * 1024**3), (
        f"scaled from {calibration_job.name} with a 25% safety allowance"
    )


def prepare_release_demonstration(
    job_dir: Path, config_path: Path, calibration_job: Path | None = None
) -> dict[str, Any]:
    """Enforce the release plan and estimate local storage before execution."""

    job_dir = job_dir.expanduser().resolve()
    config_path = config_path.expanduser().resolve()
    calibration = calibration_job.expanduser().resolve() if calibration_job else None
    request = _read_json(job_dir / "design_request.json")
    status = _read_json(job_dir / "status.json")
    if not status.get("execution_ready"):
        raise ValueError("release demonstration requires RFdiffusion2-ready PDB input")
    backbones = int(request["backbone_count"])
    per_backbone = int(request["sequences_per_backbone"])
    total = int(request["total_designs"])
    plan_ok = (
        backbones == EXPECTED_BACKBONES
        and per_backbone == EXPECTED_SEQUENCES_PER_BACKBONE
        and total == EXPECTED_SEQUENCES
    )
    estimate, estimate_method = _estimated_bytes(calibration)
    free = shutil.disk_usage(job_dir).free
    required = estimate + RESERVED_FREE_BYTES
    doctor = run_doctor(EnztraConfig.from_json(config_path))
    failed_checks = [check["name"] for check in doctor["checks"] if not check["ok"]]
    blockers: list[str] = []
    if not plan_ok:
        blockers.append(
            "release plan must be exactly 50 backbones × 10 sequences = 500"
        )
    if free < required:
        blockers.append(
            f"insufficient disk space: {free / 1024**3:.2f} GiB free; "
            f"{required / 1024**3:.2f} GiB required including reserve"
        )
    if doctor["status"] != "ready" or failed_checks:
        blockers.append("one or more required local installations failed")
    preflight = {
        "schema_version": 1,
        "enztra_version": __version__,
        "created_at": _now(),
        "job_id": request["job_id"],
        "profile": "v1_release_demonstration",
        "plan": {
            "backbones": backbones,
            "sequences_per_backbone": per_backbone,
            "total_sequences": total,
        },
        "required_plan": {
            "backbones": EXPECTED_BACKBONES,
            "sequences_per_backbone": EXPECTED_SEQUENCES_PER_BACKBONE,
            "total_sequences": EXPECTED_SEQUENCES,
        },
        "calibration_job": str(calibration) if calibration else "",
        "estimated_output_gib": round(estimate / 1024**3, 2),
        "estimation_method": estimate_method,
        "reserved_free_gib": round(RESERVED_FREE_BYTES / 1024**3, 2),
        "required_free_gib": round(required / 1024**3, 2),
        "available_free_gib": round(free / 1024**3, 2),
        "installation_status": doctor["status"],
        "failed_installation_checks": failed_checks,
        "blockers": blockers,
        "ready": not blockers,
    }
    destination = job_dir / "release/preflight.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(preflight, indent=2) + "\n", encoding="utf-8")
    return preflight


def generate_release_report(job_dir: Path) -> dict[str, Any]:
    """Create the v1 qualification record from the standard integration report."""

    job_dir = job_dir.expanduser().resolve()
    integration = generate_pilot_report(job_dir)
    request = _read_json(job_dir / "design_request.json")
    exact_plan = (
        int(request["backbone_count"]) == EXPECTED_BACKBONES
        and int(request["sequences_per_backbone"]) == EXPECTED_SEQUENCES_PER_BACKBONE
        and int(request["total_designs"]) == EXPECTED_SEQUENCES
    )
    qualified = integration["pilot_status"] == "passed" and exact_plan
    report = {
        "schema_version": 1,
        "enztra_version": __version__,
        "generated_at": _now(),
        "job_id": request["job_id"],
        "qualification_status": "passed" if qualified else "incomplete_or_failed",
        "exact_500_sequence_plan": exact_plan,
        "integration_report": integration,
        "github_release_ready": qualified,
    }
    output = job_dir / "release/qualification_report.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    markdown = [
        f"# ENZTRA v1 release qualification: {request['job_id']}", "",
        f"- Qualification: {report['qualification_status']}",
        f"- Exact 500-sequence plan: {'yes' if exact_plan else 'no'}",
        f"- GitHub release ready: {'yes' if qualified else 'no'}", "",
        "## Observed outputs", "",
    ]
    for key, value in integration["observed"].items():
        markdown.append(f"- {key.replace('_', ' ').title()}: {value}")
    (job_dir / "release/qualification_report.md").write_text(
        "\n".join(markdown) + "\n", encoding="utf-8"
    )
    return report


def run_release_demonstration(
    job_dir: Path, config_path: Path, *, calibration_job: Path | None = None,
    msa_mode: str = "single", diffusion_samples: int = 1, seed: int = 42,
) -> int:
    """Preflight, run/resume, and qualify the 500-sequence demonstration."""

    preflight = prepare_release_demonstration(job_dir, config_path, calibration_job)
    if not preflight["ready"]:
        return 2
    try:
        returncode = run_workflow(
            job_dir, config_path, msa_mode=msa_mode,
            diffusion_samples=diffusion_samples, seed=seed,
        )
    except BaseException:
        generate_release_report(job_dir)
        raise
    generate_release_report(job_dir)
    return returncode
