"""Resumable, stage-aware execution of the complete ENZTRA workflow."""

from __future__ import annotations

import json
from argparse import Namespace
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

from .boltz2 import run_boltz2
from .cli import _kinetics
from .config import EnztraConfig
from .design import pdb_chain_sequence
from .execution import run_rfdiffusion2
from .fasta import read_fasta
from .ligandmpnn import run_ligandmpnn

STAGES = ("rfdiffusion2", "ligandmpnn", "kinetics", "boltz2", "publication")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def _initial_state(job_id: str) -> dict:
    return {
        "schema_version": 1,
        "job_id": job_id,
        "workflow_status": "pending",
        "created_at": _now(),
        "updated_at": _now(),
        "stages": {
            name: {"status": "pending", "attempts": 0, "message": ""}
            for name in STAGES
        },
    }


def load_workflow_state(job_dir: Path) -> dict:
    """Load, validate, or initialize the persistent workflow state."""

    job_dir = job_dir.expanduser().resolve()
    request = _read_json(job_dir / "design_request.json")
    path = job_dir / "workflow_state.json"
    if path.is_file():
        state = _read_json(path)
        if state.get("schema_version") != 1 or not isinstance(state.get("stages"), dict):
            raise ValueError("workflow_state.json has an unsupported or invalid schema")
        for name in STAGES:
            state["stages"].setdefault(
                name, {"status": "pending", "attempts": 0, "message": ""}
            )
        return state
    state = _initial_state(request["job_id"])
    _write_json(path, state)
    return state


def _stage_outputs_valid(job_dir: Path, stage: str) -> bool:
    """Return whether durable outputs prove that a stage is complete."""

    if stage == "rfdiffusion2":
        status = _read_json(job_dir / "status.json")
        output = job_dir / "outputs/rfdiffusion2"
        return (
            status.get("rfdiffusion2_returncode") == 0
            and status.get("rfdiffusion2_validated") is True
            and (job_dir / "outputs/functional_site_mapping.csv").is_file()
            and bool(list(output.glob("*.pdb")))
            and bool(list(output.glob("*.trb")))
        )
    if stage == "ligandmpnn":
        status = _read_json(job_dir / "status.json")
        fasta = Path(status.get("ligandmpnn_designs_fasta", ""))
        if not fasta.is_absolute():
            fasta = job_dir / fasta
        return (
            status.get("ligandmpnn_returncode") == 0
            and status.get("ligandmpnn_fixed_positions_validated") is True
            and fasta.is_file()
            and bool(read_fasta(fasta))
        )
    if stage == "kinetics":
        required = ("manifest.json", "kinetics.csv", "selection.csv", "summary.json")
        if not all((job_dir / name).is_file() for name in required):
            return False
        summary = _read_json(job_dir / "summary.json")
        return isinstance(summary.get("boltz2_queue"), list)
    if stage == "boltz2":
        status_path = job_dir / "boltz2/status.json"
        if not status_path.is_file():
            return False
        status = _read_json(status_path)
        return (
            status.get("returncode") == 0
            and Path(status.get("best_models_csv", "")).is_file()
            and Path(status.get("confidence_plot", "")).is_file()
        )
    if stage == "publication":
        required = (
            "publication/tables/master_ranking.csv",
            "publication/tables/master_ranking.json",
            "publication/tables/quality_control.csv",
            "publication/figures/candidate_outcomes.svg",
            "publication/publication_summary.json",
            "publication/reproducibility_manifest.json",
            "publication/publication_report.md",
        )
        return all((job_dir / name).is_file() for name in required)
    raise ValueError(f"unknown workflow stage: {stage}")


def reconcile_workflow_state(job_dir: Path, state: dict | None = None) -> dict:
    """Reconcile recorded completion with validated files on disk."""

    job_dir = job_dir.expanduser().resolve()
    state = state or load_workflow_state(job_dir)
    upstream_complete = True
    for name in STAGES:
        record = state["stages"][name]
        if record.get("status") == "skipped":
            upstream_complete = upstream_complete and name == "boltz2"
            continue
        valid = False
        if upstream_complete:
            try:
                valid = _stage_outputs_valid(job_dir, name)
            except (OSError, ValueError, json.JSONDecodeError):
                valid = False
        if valid:
            record["status"] = "completed"
            record["message"] = "Validated existing outputs; stage will not be rerun."
        elif record.get("status") == "completed":
            record["status"] = "failed"
            record["message"] = "Recorded completion is invalid because required outputs are missing or unreadable."
        elif record.get("status") == "running":
            record["status"] = "interrupted"
            record["message"] = "A previous run stopped before this stage completed."
        upstream_complete = upstream_complete and valid
    state["updated_at"] = _now()
    _write_json(job_dir / "workflow_state.json", state)
    return state


def _reference_fasta(job_dir: Path, request: dict) -> Path:
    destination = job_dir / "inputs/workflow_reference.fasta"
    source = job_dir / request["reference_file"]
    if request["reference_format"] == "pdb":
        sequence = pdb_chain_sequence(source.read_text(encoding="utf-8"), request["protein_chain"])
    else:
        records = read_fasta(source)
        if len(records) != 1:
            raise ValueError("the prepared reference FASTA must contain exactly one record")
        sequence = records[0].sequence
    destination.write_text(f">{request['reference_id']}\n{sequence}\n", encoding="utf-8")
    return destination


def _run_kinetics(job_dir: Path, config_path: Path) -> int:
    request = _read_json(job_dir / "design_request.json")
    status = _read_json(job_dir / "status.json")
    designs = Path(status.get("ligandmpnn_designs_fasta", ""))
    if not designs.is_absolute():
        designs = job_dir / designs
    args = Namespace(
        config=config_path,
        reference_fasta=_reference_fasta(job_dir, request),
        designs_fasta=designs,
        substrate_name=request["substrate_name"],
        substrate_smiles=request["substrate_smiles"],
        job_dir=job_dir,
        prepare_only=False,
        reuse_raw=False,
    )
    return int(_kinetics(args))


def _run_stage(
    name: str, job_dir: Path, config: EnztraConfig, config_path: Path,
    msa_mode: str, diffusion_samples: int, seed: int,
) -> int:
    if name == "rfdiffusion2":
        return run_rfdiffusion2(job_dir, config)
    if name == "ligandmpnn":
        return run_ligandmpnn(job_dir, config)
    if name == "kinetics":
        return _run_kinetics(job_dir, config_path)
    if name == "boltz2":
        return run_boltz2(job_dir, config, msa_mode, diffusion_samples, seed)
    if name == "publication":
        from .publication import generate_publication_exports

        generate_publication_exports(job_dir)
        return 0
    raise ValueError(f"unknown workflow stage: {name}")


def run_workflow(
    job_dir: Path,
    config_path: Path,
    *,
    msa_mode: str = "single",
    diffusion_samples: int = 1,
    seed: int = 42,
    stage_runner: Callable[[str], int] | None = None,
) -> int:
    """Run or resume all stages, stopping safely at the first failure."""

    job_dir = job_dir.expanduser().resolve()
    config_path = config_path.expanduser().resolve()
    config = EnztraConfig.from_json(config_path)
    state = reconcile_workflow_state(job_dir)
    state["workflow_status"] = "running"
    state["updated_at"] = _now()
    _write_json(job_dir / "workflow_state.json", state)

    for name in STAGES:
        record = state["stages"][name]
        if record["status"] in {"completed", "skipped"}:
            print(f"[SKIP] {name}: {record['message']}")
            continue
        if name == "boltz2":
            queue = _read_json(job_dir / "summary.json").get("boltz2_queue", [])
            if not queue:
                record.update({
                    "status": "skipped",
                    "message": "No design passed both strict kinetic thresholds.",
                    "finished_at": _now(),
                })
                state["updated_at"] = _now()
                _write_json(job_dir / "workflow_state.json", state)
                print("[SKIP] boltz2: no strict kinetic survivors")
                continue

        record.update({
            "status": "running",
            "attempts": int(record.get("attempts", 0)) + 1,
            "started_at": _now(),
            "message": "",
        })
        state["current_stage"] = name
        state["updated_at"] = _now()
        _write_json(job_dir / "workflow_state.json", state)
        print(f"[RUN] {name} (attempt {record['attempts']})")
        try:
            returncode = (
                int(stage_runner(name)) if stage_runner is not None
                else _run_stage(name, job_dir, config, config_path, msa_mode, diffusion_samples, seed)
            )
            if returncode == 0 and stage_runner is None and not _stage_outputs_valid(job_dir, name):
                raise ValueError("stage returned success but its required outputs did not validate")
        except BaseException as error:
            record.update({
                "status": "interrupted" if isinstance(error, (KeyboardInterrupt, SystemExit)) else "failed",
                "finished_at": _now(),
                "message": f"{type(error).__name__}: {error}",
            })
            state["workflow_status"] = record["status"]
            state["updated_at"] = _now()
            _write_json(job_dir / "workflow_state.json", state)
            if isinstance(error, (KeyboardInterrupt, SystemExit)):
                raise
            print(f"[FAIL] {name}: {record['message']}")
            return 1
        if returncode != 0:
            record.update({
                "status": "failed", "returncode": returncode,
                "finished_at": _now(), "message": f"Stage exited with code {returncode}.",
            })
            state["workflow_status"] = "failed"
            state["updated_at"] = _now()
            _write_json(job_dir / "workflow_state.json", state)
            print(f"[FAIL] {name}: exit code {returncode}")
            return returncode
        record.update({
            "status": "completed", "returncode": 0,
            "finished_at": _now(), "message": "Stage completed and required outputs were validated.",
        })
        state["updated_at"] = _now()
        _write_json(job_dir / "workflow_state.json", state)
        print(f"[PASS] {name}")

    state["workflow_status"] = "completed"
    state["current_stage"] = ""
    state["completed_at"] = _now()
    state["updated_at"] = _now()
    _write_json(job_dir / "workflow_state.json", state)
    return 0


def workflow_summary(job_dir: Path) -> list[tuple[str, str, int, str]]:
    """Return compact stage status rows for terminal presentation."""

    state = reconcile_workflow_state(job_dir)
    return [
        (name, str(state["stages"][name]["status"]),
         int(state["stages"][name].get("attempts", 0)),
         str(state["stages"][name].get("message", "")))
        for name in STAGES
    ]
