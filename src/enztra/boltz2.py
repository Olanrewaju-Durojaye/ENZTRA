"""Guarded Boltz-2 structural validation for strict kinetic survivors."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import shlex
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

from .boltz2_reporting import write_best_models_csv, write_boltz2_confidence_plot
from .config import EnztraConfig

Runner = Callable[..., subprocess.CompletedProcess]
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
CONFIDENCE_FIELDS = (
    "confidence_score",
    "ptm",
    "iptm",
    "ligand_iptm",
    "protein_iptm",
    "complex_plddt",
    "complex_iplddt",
    "complex_pde",
    "complex_ipde",
)


def _safe_stem(value: str) -> str:
    stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_.-")
    if not stem:
        raise ValueError(f"design identifier cannot form a safe filename: {value!r}")
    return stem


def _load_survivors(job_dir: Path) -> tuple[list[dict[str, str]], dict]:
    summary_path = job_dir / "summary.json"
    selection_path = job_dir / "selection.csv"
    manifest_path = job_dir / "manifest.json"
    for path in (summary_path, selection_path, manifest_path):
        if not path.is_file():
            raise FileNotFoundError(f"required completed-kinetics file is missing: {path}")

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    queue = summary.get("boltz2_queue")
    if not isinstance(queue, list):
        raise ValueError("summary.json contains no valid Boltz-2 queue")
    if not queue:
        raise ValueError("no design passed both strict kinetic thresholds")
    if len(queue) != len(set(queue)):
        raise ValueError("Boltz-2 queue contains duplicate design identifiers")
    smiles = str(manifest.get("substrate_smiles", "")).strip()
    if not smiles:
        raise ValueError("kinetic manifest contains no substrate SMILES")

    with selection_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    by_id: dict[str, dict[str, str]] = {}
    for row in rows:
        design_id = str(row.get("design_id", "")).strip()
        if not design_id:
            raise ValueError("selection.csv contains an empty design identifier")
        if design_id in by_id:
            raise ValueError(f"selection.csv contains duplicate design ID: {design_id}")
        by_id[design_id] = row

    survivors: list[dict[str, str]] = []
    for design_id in queue:
        if design_id not in by_id:
            raise ValueError(f"queued survivor is absent from selection.csv: {design_id}")
        row = by_id[design_id]
        if str(row.get("passes_strict_gate", "")).lower() != "true":
            raise ValueError(f"queued design did not pass the strict gate: {design_id}")
        sequence = str(row.get("sequence", "")).replace(" ", "").upper()
        invalid = sorted(set(sequence) - AMINO_ACIDS)
        if not sequence or invalid:
            detail = f": {', '.join(invalid)}" if invalid else ""
            raise ValueError(f"invalid survivor sequence for {design_id}{detail}")
        survivors.append({"design_id": design_id, "sequence": sequence})
    return survivors, manifest


def _write_boltz_yaml(path: Path, sequence: str, smiles: str, msa_mode: str) -> None:
    protein_lines = [
        "  - protein:",
        "      id: A",
        f"      sequence: {json.dumps(sequence)}",
    ]
    if msa_mode == "single":
        protein_lines.append("      msa: empty")
    lines = [
        "version: 1",
        "sequences:",
        *protein_lines,
        "  - ligand:",
        "      id: B",
        f"      smiles: {json.dumps(smiles)}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def find_boltz2_ready_jobs(jobs_root: Path) -> list[dict[str, object]]:
    """Return completed kinetic jobs with at least one strict survivor."""

    jobs: list[dict[str, object]] = []
    if not jobs_root.is_dir():
        return jobs
    for summary_path in jobs_root.rglob("summary.json"):
        job_dir = summary_path.parent
        if not (job_dir / "selection.csv").is_file() or not (
            job_dir / "manifest.json"
        ).is_file():
            continue
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            queue = summary.get("boltz2_queue", [])
        except (OSError, json.JSONDecodeError):
            continue
        if isinstance(queue, list) and queue:
            jobs.append(
                {
                    "job_dir": job_dir,
                    "name": str(job_dir.relative_to(jobs_root)),
                    "survivor_count": len(queue),
                    "modified": summary_path.stat().st_mtime,
                }
            )
    return sorted(jobs, key=lambda item: float(item["modified"]), reverse=True)


def find_completed_boltz2_jobs(jobs_root: Path) -> list[dict[str, object]]:
    """Return jobs whose latest Boltz-2 execution completed successfully."""

    jobs: list[dict[str, object]] = []
    if not jobs_root.is_dir():
        return jobs
    for status_path in jobs_root.rglob("boltz2/status.json"):
        job_dir = status_path.parent.parent
        execution_path = job_dir / "boltz2/execution.json"
        if not execution_path.is_file():
            continue
        try:
            status = json.loads(status_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if status.get("returncode") == 0 or status.get("stage") == "boltz2_complete":
            jobs.append({
                "job_dir": job_dir,
                "name": str(job_dir.relative_to(jobs_root)),
                "survivor_count": int(status.get("survivor_count", 0)),
                "modified": status_path.stat().st_mtime,
            })
    return sorted(jobs, key=lambda item: float(item["modified"]), reverse=True)


def refresh_boltz2_reports(job_dir: Path) -> tuple[Path, Path]:
    """Rebuild tables and plots from an already completed Boltz-2 run."""

    job_dir = job_dir.expanduser().resolve()
    execution_path = job_dir / "boltz2/execution.json"
    if not execution_path.is_file():
        raise FileNotFoundError(f"Boltz-2 execution plan is missing: {execution_path}")
    plan = json.loads(execution_path.read_text(encoding="utf-8"))
    csv_path, summary_path = collect_boltz2_results(job_dir, plan)
    report = json.loads(summary_path.read_text(encoding="utf-8"))
    status_path = job_dir / "boltz2/status.json"
    status = (
        json.loads(status_path.read_text(encoding="utf-8"))
        if status_path.is_file()
        else {}
    )
    status.update({
        "stage": "boltz2_complete",
        "returncode": 0,
        "confidence_csv": str(csv_path),
        "summary": str(summary_path),
        "best_models_csv": report["best_models_csv"],
        "confidence_plot": report["confidence_plot"],
        "reports_refreshed_at": datetime.now(timezone.utc).isoformat(),
    })
    status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    return Path(report["best_models_csv"]), Path(report["confidence_plot"])


def prepare_boltz2_execution(
    job_dir: Path,
    config: EnztraConfig,
    msa_mode: str = "single",
    diffusion_samples: int = 1,
    seed: int = 42,
) -> dict:
    """Create validated survivor YAML files and a reviewable Boltz-2 command."""

    job_dir = job_dir.expanduser().resolve()
    if msa_mode not in {"single", "server"}:
        raise ValueError("MSA mode must be 'single' or 'server'")
    if not 1 <= int(diffusion_samples) <= 25:
        raise ValueError("diffusion samples must be between 1 and 25")
    survivors, manifest = _load_survivors(job_dir)
    fingerprint_payload = {
        "survivors": survivors,
        "substrate_smiles": manifest["substrate_smiles"],
        "msa_mode": msa_mode,
        "diffusion_samples": int(diffusion_samples),
        "seed": int(seed),
    }
    run_id = hashlib.sha256(
        json.dumps(fingerprint_payload, sort_keys=True).encode("utf-8")
    ).hexdigest()[:12]
    # A fingerprinted directory prevents a changed survivor queue or sampling
    # plan from accidentally consuming inputs/results left by an earlier run.
    inputs_dir = job_dir / "boltz2/inputs" / run_id
    results_dir = job_dir / "boltz2/results" / run_id
    inputs_dir.mkdir(parents=True, exist_ok=True)

    used_stems: dict[str, str] = {}
    inputs: list[dict[str, str]] = []
    for survivor in survivors:
        stem = _safe_stem(survivor["design_id"])
        if stem in used_stems and used_stems[stem] != survivor["design_id"]:
            raise ValueError(
                "survivor identifiers collide after filename normalization: "
                f"{used_stems[stem]!r} and {survivor['design_id']!r}"
            )
        used_stems[stem] = survivor["design_id"]
        yaml_path = inputs_dir / f"{stem}.yaml"
        _write_boltz_yaml(
            yaml_path,
            survivor["sequence"],
            str(manifest["substrate_smiles"]),
            msa_mode,
        )
        inputs.append(
            {
                "design_id": survivor["design_id"],
                "input_stem": stem,
                "yaml": str(yaml_path),
            }
        )

    command = [
        config.conda_executable,
        "run",
        "--no-capture-output",
        "-n",
        config.boltz2_environment,
        "boltz",
        "predict",
        str(inputs_dir),
        "--model",
        "boltz2",
        "--accelerator",
        "gpu",
        "--devices",
        "1",
        "--diffusion_samples",
        str(diffusion_samples),
        "--max_parallel_samples",
        "1",
        "--seed",
        str(seed),
        "--output_format",
        "mmcif",
        "--cache",
        str(config.boltz_cache),
        "--out_dir",
        str(results_dir),
        "--override",
    ]
    if msa_mode == "server":
        command.append("--use_msa_server")
    plan = {
        "job_dir": str(job_dir),
        "run_id": run_id,
        "survivor_count": len(inputs),
        "survivor_ids": [item["design_id"] for item in inputs],
        "inputs": inputs,
        "msa_mode": msa_mode,
        "diffusion_samples": int(diffusion_samples),
        "seed": int(seed),
        "results_dir": str(results_dir),
        "command": command,
        "command_display": shlex.join(command),
        "affinity_requested": False,
    }
    (job_dir / "boltz2/execution.json").write_text(
        json.dumps(plan, indent=2) + "\n", encoding="utf-8"
    )
    return plan


def collect_boltz2_results(job_dir: Path, plan: dict) -> tuple[Path, Path]:
    """Validate all expected confidence/structure files and consolidate metrics."""

    job_dir = job_dir.expanduser().resolve()
    results_dir = Path(plan["results_dir"])
    rows: list[dict[str, object]] = []
    best: list[dict[str, object]] = []
    samples = int(plan["diffusion_samples"])
    for item in plan["inputs"]:
        design_id = item["design_id"]
        stem = item["input_stem"]
        design_rows: list[dict[str, object]] = []
        for model_index in range(samples):
            confidence_name = f"confidence_{stem}_model_{model_index}.json"
            matches = list(results_dir.rglob(confidence_name))
            if len(matches) != 1:
                raise ValueError(
                    f"expected one {confidence_name}, found {len(matches)}"
                )
            structure_name = f"{stem}_model_{model_index}.cif"
            structures = list(results_dir.rglob(structure_name))
            if len(structures) != 1:
                raise ValueError(
                    f"expected one {structure_name}, found {len(structures)}"
                )
            data = json.loads(matches[0].read_text(encoding="utf-8"))
            missing = [field for field in CONFIDENCE_FIELDS if field not in data]
            if missing:
                raise ValueError(
                    f"{confidence_name} lacks confidence fields: {', '.join(missing)}"
                )
            row: dict[str, object] = {
                "design_id": design_id,
                "model_index": model_index,
                "structure_file": str(structures[0]),
                "confidence_file": str(matches[0]),
            }
            for field in CONFIDENCE_FIELDS:
                value = float(data[field])
                if not math.isfinite(value):
                    raise ValueError(f"{confidence_name} contains non-finite {field}")
                if field in {
                    "confidence_score", "ptm", "iptm", "ligand_iptm",
                    "protein_iptm", "complex_plddt", "complex_iplddt",
                } and not 0.0 <= value <= 1.0:
                    raise ValueError(f"{confidence_name} contains invalid {field}")
                if field in {"complex_pde", "complex_ipde"} and value < 0.0:
                    raise ValueError(f"{confidence_name} contains invalid {field}")
                row[field] = value
            rows.append(row)
            design_rows.append(row)
        ordered_models = sorted(
            design_rows,
            key=lambda row: (
                float(row["confidence_score"]),
                float(row["ligand_iptm"]),
                float(row["complex_plddt"]),
            ),
            reverse=True,
        )
        for sample_rank, row in enumerate(ordered_models, 1):
            row["sample_rank"] = sample_rank
            row["is_best_model"] = sample_rank == 1
        best.append(ordered_models[0])

    best.sort(
        key=lambda row: (
            float(row["confidence_score"]),
            float(row["ligand_iptm"]),
            float(row["complex_plddt"]),
        ),
        reverse=True,
    )
    selection: dict[str, dict[str, str]] = {}
    with (job_dir / "selection.csv").open(newline="", encoding="utf-8") as handle:
        selection = {row["design_id"]: row for row in csv.DictReader(handle)}
    for structural_rank, row in enumerate(best, 1):
        row["structural_rank"] = structural_rank
        kinetic = selection.get(str(row["design_id"]), {})
        row["kinetic_survivor_rank"] = kinetic.get("survivor_rank", "")
        row["kcat_s"] = kinetic.get("kcat_s", "")
        row["km_mm"] = kinetic.get("km_mm", "")
        row["catalytic_efficiency_s-1_mM-1"] = kinetic.get(
            "catalytic_efficiency_s-1_mM-1", ""
        )

    csv_path = job_dir / "boltz2/confidence.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "design_id",
                "model_index",
                "sample_rank",
                "is_best_model",
                *CONFIDENCE_FIELDS,
                "structure_file",
                "confidence_file",
            ),
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)

    reports_dir = job_dir / "boltz2/results/reports"
    best_csv = write_best_models_csv(reports_dir / "best_models.csv", best)
    plot_path = write_boltz2_confidence_plot(
        reports_dir / "boltz2_confidence.svg", best
    )
    summary_path = job_dir / "boltz2/summary.json"
    summary_path.write_text(
        json.dumps(
            {
                "survivor_count": len(plan["inputs"]),
                "diffusion_samples": samples,
                "msa_mode": plan["msa_mode"],
                "affinity_requested": False,
                "ranking_method": (
                    "Descending confidence_score, then ligand_iptm, then "
                    "complex_plddt; ties retain deterministic input order."
                ),
                "all_models_csv": str(csv_path),
                "best_models_csv": str(best_csv),
                "confidence_plot": str(plot_path),
                "best_models": best,
                "interpretation_note": (
                    "Boltz-2 confidence estimates structural reliability; it does "
                    "not establish catalytic activity or experimental binding affinity."
                ),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return csv_path, summary_path


def run_boltz2(
    job_dir: Path,
    config: EnztraConfig,
    msa_mode: str = "single",
    diffusion_samples: int = 1,
    seed: int = 42,
    runner: Runner = subprocess.run,
) -> int:
    """Run Boltz-2 and persist guarded completion or failure state."""

    job_dir = job_dir.expanduser().resolve()
    plan = prepare_boltz2_execution(
        job_dir, config, msa_mode, diffusion_samples, seed
    )
    status_path = job_dir / "boltz2/status.json"
    status = {
        "stage": "boltz2_running",
        "started_at": datetime.now(timezone.utc).isoformat(),
        "survivor_count": plan["survivor_count"],
    }
    status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    try:
        completed = runner(plan["command"], cwd=job_dir, check=False)
        returncode = int(completed.returncode)
        confidence_csv = summary_path = None
        if returncode == 0:
            try:
                confidence_csv, summary_path = collect_boltz2_results(job_dir, plan)
                report = json.loads(summary_path.read_text(encoding="utf-8"))
                best_models_csv = report["best_models_csv"]
                confidence_plot = report["confidence_plot"]
            except (FileNotFoundError, ValueError, json.JSONDecodeError) as error:
                status["validation_error"] = str(error)
                returncode = 2
    except BaseException:
        status.update(
            {
                "stage": "boltz2_interrupted",
                "finished_at": datetime.now(timezone.utc).isoformat(),
            }
        )
        status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
        raise
    status.update(
        {
            "stage": "boltz2_complete" if returncode == 0 else "boltz2_failed",
            "returncode": returncode,
            "confidence_csv": str(confidence_csv) if confidence_csv else "",
            "summary": str(summary_path) if summary_path else "",
            "best_models_csv": best_models_csv if returncode == 0 else "",
            "confidence_plot": confidence_plot if returncode == 0 else "",
            "finished_at": datetime.now(timezone.utc).isoformat(),
        }
    )
    status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    if returncode == 0:
        print(f"Boltz-2 best-model table: {best_models_csv}")
        print(f"Boltz-2 confidence plot: {confidence_plot}")
    return returncode
