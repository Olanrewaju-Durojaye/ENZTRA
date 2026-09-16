"""Consolidated ranking, quality control, and publication exports."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from . import __version__
from .boltz2_reporting import write_boltz2_confidence_plot
from .models import Candidate, ReferenceThresholds
from .plotting import write_kinetic_selection_plot
from .selection import evaluate_candidates

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
STRUCTURAL_THRESHOLDS = {
    "confidence_score": 0.70,
    "ligand_iptm": 0.70,
    "complex_plddt": 0.70,
}


def _read_json(path: Path, required: bool = True) -> dict[str, Any]:
    if not path.is_file():
        if required:
            raise FileNotFoundError(f"required result file is missing: {path}")
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _read_csv(path: Path, required: bool = True) -> list[dict[str, str]]:
    if not path.is_file():
        if required:
            raise FileNotFoundError(f"required result file is missing: {path}")
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _number(row: dict[str, str], field: str) -> float | None:
    value = str(row.get(field, "")).strip()
    if not value:
        return None
    try:
        number = float(value)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def _truth(value: object) -> bool:
    return str(value).strip().lower() == "true"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _source_metadata(job_dir: Path) -> dict[str, dict[str, str]]:
    rows = _read_csv(job_dir / "imported_inputs/source_metadata.csv", required=False)
    return {row.get("design_id", ""): row for row in rows if row.get("design_id")}


def _best_models(job_dir: Path) -> dict[str, dict[str, str]]:
    candidates = [
        job_dir / "boltz2/results/reports/best_models.csv",
        job_dir / "boltz2/best_models.csv",
    ]
    status = _read_json(job_dir / "boltz2/status.json", required=False)
    reported = str(status.get("best_models_csv", "")).strip()
    if reported:
        path = Path(reported)
        candidates.insert(0, path if path.is_absolute() else job_dir / path)
    for path in candidates:
        rows = _read_csv(path, required=False)
        if rows:
            return {row["design_id"]: row for row in rows}
    return {}


def _qc_flags(
    row: dict[str, str], structure: dict[str, str] | None, boltz_expected: bool,
    job_dir: Path, km_uncertainty_review_threshold: float | None,
) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []
    sequence = str(row.get("sequence", "")).replace(" ", "").upper()
    invalid = sorted(set(sequence) - AMINO_ACIDS)
    if not sequence:
        errors.append("missing_sequence")
    elif invalid:
        errors.append("invalid_amino_acid_symbols")
    for field in ("kcat_s", "km_mm", "catalytic_efficiency_s-1_mM-1"):
        value = _number(row, field)
        if value is None or value <= 0:
            errors.append(f"invalid_{field}")
    km_uncertainty = _number(row, "km_sd_total")
    if km_uncertainty is None:
        warnings.append("missing_km_uncertainty")
    elif (
        km_uncertainty_review_threshold is not None
        and km_uncertainty > km_uncertainty_review_threshold
    ):
        warnings.append("km_uncertainty_above_review_threshold")
    strict = _truth(row.get("passes_strict_gate"))
    if strict and not row.get("survivor_rank"):
        errors.append("missing_kinetic_survivor_rank")
    if strict and boltz_expected and structure is None:
        errors.append("missing_boltz2_best_model")
    elif strict and not boltz_expected and structure is None:
        warnings.append("structural_validation_not_completed")
    if structure:
        for field, threshold in STRUCTURAL_THRESHOLDS.items():
            value = _number(structure, field)
            if value is None:
                errors.append(f"missing_{field}")
            elif value < threshold:
                warnings.append(f"low_{field}")
        structure_file = str(structure.get("structure_file", "")).strip()
        structure_path = Path(structure_file) if structure_file else Path()
        if structure_file and not structure_path.is_absolute():
            structure_path = job_dir / structure_path
        if not structure_file or not structure_path.is_file():
            errors.append("missing_structure_file")
    return errors, warnings


def _priority_score(row: dict[str, str], structure: dict[str, str]) -> float | None:
    kcat_fold = _number(row, "kcat_fold_change")
    km_fold = _number(row, "km_fold_change")
    confidence = _number(structure, "confidence_score")
    ligand_iptm = _number(structure, "ligand_iptm")
    complex_plddt = _number(structure, "complex_plddt")
    values = (kcat_fold, km_fold, confidence, ligand_iptm, complex_plddt)
    if any(value is None or value <= 0 for value in values):
        return None
    # Dimensionless geometric mean; lower Km becomes a benefit through 1/Km-fold.
    product = kcat_fold * (1.0 / km_fold) * confidence * ligand_iptm * complex_plddt
    return product ** (1.0 / 5.0)


MASTER_FIELDS = (
    "final_priority_rank", "publication_label", "combined_priority_score", "design_id",
    "original_header", "sequence", "decision", "passes_strict_gate",
    "kinetic_survivor_rank", "kcat_s", "km_mm",
    "km_sd_total", "km_sd_aleatoric", "km_sd_epistemic",
    "catalytic_efficiency_s-1_mM-1", "kcat_fold_change", "km_fold_change",
    "efficiency_fold_change", "structural_rank", "model_index",
    "confidence_score", "ligand_iptm", "complex_plddt", "structure_file",
    "qc_status", "qc_errors", "qc_warnings",
)


def _write_csv(path: Path, fieldnames: tuple[str, ...], rows: list[dict[str, Any]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return path


def _write_outcome_svg(path: Path, accepted: int, rejected: int, qc_failed: int) -> Path:
    maximum = max(accepted, rejected, qc_failed, 1)
    data = (("Strict survivors", accepted, "#13795b"),
            ("Rejected", rejected, "#c65b4b"),
            ("QC failures", qc_failed, "#c38b2a"))
    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="430" viewBox="0 0 900 430" role="img">',
        '<title>ENZTRA candidate outcomes and quality control</title>',
        '<rect width="900" height="430" fill="#f4f2ea"/>',
        '<rect x="30" y="25" width="840" height="375" rx="18" fill="#fffdf7" stroke="#d9ddd5"/>',
        '<text x="65" y="70" font-family="Georgia,serif" font-size="27" font-weight="700" fill="#10231d">Candidate outcomes and quality control</text>',
        '<text x="65" y="96" font-family="sans-serif" font-size="13" fill="#64736d">Strict kinetic decisions remain separate from data-quality checks.</text>',
    ]
    for index, (label, value, colour) in enumerate(data):
        y = 145 + index * 78
        width = 600 * value / maximum
        lines.extend([
            f'<text x="65" y="{y + 23}" font-family="sans-serif" font-size="15" fill="#10231d">{label}</text>',
            f'<rect x="220" y="{y}" width="600" height="32" rx="7" fill="#e8ece7"/>',
            f'<rect x="220" y="{y}" width="{width:.2f}" height="32" rx="7" fill="{colour}"/>',
            f'<text x="{min(830, 230 + width):.2f}" y="{y + 22}" font-family="sans-serif" font-size="14" font-weight="700" fill="#10231d">{value}</text>',
        ])
    lines.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _copy_figure(source: Path, destination: Path) -> str:
    if not source.is_file():
        return ""
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return str(destination)


def _artifact_inventory(job_dir: Path, excluded_root: Path) -> list[dict[str, Any]]:
    inventory = []
    for path in sorted(job_dir.rglob("*")):
        if not path.is_file() or excluded_root in path.parents:
            continue
        inventory.append({
            "path": str(path.relative_to(job_dir)),
            "size_bytes": path.stat().st_size,
            "sha256": _sha256(path),
        })
    return inventory


def generate_publication_exports(
    job_dir: Path, km_uncertainty_review_threshold: float | None = None
) -> dict[str, Any]:
    """Build a self-contained, auditable publication export directory."""

    job_dir = job_dir.expanduser().resolve()
    if (
        km_uncertainty_review_threshold is not None
        and (not math.isfinite(km_uncertainty_review_threshold)
             or km_uncertainty_review_threshold <= 0)
    ):
        raise ValueError("Km uncertainty review threshold must be positive and finite")
    selection = _read_csv(job_dir / "selection.csv")
    summary = _read_json(job_dir / "summary.json")
    manifest = _read_json(job_dir / "manifest.json")
    if not selection:
        raise ValueError("selection.csv contains no candidate rows")
    source = _source_metadata(job_dir)
    structures = _best_models(job_dir)
    kinetics = {
        row.get("design_id", ""): row
        for row in _read_csv(job_dir / "kinetics.csv", required=False)
        if row.get("design_id")
    }
    boltz_status = _read_json(job_dir / "boltz2/status.json", required=False)
    boltz_expected = bool(summary.get("boltz2_queue")) and (
        boltz_status.get("stage") == "boltz2_complete"
        or boltz_status.get("returncode") == 0
    )
    rows: list[dict[str, Any]] = []
    qc_rows: list[dict[str, Any]] = []
    for kinetic in selection:
        design_id = kinetic.get("design_id", "")
        merged_kinetic = {**kinetic}
        merged_kinetic.update({
            field: kinetics.get(design_id, {}).get(field, "")
            for field in ("km_sd_total", "km_sd_aleatoric", "km_sd_epistemic")
        })
        structure = structures.get(design_id)
        errors, warnings = _qc_flags(
            merged_kinetic, structure, boltz_expected, job_dir,
            km_uncertainty_review_threshold,
        )
        score = _priority_score(merged_kinetic, structure) if structure and _truth(
            merged_kinetic.get("passes_strict_gate")
        ) else None
        metadata = source.get(design_id, {})
        row: dict[str, Any] = {
            **merged_kinetic,
            "kinetic_survivor_rank": merged_kinetic.get("survivor_rank", ""),
            "original_header": metadata.get("original_header", ""),
            "combined_priority_score": f"{score:.8f}" if score is not None else "",
            "qc_status": "FAIL" if errors else ("WARN" if warnings else "PASS"),
            "qc_errors": ";".join(errors),
            "qc_warnings": ";".join(warnings),
        }
        if structure:
            row.update({field: structure.get(field, "") for field in (
                "structural_rank", "model_index", "confidence_score",
                "ligand_iptm", "complex_plddt", "structure_file",
            )})
        rows.append(row)
        qc_rows.append({
            "design_id": design_id, "qc_status": row["qc_status"],
            "errors": row["qc_errors"], "warnings": row["qc_warnings"],
            "strict_kinetic_survivor": merged_kinetic.get("passes_strict_gate", "False"),
            "boltz2_expected": str(boltz_expected and _truth(merged_kinetic.get("passes_strict_gate"))),
        })

    priority = sorted(
        (row for row in rows if row["combined_priority_score"]),
        key=lambda row: float(row["combined_priority_score"]), reverse=True,
    )
    for rank, row in enumerate(priority, 1):
        row["final_priority_rank"] = rank
        row["publication_label"] = f"Priority {rank}"
    for row in rows:
        row.setdefault("final_priority_rank", "")
        row.setdefault("publication_label", "")

    output = job_dir / "publication"
    tables = output / "tables"
    figures = output / "figures"
    master_path = _write_csv(tables / "master_ranking.csv", MASTER_FIELDS, rows)
    master_json = tables / "master_ranking.json"
    master_json.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    qc_path = _write_csv(
        tables / "quality_control.csv",
        ("design_id", "qc_status", "errors", "warnings",
         "strict_kinetic_survivor", "boltz2_expected"), qc_rows,
    )
    accepted = sum(_truth(row.get("passes_strict_gate")) for row in selection)
    qc_failed = sum(row["qc_status"] == "FAIL" for row in rows)
    outcome_plot = _write_outcome_svg(
        figures / "candidate_outcomes.svg", accepted, len(selection) - accepted, qc_failed
    )
    labels = {
        str(row["design_id"]): str(row["publication_label"])
        for row in priority
    }
    reference_data = summary.get("reference", {})
    reference = ReferenceThresholds(
        str(reference_data["enzyme_id"]),
        str(reference_data.get("substrate_id", manifest.get("substrate_name", "substrate"))),
        float(reference_data["kcat_s"]),
        float(reference_data["km_mm"]),
    )
    candidates = [
        Candidate(
            str(row["design_id"]), str(row["sequence"]),
            float(row["kcat_s"]), float(row["km_mm"]),
        )
        for row in selection
    ]
    kinetic_plot = write_kinetic_selection_plot(
        figures / "kinetic_selection.svg",
        reference,
        evaluate_candidates(candidates, reference),
        display_labels=labels,
    )
    copied_figures = {
        "kinetic_selection": str(kinetic_plot),
        "boltz2_confidence": "",
    }
    if priority:
        boltz_plot = write_boltz2_confidence_plot(
            figures / "boltz2_confidence.svg",
            priority,
            display_labels=labels,
            rank_field="final_priority_rank",
        )
        copied_figures["boltz2_confidence"] = str(boltz_plot)

    qc_counts = {name: sum(row["qc_status"] == name for row in rows) for name in ("PASS", "WARN", "FAIL")}
    publication_summary = {
        "enztra_version": __version__,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "job": str(job_dir),
        "reference": summary.get("reference", {}),
        "candidate_count": len(rows),
        "strict_survivor_count": accepted,
        "structurally_ranked_count": len(priority),
        "qc_counts": qc_counts,
        "km_uncertainty_review_threshold": km_uncertainty_review_threshold,
        "priority_score_definition": "geometric_mean(kcat_fold, 1/km_fold, confidence_score, ligand_iptm, complex_plddt)",
        "important_interpretation": "Combined priority does not change the strict kinetic gate and is not experimental evidence.",
        "outputs": {
            "master_ranking": str(master_path), "master_ranking_json": str(master_json),
            "quality_control": str(qc_path),
            "candidate_outcomes": str(outcome_plot), **copied_figures,
        },
    }
    summary_path = output / "publication_summary.json"
    _write = lambda p, d: p.write_text(json.dumps(d, indent=2) + "\n", encoding="utf-8")
    _write(summary_path, publication_summary)

    commands = {}
    for name in ("rfdiffusion2_execution.json", "ligandmpnn_execution.json"):
        data = _read_json(job_dir / name, required=False)
        if data:
            commands[name.removesuffix("_execution.json")] = data.get("command_display", data.get("command", []))
    boltz_execution = _read_json(job_dir / "boltz2/execution.json", required=False)
    if boltz_execution:
        commands["boltz2"] = boltz_execution.get("command_display", boltz_execution.get("command", []))
    reproducibility = {
        "schema_version": 1,
        "enztra_version": __version__,
        "generated_at": publication_summary["generated_at"],
        "job_id": manifest.get("job_id", job_dir.name),
        "design_request": _read_json(job_dir / "design_request.json", required=False),
        "kinetic_manifest": manifest,
        "workflow_state": _read_json(job_dir / "workflow_state.json", required=False),
        "commands": commands,
        "input_and_result_inventory": _artifact_inventory(job_dir, output),
        "checksum_note": "SHA-256 values identify file changes; users do not need to calculate them manually.",
    }
    reproducibility_path = output / "reproducibility_manifest.json"
    _write(reproducibility_path, reproducibility)

    report = output / "publication_report.md"
    report.write_text(
        "\n".join([
            f"# ENZTRA publication export: {job_dir.name}", "",
            f"- Candidates evaluated: {len(rows)}",
            f"- Strict kinetic survivors: {accepted}",
            f"- Structurally ranked survivors: {len(priority)}",
            f"- QC: {qc_counts['PASS']} pass, {qc_counts['WARN']} warning, {qc_counts['FAIL']} fail", "",
            "The strict gate requires kcat above and Km below the reference. The combined priority score is used only after strict acceptance and available Boltz-2 validation; it does not rescue rejected candidates.", "",
            "All reported values are computational predictions and require experimental validation.",
        ]) + "\n", encoding="utf-8"
    )
    publication_summary["outputs"].update({
        "publication_summary": str(summary_path),
        "reproducibility_manifest": str(reproducibility_path),
        "publication_report": str(report),
    })
    _write(summary_path, publication_summary)
    return publication_summary
