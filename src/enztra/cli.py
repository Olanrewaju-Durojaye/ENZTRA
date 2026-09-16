"""Command-line entry point for the first ENZTRA milestone."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any

from .adapters import (
    combine_predictions,
    prepare_inputs,
    run_catpred,
    run_dlkcat,
)
from .config import EnztraConfig
from .doctor import run_doctor
from .fasta import read_fasta
from .models import Candidate, ReferenceThresholds
from .plotting import write_kinetic_selection_plot
from .selection import evaluate_candidates


def _optional_float(row: dict[str, str], key: str) -> float | None:
    value = row.get(key, "").strip()
    return float(value) if value else None


def _load_reference(path: Path) -> ReferenceThresholds:
    data: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    return ReferenceThresholds(
        enzyme_id=str(data["enzyme_id"]),
        substrate_id=str(data["substrate_id"]),
        kcat_s=float(data["kcat_s"]),
        km_mm=float(data["km_mm"]),
    )


def _load_candidates(path: Path) -> list[Candidate]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = csv.DictReader(handle)
        return [
            Candidate(
                design_id=row["design_id"],
                sequence=row["sequence"],
                kcat_s=float(row["kcat_s"]),
                km_mm=float(row["km_mm"]),
                kcat_uncertainty=_optional_float(row, "kcat_uncertainty"),
                km_uncertainty=_optional_float(row, "km_uncertainty"),
            )
            for row in rows
        ]


def _write_results(path: Path, results: list) -> None:
    fieldnames = [
        "design_id",
        "sequence",
        "kcat_s",
        "km_mm",
        "catalytic_efficiency_s-1_mM-1",
        "passes_kcat",
        "passes_km",
        "passes_strict_gate",
        "decision",
        "kcat_fold_change",
        "km_fold_change",
        "efficiency_fold_change",
        "survivor_rank",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            candidate = result.candidate
            writer.writerow(
                {
                    "design_id": candidate.design_id,
                    "sequence": candidate.sequence,
                    "kcat_s": candidate.kcat_s,
                    "km_mm": candidate.km_mm,
                    "catalytic_efficiency_s-1_mM-1": candidate.catalytic_efficiency,
                    "passes_kcat": result.passes_kcat,
                    "passes_km": result.passes_km,
                    "passes_strict_gate": result.passes_strict_gate,
                    "decision": result.decision,
                    "kcat_fold_change": result.kcat_fold_change,
                    "km_fold_change": result.km_fold_change,
                    "efficiency_fold_change": result.efficiency_fold_change,
                    "survivor_rank": result.rank or "",
                }
            )


def _select(args: argparse.Namespace) -> int:
    reference = _load_reference(args.reference)
    candidates = _load_candidates(args.designs)
    results = evaluate_candidates(candidates, reference)
    _write_results(args.output, results)
    survivors = sorted(
        (result for result in results if result.passes_strict_gate),
        key=lambda result: result.rank or 0,
    )
    summary = {
        "reference": {
            "enzyme_id": reference.enzyme_id,
            "substrate_id": reference.substrate_id,
            "kcat_s": reference.kcat_s,
            "km_mm": reference.km_mm,
            "catalytic_efficiency_s-1_mM-1": reference.catalytic_efficiency,
        },
        "design_count": len(results),
        "survivor_count": len(survivors),
        "survivor_ids": [result.candidate.design_id for result in survivors],
        "boltz2_queue": [result.candidate.design_id for result in survivors],
        "kinetic_plot": "results/kinetic_selection.svg",
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    plot_path = write_kinetic_selection_plot(
        args.output.parent / "results" / "kinetic_selection.svg", reference, results
    )
    print(
        f"Evaluated {len(results)} designs: {len(survivors)} strict survivors. "
        f"Results: {args.output}\nKinetic plot: {plot_path}"
    )
    return 0


def _doctor(args: argparse.Namespace) -> int:
    report = run_doctor(EnztraConfig.from_json(args.config))
    for check in report["checks"]:
        marker = "PASS" if check["ok"] else "FAIL"
        print(f"[{marker}] {check['name']}: {check['detail']}")
    print(f"\nENZTRA status: {report['status'].upper()}")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return 0 if report["status"] == "ready" else 1


def _write_kinetics(path: Path, predictions: list, reference_id: str) -> None:
    fieldnames = [
        "design_id",
        "role",
        "sequence",
        "kcat_s",
        "km_mm",
        "km_sd_total",
        "km_sd_aleatoric",
        "km_sd_epistemic",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for prediction in predictions:
            writer.writerow(
                {
                    "design_id": prediction.design_id,
                    "role": (
                        "reference" if prediction.design_id == reference_id else "design"
                    ),
                    "sequence": prediction.sequence,
                    "kcat_s": prediction.kcat_s,
                    "km_mm": prediction.km_mm,
                    "km_sd_total": prediction.km_sd_total or "",
                    "km_sd_aleatoric": prediction.km_sd_aleatoric or "",
                    "km_sd_epistemic": prediction.km_sd_epistemic or "",
                }
            )


def _kinetics(args: argparse.Namespace) -> int:
    config = EnztraConfig.from_json(args.config)
    reference_records = read_fasta(args.reference_fasta)
    if len(reference_records) != 1:
        raise ValueError("The reference FASTA must contain exactly one record")
    design_records = read_fasta(args.designs_fasta)
    records = reference_records + design_records
    identifiers = [record.record_id for record in records]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("Reference and design FASTA identifiers must be unique")

    args.job_dir.mkdir(parents=True, exist_ok=True)
    dlkcat_input, catpred_input = prepare_inputs(
        records,
        args.substrate_name,
        args.substrate_smiles,
        args.job_dir / "inputs",
    )
    manifest = {
        "reference_id": reference_records[0].record_id,
        "design_ids": [record.record_id for record in design_records],
        "substrate_name": args.substrate_name,
        "substrate_smiles": args.substrate_smiles,
        "inputs": {
            "dlkcat": str(dlkcat_input),
            "catpred": str(catpred_input),
        },
    }
    (args.job_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    if args.prepare_only:
        print(f"Prepared {len(records)} enzyme–substrate pairs in {args.job_dir}")
        return 0

    raw_dir = args.job_dir / "raw"
    dlkcat_output = raw_dir / "dlkcat_output.tsv"
    catpred_output = raw_dir / "catpred_output.csv"
    if args.reuse_raw:
        missing = [
            str(path) for path in (dlkcat_output, catpred_output) if not path.is_file()
        ]
        if missing:
            raise FileNotFoundError(
                "--reuse-raw requested, but required files are missing: "
                + ", ".join(missing)
            )
    else:
        dlkcat_output = run_dlkcat(config, dlkcat_input, dlkcat_output)
        catpred_output = run_catpred(config, catpred_input, catpred_output)
    predictions = combine_predictions(records, dlkcat_output, catpred_output)
    kinetics_path = args.job_dir / "kinetics.csv"
    _write_kinetics(kinetics_path, predictions, reference_records[0].record_id)

    reference_prediction = predictions[0]
    reference = ReferenceThresholds(
        enzyme_id=reference_prediction.design_id,
        substrate_id=args.substrate_name,
        kcat_s=reference_prediction.kcat_s,
        km_mm=reference_prediction.km_mm,
    )
    candidates = [
        Candidate(
            design_id=prediction.design_id,
            sequence=prediction.sequence,
            kcat_s=prediction.kcat_s,
            km_mm=prediction.km_mm,
            km_uncertainty=prediction.km_sd_total,
        )
        for prediction in predictions[1:]
    ]
    results = evaluate_candidates(candidates, reference)
    _write_results(args.job_dir / "selection.csv", results)
    survivors = sorted(
        (result for result in results if result.passes_strict_gate),
        key=lambda result: result.rank or 0,
    )
    summary = {
        "reference": {
            "enzyme_id": reference.enzyme_id,
            "substrate_id": reference.substrate_id,
            "kcat_s": reference.kcat_s,
            "km_mm": reference.km_mm,
            "catalytic_efficiency_s-1_mM-1": reference.catalytic_efficiency,
        },
        "design_count": len(candidates),
        "survivor_count": len(survivors),
        "survivor_ids": [result.candidate.design_id for result in survivors],
        "boltz2_queue": [result.candidate.design_id for result in survivors],
        "kinetic_plot": "results/kinetic_selection.svg",
    }
    (args.job_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    plot_path = write_kinetic_selection_plot(
        args.job_dir / "results" / "kinetic_selection.svg", reference, results
    )
    print(
        f"Predicted {len(records)} pairs and retained {len(survivors)} of "
        f"{len(candidates)} designs. Job: {args.job_dir}\n"
        f"Kinetic plot: {plot_path}"
    )
    return 0


def _serve(args: argparse.Namespace) -> int:
    """Launch the dependency-free local results dashboard."""

    from .web import serve

    return serve(args.jobs_root, args.host, args.port)


def _menu(args: argparse.Namespace) -> int:
    from .interactive import run_menu

    return run_menu(args.jobs_root, args.config)


def _prepare_design(args: argparse.Namespace) -> int:
    from .interactive import prepare_design_interactively

    return prepare_design_interactively(args.jobs_root)


def _boltz2(args: argparse.Namespace) -> int:
    """Prepare or run structural validation for strict kinetic survivors."""

    from .boltz2 import prepare_boltz2_execution, run_boltz2

    config = EnztraConfig.from_json(args.config)
    plan = prepare_boltz2_execution(
        args.job_dir,
        config,
        args.msa_mode,
        args.diffusion_samples,
        args.seed,
    )
    if args.prepare_only:
        print(
            f"Prepared {plan['survivor_count']} survivor inputs in "
            f"{args.job_dir / 'boltz2/inputs'}"
        )
        return 0
    return run_boltz2(
        args.job_dir,
        config,
        args.msa_mode,
        args.diffusion_samples,
        args.seed,
    )


def _boltz2_report(args: argparse.Namespace) -> int:
    """Regenerate reports from an existing completed Boltz-2 execution."""

    from .boltz2 import refresh_boltz2_reports

    best_csv, plot_path = refresh_boltz2_reports(args.job_dir)
    print(f"Boltz-2 best-model table: {best_csv}")
    print(f"Boltz-2 confidence plot: {plot_path}")
    return 0


def _workflow(args: argparse.Namespace) -> int:
    """Run or resume all computational stages for one prepared project."""

    from .workflow import run_workflow

    return run_workflow(
        args.job_dir,
        args.config,
        msa_mode=args.msa_mode,
        diffusion_samples=args.diffusion_samples,
        seed=args.seed,
    )


def _publication(args: argparse.Namespace) -> int:
    """Generate consolidated QC and publication-ready project exports."""

    from .publication import generate_publication_exports

    result = generate_publication_exports(
        args.job_dir, args.km_uncertainty_review_threshold
    )
    print(f"Publication exports: {args.job_dir.expanduser().resolve() / 'publication'}")
    print(f"Candidates: {result['candidate_count']}")
    print(f"Strict survivors: {result['strict_survivor_count']}")
    print(f"Structurally ranked: {result['structurally_ranked_count']}")
    return 0


def _pilot(args: argparse.Namespace) -> int:
    """Prepare or execute a controlled native ENZTRA pilot."""

    from .pilot import prepare_pilot, run_controlled_pilot

    if args.prepare_only:
        report = prepare_pilot(args.job_dir, args.config)
        print(f"Pilot preflight: {args.job_dir / 'pilot/preflight.json'}")
        print(f"Planned designs: {report['planned_sequences']}")
        print(f"Status: {'READY' if report['ready'] else 'BLOCKED'}")
        return 0 if report["ready"] else 2
    return run_controlled_pilot(
        args.job_dir, args.config, msa_mode=args.msa_mode,
        diffusion_samples=args.diffusion_samples, seed=args.seed,
    )


def _demonstration(args: argparse.Namespace) -> int:
    """Preflight or execute the exact v1 release demonstration."""

    from .demonstration import (
        prepare_release_demonstration,
        run_release_demonstration,
    )

    if args.prepare_only:
        report = prepare_release_demonstration(
            args.job_dir, args.config, args.calibration_job
        )
        print(f"Release preflight: {args.job_dir / 'release/preflight.json'}")
        plan = report["plan"]
        print(
            "Plan: "
            f"{plan['backbones']} backbones × "
            f"{plan['sequences_per_backbone']} sequences = "
            f"{plan['total_sequences']} designs"
        )
        print(f"Estimated output: {report['estimated_output_gib']:.2f} GiB")
        print(f"Available disk: {report['available_free_gib']:.2f} GiB")
        for blocker in report["blockers"]:
            print(f"BLOCKER: {blocker}")
        print(f"Status: {'READY' if report['ready'] else 'BLOCKED'}")
        return 0 if report["ready"] else 2
    return run_release_demonstration(
        args.job_dir,
        args.config,
        calibration_job=args.calibration_job,
        msa_mode=args.msa_mode,
        diffusion_samples=args.diffusion_samples,
        seed=args.seed,
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="enztra",
        description="Reference-guided enzyme redesign and triage",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    workflow = subparsers.add_parser(
        "workflow", help="Run or resume the complete ENZTRA workflow"
    )
    workflow.add_argument("--config", type=Path, required=True)
    workflow.add_argument("--job-dir", type=Path, required=True)
    workflow.add_argument(
        "--msa-mode", choices=("single", "server"), default="single"
    )
    workflow.add_argument(
        "--diffusion-samples", type=int, choices=range(1, 26), default=1
    )
    workflow.add_argument("--seed", type=int, default=42)
    workflow.set_defaults(handler=_workflow)

    publication = subparsers.add_parser(
        "publication-export",
        help="Generate consolidated rankings, QC, and publication files",
    )
    publication.add_argument("--job-dir", type=Path, required=True)
    publication.add_argument(
        "--km-uncertainty-review-threshold",
        type=float,
        help="Optional CatPred Km SD threshold for review warnings",
    )
    publication.set_defaults(handler=_publication)

    pilot = subparsers.add_parser(
        "pilot", help="Preflight and run/resume a controlled native pilot"
    )
    pilot.add_argument("--config", type=Path, required=True)
    pilot.add_argument("--job-dir", type=Path, required=True)
    pilot.add_argument("--msa-mode", choices=("single", "server"), default="single")
    pilot.add_argument("--diffusion-samples", type=int, choices=range(1, 26), default=1)
    pilot.add_argument("--seed", type=int, default=42)
    pilot.add_argument("--prepare-only", action="store_true")
    pilot.set_defaults(handler=_pilot)

    demonstration = subparsers.add_parser(
        "demonstrate",
        help="Preflight and run/resume the exact 500-sequence v1 qualification",
    )
    demonstration.add_argument("--config", type=Path, required=True)
    demonstration.add_argument("--job-dir", type=Path, required=True)
    demonstration.add_argument("--calibration-job", type=Path)
    demonstration.add_argument(
        "--msa-mode", choices=("single", "server"), default="single"
    )
    demonstration.add_argument(
        "--diffusion-samples", type=int, choices=range(1, 26), default=1
    )
    demonstration.add_argument("--seed", type=int, default=42)
    demonstration.add_argument("--prepare-only", action="store_true")
    demonstration.set_defaults(handler=_demonstration)

    select = subparsers.add_parser("select", help="Apply the strict kinetic gate")
    select.add_argument("--reference", type=Path, required=True)
    select.add_argument("--designs", type=Path, required=True)
    select.add_argument("--output", type=Path, required=True)
    select.add_argument("--summary", type=Path, required=True)
    select.set_defaults(handler=_select)

    doctor = subparsers.add_parser("doctor", help="Check external tool installations")
    doctor.add_argument("--config", type=Path, required=True)
    doctor.add_argument("--output", type=Path)
    doctor.set_defaults(handler=_doctor)

    kinetics = subparsers.add_parser(
        "kinetics",
        help="Batch DLKcat/CatPred prediction followed by strict selection",
    )
    kinetics.add_argument("--config", type=Path, required=True)
    kinetics.add_argument("--reference-fasta", type=Path, required=True)
    kinetics.add_argument("--designs-fasta", type=Path, required=True)
    kinetics.add_argument("--substrate-name", required=True)
    kinetics.add_argument("--substrate-smiles", required=True)
    kinetics.add_argument("--job-dir", type=Path, required=True)
    kinetics.add_argument(
        "--prepare-only",
        action="store_true",
        help="Validate FASTA and write model inputs without running predictions",
    )
    kinetics.add_argument(
        "--reuse-raw",
        action="store_true",
        help="Skip model execution and combine existing raw outputs in the job directory",
    )
    kinetics.set_defaults(handler=_kinetics)

    serve = subparsers.add_parser(
        "serve", help="Open the local browser dashboard for completed jobs"
    )
    serve.add_argument(
        "--jobs-root",
        type=Path,
        default=Path("jobs"),
        help="Directory containing ENZTRA job folders (default: jobs)",
    )
    serve.add_argument("--host", default="127.0.0.1")
    serve.add_argument("--port", type=int, default=8000)
    serve.set_defaults(handler=_serve)

    menu = subparsers.add_parser("menu", help="Open the guided ENZTRA terminal menu")
    menu.add_argument("--jobs-root", type=Path, default=Path("jobs"))
    menu.add_argument("--config", type=Path, default=Path("enztra.config.json"))
    menu.set_defaults(handler=_menu)

    prepare_design = subparsers.add_parser(
        "prepare-design", help="Prepare a design project through guided questions"
    )
    prepare_design.add_argument("--jobs-root", type=Path, default=Path("jobs"))
    prepare_design.set_defaults(handler=_prepare_design)

    boltz2 = subparsers.add_parser(
        "boltz2", help="Run Boltz-2 structural validation on strict survivors"
    )
    boltz2.add_argument("--config", type=Path, required=True)
    boltz2.add_argument("--job-dir", type=Path, required=True)
    boltz2.add_argument(
        "--msa-mode", choices=("single", "server"), default="single"
    )
    boltz2.add_argument("--diffusion-samples", type=int, default=1)
    boltz2.add_argument("--seed", type=int, default=42)
    boltz2.add_argument("--prepare-only", action="store_true")
    boltz2.set_defaults(handler=_boltz2)

    report = subparsers.add_parser(
        "boltz2-report",
        help="Regenerate tables and plots from completed Boltz-2 outputs",
    )
    report.add_argument("--job-dir", type=Path, required=True)
    report.set_defaults(handler=_boltz2_report)
    return parser


def main() -> int:
    if len(sys.argv) == 1:
        from .interactive import run_menu

        return run_menu()
    args = _parser().parse_args()
    return args.handler(args)


if __name__ == "__main__":
    raise SystemExit(main())
