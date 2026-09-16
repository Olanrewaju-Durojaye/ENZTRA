"""Adapters around the official DLKcat and CatPred command-line workflows."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import re
import shutil
import subprocess
from typing import Callable, Sequence

from .config import EnztraConfig
from .fasta import FastaRecord


Runner = Callable[..., subprocess.CompletedProcess[str]]


@dataclass(frozen=True, slots=True)
class KineticPrediction:
    design_id: str
    sequence: str
    kcat_s: float
    km_mm: float
    km_sd_total: float | None
    km_sd_aleatoric: float | None
    km_sd_epistemic: float | None


def _run(command: Sequence[str], cwd: Path, runner: Runner) -> None:
    runner(
        list(command),
        cwd=cwd,
        check=True,
        text=True,
    )


def prepare_inputs(
    records: Sequence[FastaRecord],
    substrate_name: str,
    substrate_smiles: str,
    job_dir: Path,
) -> tuple[Path, Path]:
    """Write equivalent batch inputs for DLKcat and CatPred."""

    job_dir.mkdir(parents=True, exist_ok=True)
    dlkcat_input = job_dir / "dlkcat_input.tsv"
    job_key = re.sub(r"[^A-Za-z0-9_.-]+", "_", job_dir.parent.name).strip("_")
    catpred_input = job_dir / f"{job_key or 'enztra_job'}_catpred.csv"

    with dlkcat_input.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["Substrate Name", "Substrate SMILES", "Protein Sequence"],
            delimiter="\t",
        )
        writer.writeheader()
        for record in records:
            writer.writerow(
                {
                    "Substrate Name": substrate_name,
                    "Substrate SMILES": substrate_smiles,
                    "Protein Sequence": record.sequence,
                }
            )

    with catpred_input.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["Substrate", "SMILES", "sequence", "pdbpath"],
        )
        writer.writeheader()
        for record in records:
            writer.writerow(
                {
                    "Substrate": substrate_name,
                    "SMILES": substrate_smiles,
                    "sequence": record.sequence,
                    "pdbpath": f"{record.record_id}.pdb",
                }
            )
    return dlkcat_input, catpred_input


def run_dlkcat(
    config: EnztraConfig,
    input_path: Path,
    destination: Path,
    runner: Runner = subprocess.run,
) -> Path:
    """Run official DLKcat and immediately preserve its fixed-name output."""

    command = [
        config.conda_executable,
        "run",
        "--no-capture-output",
        "-n",
        config.dlkcat_environment,
        "python",
        "prediction_for_input.py",
        str(input_path.resolve()),
    ]
    _run(command, config.dlkcat_example_dir, runner)
    source = config.dlkcat_example_dir / "output.tsv"
    if not source.is_file():
        raise FileNotFoundError(f"DLKcat did not create {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return destination


def run_catpred(
    config: EnztraConfig,
    input_path: Path,
    destination: Path,
    runner: Runner = subprocess.run,
) -> Path:
    """Run CatPred's official GPU batch entry point and preserve its output."""

    command = [
        config.conda_executable,
        "run",
        "--no-capture-output",
        "-n",
        config.catpred_environment,
        "python",
        "demo_run.py",
        "--parameter",
        "km",
        "--input_file",
        str(input_path.resolve()),
        "--use_gpu",
        "--checkpoint_dir",
        str(config.catpred_km_checkpoint),
    ]
    _run(command, config.catpred_root, runner)
    source = config.catpred_results_root / f"{input_path.stem}_input_output.csv"
    if not source.is_file():
        raise FileNotFoundError(f"CatPred did not create {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return destination


def _float_or_none(row: dict[str, str], key: str) -> float | None:
    value = row.get(key, "").strip()
    return float(value) if value else None


def combine_predictions(
    records: Sequence[FastaRecord],
    dlkcat_output: Path,
    catpred_output: Path,
) -> list[KineticPrediction]:
    """Combine tool outputs with strict row-count and identity validation."""

    with dlkcat_output.open(newline="", encoding="utf-8") as handle:
        dlkcat_rows = list(csv.DictReader(handle, delimiter="\t"))
    with catpred_output.open(newline="", encoding="utf-8") as handle:
        catpred_rows = list(csv.DictReader(handle))
    expected = len(records)
    if len(dlkcat_rows) != expected or len(catpred_rows) != expected:
        raise ValueError(
            "Prediction row-count mismatch: "
            f"expected={expected}, DLKcat={len(dlkcat_rows)}, CatPred={len(catpred_rows)}"
        )

    predictions: list[KineticPrediction] = []
    for record, dlkcat_row, catpred_row in zip(
        records, dlkcat_rows, catpred_rows, strict=True
    ):
        expected_key = f"{record.record_id}.pdb"
        if catpred_row.get("pdbpath") != expected_key:
            raise ValueError(
                f"CatPred row identity mismatch for {record.record_id}: "
                f"received {catpred_row.get('pdbpath')!r}"
            )
        if dlkcat_row.get("Protein Sequence") != record.sequence:
            raise ValueError(f"DLKcat sequence mismatch for {record.record_id}")
        predictions.append(
            KineticPrediction(
                design_id=record.record_id,
                sequence=record.sequence,
                kcat_s=float(dlkcat_row["Kcat value (1/s)"]),
                km_mm=float(catpred_row["Prediction_(mM)"]),
                km_sd_total=_float_or_none(catpred_row, "SD_total"),
                km_sd_aleatoric=_float_or_none(catpred_row, "SD_aleatoric"),
                km_sd_epistemic=_float_or_none(catpred_row, "SD_epistemic"),
            )
        )
    return predictions
