"""Guarded LigandMPNN sequence design for completed RFdiffusion2 jobs."""

from __future__ import annotations

import csv
import json
import re
import shlex
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

from .config import EnztraConfig

Runner = Callable[..., subprocess.CompletedProcess]
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def _position_only(token: str) -> str:
    match = re.fullmatch(r"([A-Za-z0-9]):[A-Z]?(\d+)", token)
    if not match:
        raise ValueError(f"invalid mapped scaffold residue: {token}")
    return f"{match.group(1)}{match.group(2)}"


def validate_fixed_positions(
    job_dir: Path, expected_backbones: list[str] | None = None
) -> None:
    """Require LigandMPNN's fixed positions to equal the TRB-derived mapping."""

    mapping_path = job_dir / "outputs/functional_site_mapping.csv"
    if not mapping_path.is_file():
        raise FileNotFoundError("functional-site mapping CSV is missing")
    expected: dict[str, set[str]] = {}
    with mapping_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            expected.setdefault(row["backbone"], set()).add(
                _position_only(row["scaffold_residue"])
            )

    fixed_files = sorted(
        (job_dir / "outputs/rfdiffusion2/ligmpnn").glob("pdbs_position_fixed_*.jsonl")
    )
    if not fixed_files:
        raise FileNotFoundError("LigandMPNN fixed-position output is missing")
    actual: dict[str, set[str]] = {}
    for path in fixed_files:
        for line in path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            for pdb, positions in json.loads(line).items():
                actual[Path(pdb).stem] = set(positions)
    if expected_backbones is not None:
        selected = set(expected_backbones)
        expected = {key: value for key, value in expected.items() if key in selected}
        actual = {key: value for key, value in actual.items() if key in selected}
        missing = selected - set(expected)
        if missing:
            raise ValueError(
                "functional-site mapping is missing backbones: "
                + ", ".join(sorted(missing))
            )
    if expected != actual:
        raise ValueError(
            f"LigandMPNN fixed-position mismatch: expected {expected}, found {actual}"
        )


def _backbone_pairs(output_dir: Path) -> list[tuple[Path, Path]]:
    """Return completed backbone PDB/TRB pairs, excluding trajectories."""

    pairs: list[tuple[Path, Path]] = []
    for pdb in sorted(output_dir.glob("*.pdb")):
        trb = pdb.with_suffix(".trb")
        if trb.is_file():
            pairs.append((pdb, trb))
    return pairs


def prepare_ligandmpnn_execution(job_dir: Path, config: EnztraConfig) -> dict:
    """Create a reviewable LigandMPNN plan for completed backbones."""

    job_dir = job_dir.expanduser().resolve()
    request_path = job_dir / "design_request.json"
    status_path = job_dir / "status.json"
    for path in (request_path, status_path):
        if not path.is_file():
            raise FileNotFoundError(f"required prepared-job file is missing: {path}")
    request = json.loads(request_path.read_text(encoding="utf-8"))
    status = json.loads(status_path.read_text(encoding="utf-8"))
    if status.get("rfdiffusion2_returncode") != 0 or not status.get("rfdiffusion2_validated"):
        raise ValueError(
            "RFdiffusion2 protein-guidepost validation has not completed successfully"
        )

    image = config.rfdiffusion2_root / "rf_diffusion/exec/bakerlab_rf_diffusion_aa.sif"
    pipeline = config.rfdiffusion2_root / "rf_diffusion/benchmark/pipeline.py"
    if not image.is_file() or not pipeline.is_file():
        raise FileNotFoundError("RFdiffusion2 container or pipeline.py was not found")

    output_dir = job_dir / "outputs/rfdiffusion2"
    mapping_path = job_dir / "outputs/functional_site_mapping.csv"
    if not mapping_path.is_file():
        raise FileNotFoundError(
            "validated functional-site mapping is missing; rerun RFdiffusion2"
        )
    backbones = _backbone_pairs(output_dir)
    if not backbones:
        raise FileNotFoundError(
            f"no completed RFdiffusion2 backbone PDB/TRB pairs were found in {output_dir}"
        )
    sequences_per_backbone = int(request["sequences_per_backbone"])
    command = [
        config.apptainer_executable,
        "exec",
        "--nv",
        "--env",
        "MKL_THREADING_LAYER=GNU",
        str(image),
        str(pipeline),
        "--config-name=open_source_demo",
        f"outdir={output_dir}",
        "start_step=mpnn",
        "stop_step=thread_mpnn",
        "use_ligand=True",
        "in_proc=True",
        "mpnn.v2=True",
        "mpnn.ligand_present_for_all=True",
        f"mpnn.num_seq_per_target={sequences_per_backbone}",
        "mpnn.slurm.submit=True",
        "mpnn.slurm.in_proc=True",
        "mpnn.preprocessing_slurm.in_proc=True",
    ]
    plan = {
        "job_id": request["job_id"],
        "backbone_count": len(backbones),
        "sequences_per_backbone": sequences_per_backbone,
        "total_sequences": len(backbones) * sequences_per_backbone,
        "output_dir": str(output_dir),
        "backbones": [str(pdb) for pdb, _ in backbones],
        "backbone_ids": [pdb.stem for pdb, _ in backbones],
        "command": command,
        "command_display": shlex.join(command),
    }
    (job_dir / "ligandmpnn_execution.json").write_text(
        json.dumps(plan, indent=2) + "\n", encoding="utf-8"
    )
    return plan


def _read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = ""
    sequence: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header:
                records.append((header, "".join(sequence)))
            header, sequence = line[1:], []
        else:
            sequence.append(line)
    if header:
        records.append((header, "".join(sequence)))
    return records


def collect_ligandmpnn_sequences(
    job_dir: Path,
    expected_backbones: list[str],
    sequences_per_backbone: int,
) -> Path:
    """Create a clean, uniquely named FASTA from LigandMPNN output files."""

    job_dir = job_dir.expanduser().resolve()
    output_dir = job_dir / "outputs/rfdiffusion2"
    records: list[tuple[str, str]] = []
    for backbone in expected_backbones:
        fasta = output_dir / "ligmpnn" / "seqs" / f"{backbone}.fa"
        if not fasta.is_file():
            raise FileNotFoundError(
                f"LigandMPNN FASTA is missing for backbone {backbone}: {fasta}"
            )
        source_records = _read_fasta(fasta)
        # LigandMPNN normally writes the input/native record first.
        if len(source_records) == sequences_per_backbone + 1:
            source_records = source_records[1:]
        if len(source_records) != sequences_per_backbone:
            raise ValueError(
                "LigandMPNN sequence-count mismatch for "
                f"{backbone}: expected {sequences_per_backbone}, "
                f"found {len(source_records)}"
            )
        for _, sequence in source_records:
            sequence = sequence.replace("/", "").upper()
            if not sequence or not set(sequence) <= AMINO_ACIDS:
                invalid = ", ".join(sorted(set(sequence) - AMINO_ACIDS))
                raise ValueError(
                    f"LigandMPNN produced an invalid sequence for {backbone}"
                    + (f": {invalid}" if invalid else "")
                )
            records.append((backbone, sequence))

    destination = job_dir / "outputs/ligandmpnn/designs.fasta"
    destination.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    counters: dict[str, int] = {}
    for backbone, sequence in records:
        counters[backbone] = counters.get(backbone, 0) + 1
        safe_backbone = re.sub(r"[^A-Za-z0-9_.-]+", "_", backbone).strip("_")
        lines.extend(
            [f">{safe_backbone}_seq{counters[backbone]:04d}", sequence]
        )
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return destination


def run_ligandmpnn(
    job_dir: Path,
    config: EnztraConfig,
    runner: Runner = subprocess.run,
) -> int:
    """Run LigandMPNN and persist completion, failure, and FASTA location."""

    job_dir = job_dir.expanduser().resolve()
    plan = prepare_ligandmpnn_execution(job_dir, config)
    status_path = job_dir / "status.json"
    status = json.loads(status_path.read_text(encoding="utf-8"))
    status.update({"stage": "ligandmpnn_running", "ligandmpnn_started_at": datetime.now(timezone.utc).isoformat()})
    status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    try:
        completed = runner(plan["command"], cwd=config.rfdiffusion2_root, check=False)
        returncode = int(completed.returncode)
        designs_fasta = None
        if returncode == 0:
            try:
                validate_fixed_positions(job_dir, plan["backbone_ids"])
                designs_fasta = collect_ligandmpnn_sequences(
                    job_dir,
                    plan["backbone_ids"],
                    plan["sequences_per_backbone"],
                )
            except (ValueError, FileNotFoundError) as error:
                status["ligandmpnn_validation_error"] = str(error)
                returncode = 2
    except BaseException:
        status.update({"stage": "ligandmpnn_interrupted", "ligandmpnn_finished_at": datetime.now(timezone.utc).isoformat()})
        status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
        raise
    status.update(
        {
            "stage": "ligandmpnn_complete" if returncode == 0 else "ligandmpnn_failed",
            "ligandmpnn_returncode": returncode,
            "ligandmpnn_designs_fasta": str(designs_fasta) if designs_fasta else "",
            "ligandmpnn_sequence_count": plan["total_sequences"] if returncode == 0 else 0,
            "ligandmpnn_fixed_positions_validated": returncode == 0,
            "ligandmpnn_finished_at": datetime.now(timezone.utc).isoformat(),
        }
    )
    status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    return returncode
