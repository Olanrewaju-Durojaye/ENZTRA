"""Guarded local RFdiffusion2 backbone execution for prepared ENZTRA jobs."""

from __future__ import annotations

import json
import shlex
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Sequence

from .config import EnztraConfig

Runner = Callable[..., subprocess.CompletedProcess]


def _position(token: str) -> str:
    chain, residue = token.split(":", 1)
    residue = residue[1:] if residue[0].isalpha() else residue
    return f"{chain}{residue}"


def prepare_rfdiffusion2_execution(job_dir: Path, config: EnztraConfig) -> dict:
    """Create the official custom-benchmark input and a reviewable command."""

    job_dir = job_dir.expanduser().resolve()
    request_path = job_dir / "design_request.json"
    status_path = job_dir / "status.json"
    preflight_path = job_dir / "rfdiffusion2_spec.json"
    for path in (request_path, status_path, preflight_path):
        if not path.is_file():
            raise FileNotFoundError(f"required prepared-job file is missing: {path}")
    request = json.loads(request_path.read_text(encoding="utf-8"))
    status = json.loads(status_path.read_text(encoding="utf-8"))
    if not status.get("execution_ready"):
        raise ValueError("this project has not passed RFdiffusion2 preflight")

    image = config.rfdiffusion2_root / "rf_diffusion/exec/bakerlab_rf_diffusion_aa.sif"
    inference = config.rfdiffusion2_root / "rf_diffusion/run_inference.py"
    if not image.is_file() or not inference.is_file():
        raise FileNotFoundError("RFdiffusion2 container or run_inference.py was not found")

    input_pdb = job_dir / "inputs/rfdiffusion2_input.pdb"
    if not input_pdb.is_file():
        raise FileNotFoundError(f"cleaned RFdiffusion2 input is missing: {input_pdb}")
    positions = [_position(token) for token in request["functional_site_residues"]]
    # In guidepost mode RFdiffusion2 discards the individual inter-motif gaps
    # and retains their combined scaffold length. Put the requested total
    # length in one scaffold segment so ranges such as 150-180 remain valid.
    contigs: list[str] = [request["design_length"]]
    contigs.extend(f"{position}-{position[1:]}" for position in positions)
    atom_map = {
        position.replace(":", ""): ",".join(atoms)
        for position, atoms in request["functional_site_atoms"].items()
    }
    # This is passed as one argv item directly to Hydra.  Avoid the benchmark
    # wrapper's second shell-parsing layer, which can silently strip the atom
    # dictionary and therefore discard every protein guidepost.
    atom_literal = "{" + ",".join(
        f'{position}:"{atoms}"' for position, atoms in atom_map.items()
    ) + "}"
    structural_ligand = "LIG" if request.get("structure_ligand_selector") else request["ligand_code"]
    output_dir = job_dir / "outputs/rfdiffusion2"
    output_prefix = output_dir / f"run_{request['job_id']}_cond0"
    command = [
        config.apptainer_executable,
        "exec",
        "--nv",
        "--env",
        "MKL_THREADING_LAYER=GNU",
        str(image),
        str(inference),
        "--config-name=aa",
        "inference.deterministic=True",
        "inference.ckpt_path=REPO_ROOT/rf_diffusion/model_weights/RFD_173.pt",
        "inference.seed_offset=43",
        f"inference.input_pdb={input_pdb}",
        f"inference.output_prefix={output_prefix}",
        f"inference.ligand={structural_ligand}",
        # RFdiffusion2 expects a list containing one comma-delimited contig
        # string. Without the inner quotes Hydra creates one list item per
        # segment and ContigMap silently reads only the requested length.
        f'contigmap.contigs=["{",".join(contigs)}"]',
        f"contigmap.length={request['design_length']}",
        "inference.contig_as_guidepost=True",
        f"contigmap.contig_atoms={atom_literal}",
        "inference.write_trajectory=True",
        "inference.write_trb_indep=True",
        "inference.write_trb_trajectory=True",
        f"inference.num_designs={request['backbone_count']}",
        "inference.design_startnum=0",
        "hydra.job_logging.root.level=WARN",
    ]
    validator = Path(__file__).with_name("trb_validation.py")
    mapping_path = job_dir / "outputs/functional_site_mapping.csv"
    validation_command = [
        config.apptainer_executable, "exec", "--env",
        f"PYTHONPATH={config.rfdiffusion2_root}", str(image), "python",
        str(validator), str(output_dir), str(len(positions)),
        str(request["backbone_count"]), str(request_path), str(mapping_path),
    ]
    plan = {
        "job_id": request["job_id"],
        "backbone_count": request["backbone_count"],
        "output_dir": str(output_dir),
        "command": command,
        "command_display": shlex.join(command),
        "validation_command": validation_command,
        "functional_site_mapping": str(mapping_path),
    }
    (job_dir / "rfdiffusion2_execution.json").write_text(
        json.dumps(plan, indent=2) + "\n", encoding="utf-8"
    )
    return plan


def run_rfdiffusion2(
    job_dir: Path,
    config: EnztraConfig,
    runner: Runner = subprocess.run,
) -> int:
    """Run a reviewed plan locally and persist completion or failure state."""

    job_dir = job_dir.expanduser().resolve()
    plan = prepare_rfdiffusion2_execution(job_dir, config)
    status_path = job_dir / "status.json"
    status = json.loads(status_path.read_text(encoding="utf-8"))
    status.update({"stage": "rfdiffusion2_running", "started_at": datetime.now(timezone.utc).isoformat()})
    status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    try:
        completed = runner(plan["command"], cwd=config.rfdiffusion2_root, check=False)
        returncode = int(completed.returncode)
        validation_returncode = None
        if returncode == 0:
            validated = runner(
                plan["validation_command"], cwd=config.rfdiffusion2_root, check=False
            )
            validation_returncode = int(validated.returncode)
            if validation_returncode != 0:
                returncode = validation_returncode
    except BaseException:
        status.update({"stage": "rfdiffusion2_interrupted", "finished_at": datetime.now(timezone.utc).isoformat()})
        status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
        raise
    status.update(
        {
            "stage": "rfdiffusion2_complete" if returncode == 0 else "rfdiffusion2_failed",
            "rfdiffusion2_returncode": returncode,
            "rfdiffusion2_validated": returncode == 0,
            "rfdiffusion2_validation_returncode": validation_returncode,
            "functional_site_mapping": plan["functional_site_mapping"] if returncode == 0 else "",
            "rfdiffusion2_output_dir": plan["output_dir"],
            "finished_at": datetime.now(timezone.utc).isoformat(),
        }
    )
    status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    return returncode


def export_existing_mapping(
    job_dir: Path,
    config: EnztraConfig,
    runner: Runner = subprocess.run,
) -> Path:
    """Validate existing TRBs and export their functional-site mapping only."""

    job_dir = job_dir.expanduser().resolve()
    plan = prepare_rfdiffusion2_execution(job_dir, config)
    completed = runner(
        plan["validation_command"], cwd=config.rfdiffusion2_root, check=False
    )
    if int(completed.returncode) != 0:
        raise ValueError("existing RFdiffusion2 guidepost validation failed")
    mapping = Path(plan["functional_site_mapping"])
    if not mapping.is_file():
        raise FileNotFoundError("guidepost validation did not create the mapping CSV")
    status_path = job_dir / "status.json"
    status = json.loads(status_path.read_text(encoding="utf-8"))
    status.update({
        "rfdiffusion2_validated": True,
        "rfdiffusion2_validation_returncode": 0,
        "functional_site_mapping": str(mapping),
    })
    fixed_files = list(
        (job_dir / "outputs/rfdiffusion2/ligmpnn").glob(
            "pdbs_position_fixed_*.jsonl"
        )
    )
    if fixed_files:
        from .ligandmpnn import validate_fixed_positions

        validate_fixed_positions(job_dir)
        status["ligandmpnn_fixed_positions_validated"] = True
    status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    return mapping
