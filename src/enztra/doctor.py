"""Read-only installation checks for ENZTRA's external scientific tools."""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import subprocess

from .config import EnztraConfig


def _command_exists(command: str) -> bool:
    return Path(command).expanduser().is_file() or shutil.which(command) is not None


def run_doctor(config: EnztraConfig) -> dict:
    checks: list[dict[str, object]] = []

    def add(name: str, ok: bool, detail: str) -> None:
        checks.append({"name": name, "ok": ok, "detail": detail})

    add("conda", _command_exists(config.conda_executable), config.conda_executable)
    add(
        "apptainer",
        _command_exists(config.apptainer_executable),
        config.apptainer_executable,
    )
    paths = {
        "RFdiffusion2 root": config.rfdiffusion2_root,
        "RFdiffusion2 container": config.rfdiffusion2_root
        / "rf_diffusion/exec/bakerlab_rf_diffusion_aa.sif",
        "DLKcat predictor": config.dlkcat_example_dir / "prediction_for_input.py",
        "CatPred predictor": config.catpred_root / "demo_run.py",
        "CatPred Km checkpoint": config.catpred_km_checkpoint,
        "CatPred processed results": config.catpred_results_root,
        "Boltz cache": config.boltz_cache,
    }
    for name, path in paths.items():
        add(name, path.exists(), str(path))

    environments: set[str] = set()
    if _command_exists(config.conda_executable):
        completed = subprocess.run(
            [config.conda_executable, "env", "list", "--json"],
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode == 0:
            environments = {Path(path).name for path in json.loads(completed.stdout)["envs"]}
    for environment in (
        config.dlkcat_environment,
        config.catpred_environment,
        config.boltz2_environment,
    ):
        add(f"Conda environment: {environment}", environment in environments, environment)

    nvidia_smi = shutil.which("nvidia-smi")
    if nvidia_smi is None:
        add("NVIDIA GPU", False, "nvidia-smi was not found")
    else:
        gpu = subprocess.run(
            [nvidia_smi, "--query-gpu=name,memory.total", "--format=csv,noheader"],
            check=False,
            capture_output=True,
            text=True,
        )
        detail = gpu.stdout.strip() or gpu.stderr.strip() or "GPU query failed"
        add("NVIDIA GPU", gpu.returncode == 0, detail)

    return {
        "status": "ready" if all(bool(check["ok"]) for check in checks) else "incomplete",
        "checks": checks,
    }
