"""Installation paths and isolated environment names for external tools."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path


def _path(value: str) -> Path:
    return Path(value).expanduser().resolve()


@dataclass(frozen=True, slots=True)
class EnztraConfig:
    conda_executable: str
    apptainer_executable: str
    rfdiffusion2_root: Path
    dlkcat_root: Path
    catpred_root: Path
    catpred_data_root: Path
    catpred_results_root: Path
    boltz_cache: Path
    dlkcat_environment: str = "dlkcat"
    catpred_environment: str = "catpred"
    boltz2_environment: str = "boltz2"

    @classmethod
    def from_json(cls, path: Path) -> "EnztraConfig":
        data = json.loads(path.read_text(encoding="utf-8"))
        return cls(
            conda_executable=data.get("conda_executable", "conda"),
            apptainer_executable=data.get("apptainer_executable", "apptainer"),
            rfdiffusion2_root=_path(data["rfdiffusion2_root"]),
            dlkcat_root=_path(data["dlkcat_root"]),
            catpred_root=_path(data["catpred_root"]),
            catpred_data_root=_path(data["catpred_data_root"]),
            catpred_results_root=_path(
                data.get(
                    "catpred_results_root",
                    str(Path(data["catpred_root"]).expanduser().parent / "results"),
                )
            ),
            boltz_cache=_path(data.get("boltz_cache", "~/.boltz")),
            dlkcat_environment=data.get("dlkcat_environment", "dlkcat"),
            catpred_environment=data.get("catpred_environment", "catpred"),
            boltz2_environment=data.get("boltz2_environment", "boltz2"),
        )

    @property
    def dlkcat_example_dir(self) -> Path:
        return self.dlkcat_root / "DeeplearningApproach" / "Code" / "example"

    @property
    def catpred_km_checkpoint(self) -> Path:
        return self.catpred_data_root / "pretrained" / "production" / "km"
