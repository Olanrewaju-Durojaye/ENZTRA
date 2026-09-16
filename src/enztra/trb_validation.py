"""Container-side validation of RFdiffusion2 protein-guidepost mappings."""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

import numpy as np


THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def _unwrap(value):
    if getattr(value, "shape", None) == ():
        return value.item()
    return value


def _reference_identities(request: dict) -> dict[tuple[str, int], str]:
    identities: dict[tuple[str, int], str] = {}
    for line in request.get("reference_content", "").splitlines():
        if line.startswith("ATOM  ") and len(line) > 26 and line[22:26].strip():
            identities[(line[21].strip(), int(line[22:26]))] = THREE_TO_ONE.get(
                line[17:20].strip().upper(), ""
            )
    for token in request.get("functional_site_residues", []):
        match = re.fullmatch(r"([A-Za-z0-9]):([A-Z])?(\d+)", token)
        if match and match.group(2):
            identities[(match.group(1), int(match.group(3)))] = match.group(2)
    return identities


def validate(
    directory: Path,
    expected_residues: int,
    expected_designs: int,
    request_path: Path | None = None,
    mapping_path: Path | None = None,
) -> None:
    trbs = sorted(directory.glob("*.trb"))
    if len(trbs) != expected_designs:
        raise ValueError(
            f"expected {expected_designs} backbone TRB files, found {len(trbs)}"
        )
    request = json.loads(request_path.read_text()) if request_path else {}
    identities = _reference_identities(request)
    rows: list[dict[str, str]] = []
    for path in trbs:
        data = np.load(path, allow_pickle=True)
        mapping = _unwrap(data["atomize_indices2atomname"])
        mapped = len(mapping)
        if mapped != expected_residues:
            raise ValueError(
                f"{path.name}: expected {expected_residues} mapped protein "
                f"guideposts, found {mapped}"
            )
        references = list(_unwrap(data["con_ref_pdb_idx"]))
        scaffolds = list(_unwrap(data["con_hal_pdb_idx"]))
        scaffold_indices = list(_unwrap(data["con_hal_idx0"]))
        if not (len(references) == len(scaffolds) == len(scaffold_indices) == mapped):
            raise ValueError(f"{path.name}: guidepost mapping arrays are inconsistent")
        for reference, scaffold, scaffold_idx in zip(
            references, scaffolds, scaffold_indices
        ):
            ref_chain, ref_number = str(reference[0]), int(reference[1])
            out_chain, out_number = str(scaffold[0]), int(scaffold[1])
            identity = identities.get((ref_chain, ref_number), "")
            atoms = mapping.get(int(scaffold_idx), [])
            rows.append({
                "backbone": path.stem,
                "reference_residue": f"{ref_chain}:{identity}{ref_number}",
                "scaffold_residue": f"{out_chain}:{identity}{out_number}",
                "fixed_atoms": ",".join(atoms),
            })
    if mapping_path:
        mapping_path.parent.mkdir(parents=True, exist_ok=True)
        with mapping_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=(
                "backbone", "reference_residue", "scaffold_residue", "fixed_atoms"
            ))
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    try:
        validate(
            Path(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]),
            Path(sys.argv[4]), Path(sys.argv[5]),
        )
    except Exception as error:
        print(f"RFdiffusion2 guidepost validation FAILED: {error}", file=sys.stderr)
        raise SystemExit(2)
    print(f"RFdiffusion2 guidepost validation passed. Mapping: {sys.argv[5]}")
