"""Validation and persistence for guided enzyme design requests."""

from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

JOB_ID = re.compile(r"^[a-z0-9][a-z0-9_-]{1,62}[a-z0-9]$")
RESIDUE = re.compile(
    r"^(?P<chain>[A-Za-z0-9]):(?P<amino_acid>[ACDEFGHIKLMNPQRSTVWY])?"
    r"(?P<number>[1-9][0-9]*)(?P<icode>[A-Za-z]?)$"
)
LIGAND_SELECTOR = re.compile(
    r"^(?P<chain>[A-Za-z0-9]):(?P<number>-?[0-9]+)(?P<icode>[A-Za-z]?)$"
)
DESIGN_LENGTH = re.compile(r"^(?P<minimum>[1-9][0-9]*)(?:-(?P<maximum>[1-9][0-9]*))?$")
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
}


def pdb_chain_sequence(content: str, chain: str) -> str:
    """Extract one standard protein sequence from PDB ATOM records."""

    chain = chain.strip().upper()
    residues: dict[tuple[int, str], str] = {}
    for line in content.splitlines():
        if (
            not line.startswith("ATOM  ")
            or len(line) <= 26
            or line[21].strip().upper() != chain
            or not line[22:26].strip()
        ):
            continue
        key = (int(line[22:26]), line[26].strip())
        identity = THREE_TO_ONE.get(line[17:20].strip().upper())
        if identity is None:
            raise ValueError(
                f"PDB chain {chain} contains unsupported residue "
                f"{line[17:20].strip()} at {line[22:27].strip()}"
            )
        residues.setdefault(key, identity)
    if not residues:
        raise ValueError(f"PDB contains no protein residues for chain {chain}")
    return "".join(residues[key] for key in sorted(residues))


@dataclass(frozen=True, slots=True)
class DesignRequest:
    job_id: str
    reference_id: str
    reference_format: str
    reference_content: str
    substrate_name: str
    substrate_smiles: str
    ligand_code: str
    functional_site_residues: tuple[str, ...]
    protein_chain: str
    design_length: str
    backbone_count: int
    sequences_per_backbone: int
    functional_site_atoms: dict[str, tuple[str, ...]] | None = None
    atom_preservation_mode: str = "not_applicable"
    ori_mode: str = "not_applicable"
    ori_coordinates: tuple[float, float, float] | None = None
    validation_warnings: tuple[str, ...] = ()
    validation_notes: tuple[str, ...] = ()
    structure_ligand_selector: str = ""
    notes: str = ""

    @property
    def total_designs(self) -> int:
        return self.backbone_count * self.sequences_per_backbone

    @property
    def catalytic_residues(self) -> tuple[str, ...]:
        """Compatibility alias for specifications created before v0.4.2."""

        return self.functional_site_residues


def _required(data: dict[str, Any], key: str, maximum: int = 10_000) -> str:
    value = str(data.get(key, "")).strip()
    if not value:
        raise ValueError(f"{key.replace('_', ' ')} is required")
    if len(value) > maximum:
        raise ValueError(f"{key.replace('_', ' ')} is too long")
    return value


def _positive_int(data: dict[str, Any], key: str, maximum: int) -> int:
    try:
        value = int(data.get(key))
    except (TypeError, ValueError) as error:
        raise ValueError(f"{key.replace('_', ' ')} must be a whole number") from error
    if not 1 <= value <= maximum:
        raise ValueError(f"{key.replace('_', ' ')} must be between 1 and {maximum}")
    return value


def parse_design_length(value: Any, reference_length: int) -> str:
    """Validate an exact designed-protein length or inclusive length range."""

    normalised = str(value or reference_length).strip().replace(" ", "")
    match = DESIGN_LENGTH.fullmatch(normalised)
    if not match:
        raise ValueError("designed protein length must look like 200 or 180-220")
    minimum = int(match.group("minimum"))
    maximum = int(match.group("maximum") or minimum)
    if minimum > maximum:
        raise ValueError("designed protein length range must run from lower to higher")
    if maximum > 2000:
        raise ValueError("designed protein length may not exceed 2000 residues")
    return str(minimum) if minimum == maximum else f"{minimum}-{maximum}"


def _fasta_sequence(content: str) -> str:
    lines = [line.strip() for line in content.splitlines() if line.strip()]
    if not lines or not lines[0].startswith(">"):
        raise ValueError("FASTA input must begin with a > header")
    if sum(line.startswith(">") for line in lines) != 1:
        raise ValueError("Reference FASTA must contain exactly one protein")
    sequence = "".join(line for line in lines[1:] if not line.startswith(">")).upper()
    invalid = sorted(set(sequence) - AMINO_ACIDS)
    if not sequence or invalid:
        detail = f": {', '.join(invalid)}" if invalid else ""
        raise ValueError(f"Reference FASTA contains invalid amino-acid symbols{detail}")
    return sequence


def _pdb_residues(content: str, chain: str) -> dict[str, str]:
    atoms = [line for line in content.splitlines() if line.startswith("ATOM  ")]
    if not atoms:
        raise ValueError("PDB input contains no protein ATOM records")
    chains = {line[21].strip() for line in atoms if len(line) > 21 and line[21].strip()}
    if chain not in chains:
        available = ", ".join(sorted(chains)) or "unnamed"
        raise ValueError(f"Protein chain {chain!r} was not found; available chains: {available}")
    residues: dict[str, str] = {}
    for line in atoms:
        if len(line) <= 26 or line[21].strip() != chain or not line[22:26].strip():
            continue
        position = f"{chain}:{line[22:26].strip()}{line[26].strip()}"
        residues[position] = THREE_TO_ONE.get(line[17:20].strip().upper(), "X")
    return residues


def _pdb_residue_atoms(content: str, chain: str) -> dict[str, set[str]]:
    atoms: dict[str, set[str]] = {}
    for line in content.splitlines():
        if (
            not line.startswith("ATOM  ")
            or len(line) <= 26
            or line[21].strip() != chain
            or not line[22:26].strip()
        ):
            continue
        position = f"{chain}:{line[22:26].strip()}{line[26].strip()}"
        atom = line[12:16].strip().upper()
        if atom and not atom.startswith("H"):
            atoms.setdefault(position, set()).add(atom)
    return atoms


def _residue_position(token: str) -> str:
    match = RESIDUE.fullmatch(token)
    return f"{match.group('chain')}:{match.group('number')}{match.group('icode')}"


def _parse_functional_site_atoms(
    value: Any,
    residues: tuple[str, ...],
    available: dict[str, set[str]] | None,
) -> dict[str, tuple[str, ...]] | None:
    if value in (None, {}, ""):
        return None
    if available is None or not isinstance(value, dict):
        raise ValueError("functional-site atom selections require a PDB reference")
    expected_positions = [_residue_position(token) for token in residues]
    expected_position_set = set(expected_positions)
    normalised_value = {str(key).strip().upper(): atoms for key, atoms in value.items()}
    supplied_positions = set(normalised_value)
    if supplied_positions != expected_position_set:
        raise ValueError("select atoms for every functional-site residue, and only those residues")
    parsed: dict[str, tuple[str, ...]] = {}
    for position in expected_positions:
        if position not in available:
            raise ValueError(f"functional-site residue was not found in the PDB: {position}")
        raw = normalised_value[position]
        if isinstance(raw, str):
            selected = tuple(
                atom.strip().upper() for atom in raw.split(",") if atom.strip()
            )
        else:
            selected = tuple(str(atom).strip().upper() for atom in raw if str(atom).strip())
        if selected == ("ALL",):
            selected = tuple(sorted(available[position]))
        if not selected:
            raise ValueError(f"select at least one fixed atom for {position}")
        unknown = sorted(set(selected) - available[position])
        if unknown:
            raise ValueError(
                f"unknown atom(s) for {position}: {', '.join(unknown)}; available: "
                + ", ".join(sorted(available[position]))
            )
        parsed[position] = tuple(dict.fromkeys(selected))
    return parsed


def _pdb_ligands(content: str) -> dict[str, set[str]]:
    """Return residue-name to unique chain/residue selectors for PDB ligands."""

    ligands: dict[str, set[str]] = {}
    for line in content.splitlines():
        if not line.startswith("HETATM") or len(line) <= 26:
            continue
        code = line[17:20].strip().upper()
        number = line[22:26].strip()
        if not code or not number:
            continue
        chain = line[21].strip() or "_"
        selector = f"{chain}:{number}{line[26].strip()}"
        ligands.setdefault(code, set()).add(selector)
    return ligands


def _pdb_ori_coordinates(content: str) -> tuple[float, float, float] | None:
    for line in content.splitlines():
        if (
            line.startswith("HETATM")
            and len(line) >= 54
            and line[12:16].strip().upper() == "ORI"
            and line[17:20].strip().upper() == "ORI"
        ):
            try:
                return (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError as error:
                raise ValueError("the existing ORI token has invalid coordinates") from error
    return None


def _ca_centroid(content: str, chain: str) -> tuple[float, float, float]:
    coordinates: list[tuple[float, float, float]] = []
    for line in content.splitlines():
        if (
            line.startswith("ATOM  ")
            and len(line) >= 54
            and line[21].strip() == chain
            and line[12:16].strip().upper() == "CA"
            and line[16].strip() in {"", "A"}
        ):
            try:
                coordinates.append(
                    (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                )
            except ValueError as error:
                raise ValueError("the PDB contains invalid C-alpha coordinates") from error
    if not coordinates:
        raise ValueError("the selected PDB chain contains no C-alpha coordinates")
    count = len(coordinates)
    return tuple(sum(point[index] for point in coordinates) / count for index in range(3))


def _parse_ori_coordinates(value: Any) -> tuple[float, float, float]:
    if isinstance(value, str):
        parts = value.replace(",", " ").split()
    else:
        try:
            parts = list(value)
        except TypeError as error:
            raise ValueError("ORI coordinates must contain X, Y, and Z") from error
    if len(parts) != 3:
        raise ValueError("ORI coordinates must contain exactly X, Y, and Z")
    try:
        coordinates = tuple(float(part) for part in parts)
    except (TypeError, ValueError) as error:
        raise ValueError("ORI coordinates must be numbers") from error
    if any(not (-10000 < coordinate < 10000) for coordinate in coordinates):
        raise ValueError("ORI coordinates are outside the supported PDB range")
    return coordinates


def _clean_rfdiffusion_pdb(request: DesignRequest) -> str:
    """Keep the selected protein, intended ligand, and one validated ORI token."""

    retained: list[str] = []
    ligand_selector = request.structure_ligand_selector
    structural_code = "LIG" if ligand_selector else request.ligand_code
    for line in request.reference_content.splitlines():
        if line.startswith("ATOM  ") and len(line) > 21 and line[21].strip() == request.protein_chain:
            retained.append(line)
        elif line.startswith("HETATM") and len(line) > 26:
            code = line[17:20].strip().upper()
            selector = f"{line[21].strip() or '_'}:{line[22:26].strip()}{line[26].strip()}"
            if code == structural_code and (not ligand_selector or selector == ligand_selector):
                retained.append(line)
    x, y, z = request.ori_coordinates
    serial = max(
        (int(line[6:11]) for line in retained if len(line) >= 11 and line[6:11].strip().isdigit()),
        default=0,
    ) + 1
    retained.append(
        f"HETATM{serial:5d}  ORI ORI Z   1    {x:8.3f}{y:8.3f}{z:8.3f}"
        "  0.00  0.00          ORI"
    )
    retained.extend(["TER", "END"])
    return "\n".join(retained) + "\n"


def _normalise_ligand_selector(value: Any) -> str:
    selector = str(value or "").strip().upper()
    if selector and not LIGAND_SELECTOR.fullmatch(selector):
        raise ValueError("structure ligand selector must look like B:1 or B:501")
    return selector


def parse_functional_site_residues(value: str, chain: str) -> tuple[str, ...]:
    """Parse identity-aware residues such as A:H17 and legacy A:17."""

    residues = tuple(token.strip().upper() for token in value.split(",") if token.strip())
    if not residues or any(not RESIDUE.fullmatch(token) for token in residues):
        raise ValueError("functional-site residues must look like A:H17, A:R20, or A:92")
    if any(RESIDUE.fullmatch(token).group("chain") != chain for token in residues):
        raise ValueError("every functional-site residue must use the selected protein chain")
    positions = [
        f"{match.group('chain')}:{match.group('number')}{match.group('icode')}"
        for token in residues
        for match in [RESIDUE.fullmatch(token)]
    ]
    if len(positions) != len(set(positions)):
        raise ValueError("functional-site residues must not contain duplicates")
    return residues


def _identity_warnings(
    residues: tuple[str, ...], sequence: str | None, pdb_residues: dict[str, str] | None
) -> list[str]:
    warnings: list[str] = []
    for token in residues:
        match = RESIDUE.fullmatch(token)
        expected = match.group("amino_acid")
        number = int(match.group("number"))
        position = f"{match.group('chain')}:{match.group('number')}{match.group('icode')}"
        if sequence is not None:
            if match.group("icode"):
                raise ValueError(
                    f"FASTA positions cannot contain PDB insertion codes: {token}"
                )
            if number > len(sequence):
                raise ValueError(f"functional-site residue exceeds the FASTA sequence length: {token}")
            observed = sequence[number - 1]
        else:
            if position not in pdb_residues:
                raise ValueError(f"functional-site residue was not found in the selected PDB chain: {token}")
            observed = pdb_residues[position]
        if expected and observed != expected:
            warnings.append(
                f"Expected {expected} at {position}, but the reference contains {observed}. "
                "Check for an engineered or inactive mutation."
            )
    return warnings


def parse_design_request(data: dict[str, Any]) -> DesignRequest:
    """Validate an untrusted request and return a typed specification."""

    job_id = _required(data, "job_id", 64).lower().replace(" ", "-")
    if not JOB_ID.fullmatch(job_id):
        raise ValueError("job ID must use 3–64 lowercase letters, numbers, hyphens, or underscores")
    reference_format = _required(data, "reference_format", 10).lower()
    if reference_format not in {"fasta", "pdb"}:
        raise ValueError("reference format must be fasta or pdb")
    reference_content = _required(data, "reference_content", 20_000_000)
    protein_chain = _required(data, "protein_chain", 1).upper()
    if not protein_chain.isalnum():
        raise ValueError("protein chain must be one letter or number")
    ligand_code = _required(data, "ligand_code", 3).upper()
    if not ligand_code.isalnum():
        raise ValueError("ligand code must contain 1–3 letters or numbers")

    if reference_format == "fasta":
        sequence = _fasta_sequence(reference_content)
        pdb_residues = None
        pdb_atoms = None
    else:
        sequence = None
        pdb_residues = _pdb_residues(reference_content, protein_chain)
        pdb_atoms = _pdb_residue_atoms(reference_content, protein_chain)
    reference_length = len(sequence) if sequence is not None else len(pdb_residues)

    residue_value = data.get("functional_site_residues")
    if residue_value is None:
        residue_value = data.get("catalytic_residues", "")
    residues = parse_functional_site_residues(str(residue_value), protein_chain)
    functional_site_atoms = _parse_functional_site_atoms(
        data.get("functional_site_atoms"), residues, pdb_atoms
    )
    atom_preservation_mode = str(
        data.get(
            "atom_preservation_mode",
            "custom" if functional_site_atoms else "not_applicable",
        )
    ).strip().lower()
    if atom_preservation_mode not in {"all_heavy", "custom", "not_applicable"}:
        raise ValueError("atom preservation mode must be all_heavy or custom")
    if functional_site_atoms and atom_preservation_mode == "not_applicable":
        raise ValueError("atom preservation mode is required for PDB atom selections")
    if not functional_site_atoms and atom_preservation_mode != "not_applicable":
        raise ValueError("atom preservation mode requires PDB atom selections")
    warnings = _identity_warnings(residues, sequence, pdb_residues)
    notes: list[str] = []
    structure_ligand_selector = _normalise_ligand_selector(
        data.get("structure_ligand_selector")
    )
    if reference_format == "pdb":
        pdb_ligands = _pdb_ligands(reference_content)
        if ligand_code == "LIG" and pdb_ligands.get("LIG"):
            raise ValueError(
                "LIG is Boltz's generic structure label, not the substrate's "
                "scientific CCD code; enter the original code (for example IHP)"
            )
        if ligand_code not in pdb_ligands:
            generic_ligands = sorted(pdb_ligands.get("LIG", set()))
            if len(generic_ligands) == 1:
                structure_ligand_selector = generic_ligands[0]
                notes.append(
                    f"The structure uses generic ligand label LIG at "
                    f"{structure_ligand_selector}; it is mapped to {ligand_code}. "
                    "The stereochemical SMILES remains the authoritative substrate identity."
                )
            elif len(generic_ligands) > 1:
                choices = ", ".join(generic_ligands)
                if not structure_ligand_selector:
                    raise ValueError(
                        "multiple LIG residues were found; provide a structure ligand "
                        f"selector ({choices})"
                    )
                if structure_ligand_selector not in generic_ligands:
                    raise ValueError(
                        "structure ligand selector does not identify a LIG residue; "
                        f"choose one of: {choices}"
                    )
                notes.append(
                    f"Generic ligand LIG at {structure_ligand_selector} is mapped to "
                    f"{ligand_code}. The stereochemical SMILES remains the authoritative "
                    "substrate identity."
                )
            else:
                warnings.append(
                    f"Ligand CCD code {ligand_code} was not found in PDB HETATM records. "
                    "Confirm that the intended substrate is present or supply an appropriate complex."
                )

    if reference_format == "pdb":
        existing_ori = _pdb_ori_coordinates(reference_content)
        ori_mode = str(data.get("ori_mode", "existing" if existing_ori else "ca_centroid")).strip().lower()
        if ori_mode == "existing":
            if existing_ori is None:
                raise ValueError("ORI mode existing was selected, but the PDB contains no ORI token")
            ori_coordinates = existing_ori
        elif ori_mode == "ca_centroid":
            ori_coordinates = _ca_centroid(reference_content, protein_chain)
        elif ori_mode == "custom":
            ori_coordinates = _parse_ori_coordinates(data.get("ori_coordinates"))
        else:
            raise ValueError("ORI mode must be existing, ca_centroid, or custom")
    else:
        ori_mode = "not_applicable"
        ori_coordinates = None

    request = DesignRequest(
        job_id=job_id,
        reference_id=_required(data, "reference_id", 100),
        reference_format=reference_format,
        reference_content=reference_content,
        substrate_name=_required(data, "substrate_name", 200),
        substrate_smiles=_required(data, "substrate_smiles", 5000),
        ligand_code=ligand_code,
        functional_site_residues=residues,
        protein_chain=protein_chain,
        design_length=parse_design_length(data.get("design_length"), reference_length),
        backbone_count=_positive_int(data, "backbone_count", 1000),
        sequences_per_backbone=_positive_int(data, "sequences_per_backbone", 100),
        functional_site_atoms=functional_site_atoms,
        atom_preservation_mode=atom_preservation_mode,
        ori_mode=ori_mode,
        ori_coordinates=ori_coordinates,
        validation_warnings=tuple(warnings),
        validation_notes=tuple(notes),
        structure_ligand_selector=structure_ligand_selector,
        notes=str(data.get("notes", "")).strip()[:5000],
    )
    if request.total_designs > 10_000:
        raise ValueError("the job may request at most 10,000 total sequences")
    return request


def create_design_job(jobs_root: Path, request: DesignRequest) -> dict[str, Any]:
    """Persist a validated design request without launching external models."""

    job_dir = jobs_root.resolve() / request.job_id
    try:
        job_dir.mkdir(parents=True, exist_ok=False)
    except FileExistsError as error:
        raise FileExistsError(f"job {request.job_id!r} already exists") from error
    inputs = job_dir / "inputs"
    inputs.mkdir()
    reference_name = f"reference.{request.reference_format}"
    (inputs / reference_name).write_text(request.reference_content + "\n", encoding="utf-8")
    specification = asdict(request)
    specification.pop("reference_content")
    specification["functional_site_residues"] = list(request.functional_site_residues)
    specification["validation_warnings"] = list(request.validation_warnings)
    specification["validation_notes"] = list(request.validation_notes)
    if request.functional_site_atoms is not None:
        specification["functional_site_atoms"] = {
            residue: list(atoms)
            for residue, atoms in request.functional_site_atoms.items()
        }
    specification["total_designs"] = request.total_designs
    specification["reference_file"] = f"inputs/{reference_name}"
    specification["created_at"] = datetime.now(timezone.utc).isoformat()
    (job_dir / "design_request.json").write_text(
        json.dumps(specification, indent=2) + "\n", encoding="utf-8"
    )
    execution_ready = request.reference_format == "pdb" and bool(request.functional_site_atoms)
    if execution_ready:
        cleaned_reference = inputs / "rfdiffusion2_input.pdb"
        cleaned_reference.write_text(_clean_rfdiffusion_pdb(request), encoding="utf-8")
        motif_residues = [position.replace(":", "") for position in request.functional_site_atoms]
        structural_ligand = "LIG" if request.structure_ligand_selector else request.ligand_code
        fixed_atoms = {
            position.replace(":", ""): ",".join(atoms)
            for position, atoms in request.functional_site_atoms.items()
        }
        fixed_atoms[structural_ligand] = ""
        rfdiffusion_spec = {
            request.job_id: {
                "input": str(cleaned_reference.resolve()),
                "ligand": structural_ligand,
                "unindex": ",".join(motif_residues),
                "length": request.design_length,
                "select_fixed_atoms": fixed_atoms,
                "allow_ligand_on_existing_chain": True,
            }
        }
        (job_dir / "rfdiffusion2_spec.json").write_text(
            json.dumps(rfdiffusion_spec, indent=2) + "\n", encoding="utf-8"
        )
    status = {
        "job_id": request.job_id,
        "stage": "prepared",
        "message": "Design request validated and saved; model execution has not started.",
        "total_designs": request.total_designs,
        "warning_count": len(request.validation_warnings),
        "execution_ready": execution_ready,
    }
    (job_dir / "status.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    return status


def list_design_jobs(jobs_root: Path) -> list[dict[str, Any]]:
    """List guided design requests, newest first."""

    if not jobs_root.is_dir():
        return []
    jobs = []
    for path in jobs_root.rglob("design_request.json"):
        try:
            request = json.loads(path.read_text(encoding="utf-8"))
            status = json.loads((path.parent / "status.json").read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        jobs.append({
            **status,
            "substrate_name": request.get("substrate_name", ""),
            "backbone_count": request.get("backbone_count", 0),
            "sequences_per_backbone": request.get("sequences_per_backbone", 0),
            "created_at": request.get("created_at", ""),
        })
    return sorted(jobs, key=lambda job: job["created_at"], reverse=True)
