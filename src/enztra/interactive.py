"""BOLTRA-style guided terminal workflow for ENZTRA."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Callable

from .design import (
    JOB_ID,
    create_design_job,
    list_design_jobs,
    parse_design_request,
    parse_design_length,
    parse_functional_site_residues,
)


def _ask(prompt: str, validate: Callable[[str], str], default: str | None = None) -> str:
    while True:
        suffix = f" [{default}]" if default is not None else ""
        answer = input(f"{prompt}{suffix}: ").strip()
        if not answer and default is not None:
            answer = default
        try:
            return validate(answer)
        except ValueError as error:
            print(f"Invalid entry: {error}")


def _not_empty(value: str) -> str:
    if not value:
        raise ValueError("a value is required")
    return value


def _path(value: str) -> str:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise ValueError(f"file not found: {path}")
    if path.suffix.lower() not in {".pdb", ".fasta", ".fa", ".faa"}:
        raise ValueError("use a .pdb, .fasta, .fa, or .faa file")
    if path.stat().st_size > 20_000_000:
        raise ValueError("reference file exceeds 20 MB")
    return str(path)


def _sequence_path(value: str) -> str:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise ValueError(f"file not found: {path}")
    if path.suffix.lower() not in {".fasta", ".fa", ".faa", ".txt"}:
        raise ValueError("use a .fasta, .fa, .faa, or FASTA-formatted .txt file")
    if path.stat().st_size > 100_000_000:
        raise ValueError("sequence file exceeds 100 MB")
    return str(path)


def _one_character(value: str) -> str:
    value = value.upper()
    if len(value) != 1 or not value.isalnum():
        raise ValueError("enter one chain letter or number")
    return value


def _job_id(value: str) -> str:
    value = value.lower().replace(" ", "-")
    if not JOB_ID.fullmatch(value):
        raise ValueError(
            "use 3–64 lowercase letters, numbers, hyphens, or underscores"
        )
    return value


def _ligand_code(value: str) -> str:
    value = value.upper()
    if not 1 <= len(value) <= 3 or not value.isalnum():
        raise ValueError("enter 1–3 letters or numbers")
    return value


def _substrate_ccd_code(generic_boltz_ligand: bool) -> Callable[[str], str]:
    def validate(value: str) -> str:
        code = _ligand_code(value)
        if generic_boltz_ligand and code == "LIG":
            raise ValueError(
                "LIG is Boltz's generic PDB label; enter the substrate's original "
                "CCD code, for example IHP"
            )
        return code

    return validate


def _residues(chain: str) -> Callable[[str], str]:
    def validate(value: str) -> str:
        return ", ".join(parse_functional_site_residues(value, chain))

    return validate


def _pdb_chain(chains: list[str]) -> Callable[[str], str]:
    def validate(value: str) -> str:
        chain = _one_character(value)
        if chain not in chains:
            raise ValueError(f"choose one of: {', '.join(chains)}")
        return chain

    return validate


def _integer(minimum: int, maximum: int) -> Callable[[str], str]:
    def validate(value: str) -> str:
        try:
            number = int(value)
        except ValueError as error:
            raise ValueError("enter a whole number") from error
        if not minimum <= number <= maximum:
            raise ValueError(f"enter a number from {minimum} to {maximum}")
        return str(number)

    return validate


def _design_length(reference_length: int) -> Callable[[str], str]:
    def validate(value: str) -> str:
        return parse_design_length(value, reference_length)

    return validate


def _pdb_atoms_for_positions(content: str, chain: str) -> dict[str, list[str]]:
    atoms: dict[str, set[str]] = {}
    for line in content.splitlines():
        if (
            line.startswith("ATOM  ")
            and len(line) > 26
            and line[21].strip() == chain
            and line[22:26].strip()
        ):
            position = f"{chain}:{line[22:26].strip()}{line[26].strip()}"
            atom = line[12:16].strip().upper()
            if atom and not atom.startswith("H"):
                atoms.setdefault(position, set()).add(atom)
    return {position: sorted(names) for position, names in atoms.items()}


def _atom_selection(available: list[str]) -> Callable[[str], str]:
    def validate(value: str) -> str:
        names = [name.strip().upper() for name in value.split(",") if name.strip()]
        if names == ["ALL"]:
            return ",".join(available)
        if not names:
            raise ValueError("enter one or more atom names, or ALL")
        unknown = sorted(set(names) - set(available))
        if unknown:
            raise ValueError(
                f"unknown atom(s): {', '.join(unknown)}; choose from: "
                + ", ".join(available)
            )
        return ",".join(dict.fromkeys(names))

    return validate


def _atom_preservation_mode(value: str) -> str:
    if value not in {"1", "2"}:
        raise ValueError("enter 1 or 2")
    return value


def _menu_choice(choices: set[str]) -> Callable[[str], str]:
    def validate(value: str) -> str:
        if value not in choices:
            raise ValueError("enter " + " or ".join(sorted(choices)))
        return value

    return validate


def _coordinates(value: str) -> str:
    parts = value.replace(",", " ").split()
    if len(parts) != 3:
        raise ValueError("enter exactly three coordinates: X Y Z")
    try:
        tuple(float(part) for part in parts)
    except ValueError as error:
        raise ValueError("coordinates must be numbers") from error
    return " ".join(parts)


def _existing_ori(content: str) -> bool:
    return any(
        line.startswith("HETATM")
        and line[12:16].strip().upper() == "ORI"
        and line[17:20].strip().upper() == "ORI"
        for line in content.splitlines()
    )


def _confirm(prompt: str, default: bool = False) -> bool:
    marker = "Y/n" if default else "y/N"
    while True:
        answer = input(f"{prompt} [{marker}]: ").strip().lower()
        if not answer:
            return default
        if answer in {"y", "yes"}:
            return True
        if answer in {"n", "no"}:
            return False
        print("Please enter y or n.")


def _pdb_chains(content: str) -> list[str]:
    return sorted(
        {
            line[21].strip()
            for line in content.splitlines()
            if line.startswith("ATOM  ")
            and len(line) > 21
            and line[21].strip()
        }
    )


def _generic_ligand_selectors(content: str) -> list[str]:
    return sorted(
        {
            f"{line[21].strip() or '_'}:{line[22:26].strip()}{line[26].strip()}"
            for line in content.splitlines()
            if line.startswith("HETATM")
            and len(line) > 26
            and line[17:20].strip().upper() == "LIG"
            and line[22:26].strip()
        }
    )


def _ligand_selector(choices: list[str]) -> Callable[[str], str]:
    def validate(value: str) -> str:
        selector = value.strip().upper()
        if selector not in choices:
            raise ValueError(f"choose one of: {', '.join(choices)}")
        return selector

    return validate


def prepare_design_interactively(jobs_root: Path = Path("jobs")) -> int:
    """Ask one question at a time and save a validated design request."""

    print("\nPrepare a new enzyme-redesign project")
    print("--------------------------------------")
    print("This stage validates and saves the project. It does not run a model yet.\n")
    print("IMPORTANT: Experimental PDB files may contain inactive mutations,")
    print("engineered substitutions, missing residues, or unresolved loops.")
    print("Use the sequence and structure representing the intended active enzyme.\n")
    reference_path = Path(_ask("Reference PDB or FASTA path", _path))
    reference_content = reference_path.read_text(encoding="utf-8")
    reference_format = "pdb" if reference_path.suffix.lower() == ".pdb" else "fasta"
    if reference_format == "pdb":
        chains = _pdb_chains(reference_content)
        print("Detected chains:", ", ".join(chains) if chains else "none")
        chain_default = chains[0] if chains else "A"
    else:
        chain_default = "A"

    job_id = _ask("Project name", _job_id)
    reference_id = _ask("Reference enzyme ID", _not_empty, reference_path.stem)
    chain_validator = (
        _pdb_chain(chains) if reference_format == "pdb" and chains else _one_character
    )
    protein_chain = _ask("Protein chain", chain_validator, chain_default)
    if reference_format == "fasta":
        reference_length = len(
            "".join(
                line.strip()
                for line in reference_content.splitlines()
                if line.strip() and not line.startswith(">")
            )
        )
    else:
        reference_length = len(
            {
                (line[22:26].strip(), line[26].strip())
                for line in reference_content.splitlines()
                if line.startswith("ATOM  ")
                and len(line) > 26
                and line[21].strip() == protein_chain
                and line[22:26].strip()
            }
        )
    substrate_name = _ask("Descriptive substrate name", _not_empty)
    generic_ligands = (
        _generic_ligand_selectors(reference_content)
        if reference_format == "pdb"
        else []
    )
    if generic_ligands:
        print(
            "Detected Boltz generic ligand label LIG at: "
            + ", ".join(generic_ligands)
        )
        print("Enter the substrate's original scientific CCD code below, not LIG.")
    ligand_code = _ask(
        "Substrate CCD code (for example IHP)",
        _substrate_ccd_code(bool(generic_ligands)),
    )
    structure_ligand_selector = ""
    if reference_format == "pdb":
        if len(generic_ligands) > 1:
            print("Multiple generic LIG residues detected:", ", ".join(generic_ligands))
            structure_ligand_selector = _ask(
                "Structure ligand chain/residue", _ligand_selector(generic_ligands)
            )
    substrate_smiles = _ask("Substrate stereochemical SMILES", _not_empty)
    functional_site_residues = _ask(
        "Functional-site residues to preserve (for example A:H17, A:R20)",
        _residues(protein_chain),
    )
    functional_site_atoms: dict[str, str] = {}
    atom_preservation_mode = "not_applicable"
    if reference_format == "pdb":
        print("\nHow should functional-site atoms be preserved?")
        print("1. Preserve all available heavy atoms for every selected residue")
        print("2. Choose atoms separately for each residue")
        mode = _ask("Select an atom-preservation mode", _atom_preservation_mode, "1")
        atom_preservation_mode = "all_heavy" if mode == "1" else "custom"
        pdb_atoms = _pdb_atoms_for_positions(reference_content, protein_chain)
        for token in functional_site_residues.split(","):
            token = token.strip()
            match = re.fullmatch(
                r"([A-Za-z0-9]):(?:[ACDEFGHIKLMNPQRSTVWY])?([1-9][0-9]*[A-Za-z]?)",
                token,
            )
            position = f"{match.group(1)}:{match.group(2)}"
            available = pdb_atoms[position]
            if mode == "1":
                functional_site_atoms[position] = ",".join(available)
            else:
                print(f"{token} available atoms: {', '.join(available)}")
                functional_site_atoms[position] = _ask(
                    f"Fixed atoms for {token}", _atom_selection(available)
                )
        if mode == "1":
            print(
                "Selected all available heavy atoms for "
                f"{len(functional_site_atoms)} functional-site residues."
            )
        print("\nWhere should RFdiffusion2 place the scaffold-centre ORI token?")
        has_existing_ori = _existing_ori(reference_content)
        if has_existing_ori:
            print("1. Use the ORI token already present in the PDB")
            print("2. Place ORI automatically at the protein C-alpha centroid")
            print("3. Enter custom X Y Z coordinates")
            ori_choice = _ask("Select an ORI-placement mode", _menu_choice({"1", "2", "3"}), "1")
            ori_mode = {"1": "existing", "2": "ca_centroid", "3": "custom"}[ori_choice]
        else:
            print("1. Place ORI automatically at the protein C-alpha centroid")
            print("2. Enter custom X Y Z coordinates")
            ori_choice = _ask("Select an ORI-placement mode", _menu_choice({"1", "2"}), "1")
            ori_mode = "ca_centroid" if ori_choice == "1" else "custom"
        ori_coordinates = (
            _ask("ORI coordinates (X Y Z)", _coordinates)
            if ori_mode == "custom"
            else ""
        )
    else:
        ori_mode = "not_applicable"
        ori_coordinates = ""
    design_length = _ask(
        "Designed protein length or range (for example 200 or 180-220)",
        _design_length(reference_length),
        str(reference_length),
    )
    backbone_count = _ask("Number of RFdiffusion2 backbones", _integer(1, 1000), "10")
    sequences_per_backbone = _ask(
        "LigandMPNN sequences per backbone", _integer(1, 100), "4"
    )
    notes = input("Optional project notes [leave blank to skip]: ").strip()

    raw = {
        "job_id": job_id,
        "reference_id": reference_id,
        "reference_format": reference_format,
        "reference_content": reference_content,
        "substrate_name": substrate_name,
        "substrate_smiles": substrate_smiles,
        "ligand_code": ligand_code,
        "structure_ligand_selector": structure_ligand_selector,
        "functional_site_residues": functional_site_residues,
        "functional_site_atoms": functional_site_atoms,
        "atom_preservation_mode": atom_preservation_mode,
        "ori_mode": ori_mode,
        "ori_coordinates": ori_coordinates,
        "protein_chain": protein_chain,
        "design_length": design_length,
        "backbone_count": backbone_count,
        "sequences_per_backbone": sequences_per_backbone,
        "notes": notes,
    }
    try:
        request = parse_design_request(raw)
    except ValueError as error:
        print(f"\nProject validation failed: {error}")
        print("No files were written.")
        return 1

    print("\nProject review")
    print("--------------")
    print(f"Project:                {request.job_id}")
    print(f"Reference:              {request.reference_id} ({request.reference_format.upper()})")
    print(f"Protein chain:          {request.protein_chain}")
    print(f"Substrate:              {request.substrate_name} ({request.ligand_code})")
    if request.structure_ligand_selector:
        print(
            "Structure mapping:       "
            f"LIG {request.structure_ligand_selector} → CCD {request.ligand_code}"
        )
    print(
        "Functional-site residues: "
        + ", ".join(request.functional_site_residues)
    )
    print(f"Designed protein length: {request.design_length} residues")
    if request.functional_site_atoms:
        if request.atom_preservation_mode == "all_heavy":
            print("Atom-preservation mode: All available heavy atoms")
            print(
                "Residues covered:       "
                f"{len(request.functional_site_atoms)}/"
                f"{len(request.functional_site_residues)}"
            )
        else:
            print("Atom-preservation mode: Custom atoms per residue")
            print("Fixed functional atoms:")
            for position, atoms in request.functional_site_atoms.items():
                print(f"  {position}: {', '.join(atoms)}")
        ori_labels = {
            "existing": "Existing PDB ORI token",
            "ca_centroid": "Automatic protein C-alpha centroid",
            "custom": "Custom coordinates",
        }
        print(f"ORI placement:          {ori_labels[request.ori_mode]}")
        print(
            "ORI coordinates:        "
            + ", ".join(f"{coordinate:.3f}" for coordinate in request.ori_coordinates)
        )
    print(f"RFdiffusion2 backbones: {request.backbone_count}")
    print(f"Sequences per backbone: {request.sequences_per_backbone}")
    print(f"Total planned designs:  {request.total_designs}")
    if request.validation_warnings:
        print("\nREFERENCE WARNINGS")
        for warning in request.validation_warnings:
            print(f"WARNING: {warning}")
    if request.validation_notes:
        print("\nREFERENCE NOTES")
        for note in request.validation_notes:
            print(f"NOTE: {note}")
    if not _confirm(
        "Have you confirmed that the reference represents the intended active enzyme?"
    ):
        print("Cancelled. Correct or verify the reference before preparing the project.")
        return 0
    if not _confirm("Save this prepared project?"):
        print("Cancelled. No files were written.")
        return 0
    try:
        status = create_design_job(jobs_root, request)
    except FileExistsError as error:
        print(f"Could not save project: {error}")
        return 1
    print(f"\nPrepared project: {status['job_id']}")
    print(f"Saved under: {(jobs_root / request.job_id).resolve()}")
    print("Model execution has not started.")
    if status["execution_ready"]:
        print(
            "RFdiffusion2 preflight specification: "
            + str((jobs_root / request.job_id / "rfdiffusion2_spec.json").resolve())
        )
    return 0


def inspect_projects(jobs_root: Path = Path("jobs")) -> int:
    jobs = list_design_jobs(jobs_root)
    if not jobs:
        print("\nNo prepared design projects were found.")
        return 0
    print("\nPrepared design projects")
    print("------------------------")
    for index, job in enumerate(jobs, 1):
        print(
            f"{index}. {job['job_id']} — {job['stage']} — "
            f"{job['total_designs']} planned designs — {job['substrate_name']}"
        )
    return 0


def run_prepared_project_interactively(
    jobs_root: Path = Path("jobs"), config_path: Path = Path("enztra.config.json")
) -> int:
    from .config import EnztraConfig
    from .execution import prepare_rfdiffusion2_execution, run_rfdiffusion2

    if not config_path.is_file():
        print(f"\nConfiguration not found: {config_path.resolve()}")
        return 1
    eligible = [job for job in list_design_jobs(jobs_root) if job.get("execution_ready")]
    if not eligible:
        print("\nNo RFdiffusion2-ready projects were found.")
        return 0
    print("\nRFdiffusion2-ready projects")
    print("---------------------------")
    for index, job in enumerate(eligible, 1):
        print(f"{index}. {job['job_id']} — {job['total_designs']} total planned designs")
    choice = _ask("Select a project", _integer(1, len(eligible)))
    selected = eligible[int(choice) - 1]
    job_dir = jobs_root / selected["job_id"]
    config = EnztraConfig.from_json(config_path)
    try:
        plan = prepare_rfdiffusion2_execution(job_dir, config)
    except (FileNotFoundError, ValueError) as error:
        print(f"Cannot run project: {error}")
        return 1
    print("\nGPU execution review")
    print("--------------------")
    print(f"Project:            {plan['job_id']}")
    print(f"Backbones requested: {plan['backbone_count']}")
    print(f"Output directory:    {plan['output_dir']}")
    print("\nExact command:")
    print(plan["command_display"])
    print("\nThis may take several minutes per backbone and will occupy the GPU.")
    if not _confirm("Start RFdiffusion2 GPU generation now?"):
        print("Cancelled. The command preview was saved; no model was started.")
        return 0
    returncode = run_rfdiffusion2(job_dir, config)
    if returncode:
        print(f"RFdiffusion2 failed with exit code {returncode}. Project status was saved.")
        return returncode
    print(f"RFdiffusion2 backbone generation completed: {plan['output_dir']}")
    return 0


def run_ligandmpnn_interactively(
    jobs_root: Path = Path("jobs"), config_path: Path = Path("enztra.config.json")
) -> int:
    """Select a completed backbone job and run ligand-aware sequence design."""

    from .config import EnztraConfig
    from .ligandmpnn import prepare_ligandmpnn_execution, run_ligandmpnn

    if not config_path.is_file():
        print(f"\nConfiguration not found: {config_path.resolve()}")
        return 1
    eligible = [
        job for job in list_design_jobs(jobs_root)
        if job.get("rfdiffusion2_returncode") == 0
    ]
    if not eligible:
        print("\nNo projects with completed RFdiffusion2 backbones were found.")
        return 0
    print("\nLigandMPNN-ready projects")
    print("-------------------------")
    for index, job in enumerate(eligible, 1):
        print(f"{index}. {job['job_id']} — {job['total_designs']} planned sequences")
    choice = _ask("Select a project", _integer(1, len(eligible)))
    selected = eligible[int(choice) - 1]
    job_dir = jobs_root / selected["job_id"]
    config = EnztraConfig.from_json(config_path)
    try:
        plan = prepare_ligandmpnn_execution(job_dir, config)
    except (FileNotFoundError, ValueError) as error:
        print(f"Cannot run LigandMPNN: {error}")
        return 1
    print("\nLigandMPNN execution review")
    print("---------------------------")
    print(f"Project:                {plan['job_id']}")
    print(f"Completed backbones:    {plan['backbone_count']}")
    print(f"Sequences per backbone: {plan['sequences_per_backbone']}")
    print(f"Total sequences:        {plan['total_sequences']}")
    print(f"Output directory:       {plan['output_dir']}/ligmpnn")
    print("\nExact command:")
    print(plan["command_display"])
    print("\nFunctional-site positions remain fixed during sequence design.")
    if not _confirm("Start LigandMPNN sequence generation now?"):
        print("Cancelled. The command preview was saved; no model was started.")
        return 0
    try:
        returncode = run_ligandmpnn(job_dir, config)
    except (FileNotFoundError, ValueError) as error:
        print(f"LigandMPNN output validation failed: {error}")
        return 1
    if returncode:
        print(f"LigandMPNN failed with exit code {returncode}. Project status was saved.")
        return returncode
    status = json.loads((job_dir / "status.json").read_text(encoding="utf-8"))
    print(f"LigandMPNN sequence generation completed: {status['ligandmpnn_designs_fasta']}")
    return 0


def export_mapping_interactively(
    jobs_root: Path = Path("jobs"), config_path: Path = Path("enztra.config.json")
) -> int:
    """Export mappings from existing validated RFdiffusion2 TRB files."""

    from .config import EnztraConfig
    from .execution import export_existing_mapping

    if not config_path.is_file():
        print(f"\nConfiguration not found: {config_path.resolve()}")
        return 1
    eligible = [
        job for job in list_design_jobs(jobs_root)
        if job.get("rfdiffusion2_returncode") == 0
    ]
    if not eligible:
        print("\nNo completed RFdiffusion2 projects were found.")
        return 0
    print("\nCompleted RFdiffusion2 projects")
    print("--------------------------------")
    for index, job in enumerate(eligible, 1):
        print(f"{index}. {job['job_id']}")
    choice = _ask("Select a project", _integer(1, len(eligible)))
    job_dir = jobs_root / eligible[int(choice) - 1]["job_id"]
    try:
        mapping = export_existing_mapping(
            job_dir, EnztraConfig.from_json(config_path)
        )
    except (FileNotFoundError, ValueError) as error:
        print(f"Could not export mapping: {error}")
        return 1
    print(f"Functional-site mapping exported: {mapping}")
    return 0


def run_boltz2_interactively(
    jobs_root: Path = Path("jobs"), config_path: Path = Path("enztra.config.json")
) -> int:
    """Select a completed kinetic job and validate its strict survivors."""

    from .boltz2 import (
        find_boltz2_ready_jobs,
        prepare_boltz2_execution,
        run_boltz2,
    )
    from .config import EnztraConfig

    if not config_path.is_file():
        print(f"\nConfiguration not found: {config_path.resolve()}")
        return 1
    eligible = find_boltz2_ready_jobs(jobs_root)
    if not eligible:
        print("\nNo completed kinetic jobs with strict survivors were found.")
        print("Boltz-2 is intentionally restricted to the strict-survivor queue.")
        return 0
    print("\nBoltz-2-ready kinetic jobs")
    print("---------------------------")
    for index, job in enumerate(eligible, 1):
        print(
            f"{index}. {job['name']} — {job['survivor_count']} strict survivors"
        )
    choice = _ask("Select a job", _integer(1, len(eligible)))
    job_dir = Path(eligible[int(choice) - 1]["job_dir"])
    print("\nMSA mode")
    print("--------")
    print("1. Single-sequence mode (reliable locally, but lower expected accuracy)")
    print("2. ColabFold MSA server (better context, requires internet and may fail)")
    msa_choice = _ask("Select an MSA mode", _integer(1, 2), "1")
    msa_mode = "single" if msa_choice == "1" else "server"
    diffusion_samples = int(
        _ask("Number of Boltz-2 structure samples", _integer(1, 25), "1")
    )
    config = EnztraConfig.from_json(config_path)
    try:
        plan = prepare_boltz2_execution(
            job_dir, config, msa_mode, diffusion_samples, 42
        )
    except (FileNotFoundError, ValueError) as error:
        print(f"Cannot prepare Boltz-2: {error}")
        return 1
    print("\nBoltz-2 execution review")
    print("------------------------")
    print(f"Kinetic job:       {job_dir}")
    print(f"Strict survivors:  {plan['survivor_count']}")
    print(f"Structure samples: {plan['diffusion_samples']} per survivor")
    print(f"MSA mode:          {plan['msa_mode']}")
    print("Affinity head:     Not requested")
    print(f"Output directory:  {plan['results_dir']}")
    print("\nExact command:")
    print(plan["command_display"])
    print("\nBoltz-2 confidence is structural evidence, not experimental activity.")
    if not _confirm("Start Boltz-2 GPU prediction now?"):
        print("Cancelled. Inputs and command preview were saved; no model was started.")
        return 0
    returncode = run_boltz2(job_dir, config, msa_mode, diffusion_samples, 42)
    if returncode:
        print(f"Boltz-2 failed with exit code {returncode}. Status was saved.")
        return returncode
    print(f"Boltz-2 structural validation completed: {job_dir / 'boltz2/summary.json'}")
    return 0


def refresh_boltz2_reports_interactively(jobs_root: Path = Path("jobs")) -> int:
    """Regenerate reports from completed Boltz-2 files without using the GPU."""

    from .boltz2 import find_completed_boltz2_jobs, refresh_boltz2_reports

    eligible = find_completed_boltz2_jobs(jobs_root)
    if not eligible:
        print("\nNo completed Boltz-2 jobs were found.")
        return 0
    print("\nCompleted Boltz-2 jobs")
    print("-----------------------")
    for index, job in enumerate(eligible, 1):
        print(f"{index}. {job['name']} — {job['survivor_count']} validated designs")
    choice = _ask("Select a job", _integer(1, len(eligible)))
    job_dir = Path(eligible[int(choice) - 1]["job_dir"])
    try:
        best_csv, plot_path = refresh_boltz2_reports(job_dir)
    except (FileNotFoundError, ValueError, json.JSONDecodeError) as error:
        print(f"Could not refresh Boltz-2 reports: {error}")
        return 1
    print(f"Boltz-2 best-model table: {best_csv}")
    print(f"Boltz-2 confidence plot: {plot_path}")
    print("No prediction was rerun and the GPU was not used.")
    return 0


def run_kinetics_interactively(
    jobs_root: Path = Path("jobs"), config_path: Path = Path("enztra.config.json")
) -> int:
    """Import sequences from any FASTA source and prepare or run kinetics."""

    from argparse import Namespace

    from .cli import _kinetics
    from .design import pdb_chain_sequence
    from .fasta import FastaRecord, read_fasta
    from .kinetic_import import (
        normalize_external_records,
        safe_identifier,
        write_imported_kinetic_inputs,
    )

    if not config_path.is_file():
        print(f"\nConfiguration not found: {config_path.resolve()}")
        return 1
    print("\nPrepare kinetic triage from protein sequences")
    print("---------------------------------------------")
    print("Sequences may come from ENZTRA, BOLTRA, LigandMPNN, or another source.")
    print("1. Import an external FASTA or FASTA-formatted text file")
    print("2. Use a completed ENZTRA LigandMPNN project")
    source_mode = _ask("Select a sequence source", _menu_choice({"1", "2"}), "1")
    if source_mode == "1":
        source_path = Path(_ask("Design-sequence file", _sequence_path))
    else:
        eligible = []
        for job in list_design_jobs(jobs_root):
            path = jobs_root / job["job_id"] / "outputs/ligandmpnn/designs.fasta"
            if job.get("ligandmpnn_returncode") == 0 and path.is_file():
                eligible.append((job, path))
        if not eligible:
            print("No completed ENZTRA LigandMPNN sequence sets were found.")
            return 0
        for index, (job, _) in enumerate(eligible, 1):
            print(f"{index}. {job['job_id']} — {job['total_designs']} planned designs")
        selected = int(_ask("Select a project", _integer(1, len(eligible)))) - 1
        source_path = eligible[selected][1]

    reference_path = Path(_ask("Active reference PDB or FASTA path", _path))
    if reference_path.suffix.lower() == ".pdb":
        content = reference_path.read_text(encoding="utf-8")
        chains = _pdb_chains(content)
        if not chains:
            print("Reference PDB contains no protein chains.")
            return 1
        print("Detected chains:", ", ".join(chains))
        chain = _ask("Reference protein chain", _pdb_chain(chains), chains[0])
        try:
            sequence = pdb_chain_sequence(content, chain)
        except ValueError as error:
            print(f"Could not extract reference sequence: {error}")
            return 1
        reference_id = _ask("Reference enzyme ID", _not_empty, reference_path.stem)
        reference_id = safe_identifier(reference_id)
        reference = FastaRecord(reference_id, reference_id, sequence)
    else:
        try:
            references = read_fasta(reference_path)
        except ValueError as error:
            print(f"Invalid reference FASTA: {error}")
            return 1
        if len(references) != 1:
            print("Reference FASTA must contain exactly one protein.")
            return 1
        source_reference = references[0]
        reference_id = safe_identifier(source_reference.description)
        reference = FastaRecord(
            reference_id, source_reference.description, source_reference.sequence
        )

    print("\nIdentical sequence handling")
    print("---------------------------")
    print("1. Remove duplicate sequences and retain their provenance (recommended)")
    print("2. Keep every input record")
    duplicate_mode = _ask("Select duplicate handling", _menu_choice({"1", "2"}), "1")
    deduplicate = duplicate_mode == "1"
    try:
        designs, provenance = normalize_external_records(
            source_path, deduplicate, {reference.record_id}
        )
    except ValueError as error:
        print(f"Invalid design-sequence file: {error}")
        return 1
    duplicate_count = sum(bool(row["duplicate_of"]) for row in provenance)
    if not designs:
        print("No design sequences remain after import.")
        return 1

    job_id = _ask("Kinetic job name", _job_id)
    substrate_name = _ask("Descriptive substrate name", _not_empty)
    substrate_smiles = _ask("Substrate stereochemical SMILES", _not_empty)
    print("\nExecution mode")
    print("--------------")
    print("1. Prepare validated DLKcat/CatPred inputs only")
    print("2. Run DLKcat and CatPred now")
    execution_mode = _ask("Select an execution mode", _menu_choice({"1", "2"}), "1")
    job_dir = jobs_root.expanduser().resolve() / job_id
    print("\nKinetic import review")
    print("---------------------")
    print(f"Source records:       {len(provenance)}")
    print(f"Identical duplicates: {duplicate_count}")
    print(f"Designs to evaluate:  {len(designs)}")
    print(f"Reference:            {reference.record_id} ({len(reference.sequence)} residues)")
    print(f"Substrate:            {substrate_name}")
    print(f"Execution:            {'prepare only' if execution_mode == '1' else 'run now'}")
    print(f"Job directory:        {job_dir}")
    if not _confirm(
        "Have you confirmed that the reference represents the intended active enzyme?"
    ):
        print("Cancelled. Verify the active reference before kinetic comparison.")
        return 0
    if not _confirm("Create this kinetic job?"):
        print("Cancelled. No files were written.")
        return 0
    if job_dir.exists():
        print(f"Could not create job: directory already exists: {job_dir}")
        return 1
    job_dir.mkdir(parents=True)
    reference_fasta, designs_fasta = write_imported_kinetic_inputs(
        job_dir, reference, designs, provenance, source_path, deduplicate
    )
    args = Namespace(
        config=config_path,
        reference_fasta=reference_fasta,
        designs_fasta=designs_fasta,
        substrate_name=substrate_name,
        substrate_smiles=substrate_smiles,
        job_dir=job_dir,
        prepare_only=execution_mode == "1",
        reuse_raw=False,
    )
    try:
        returncode = _kinetics(args)
    except (FileNotFoundError, ValueError) as error:
        print(f"Kinetic workflow failed: {error}")
        return 1
    manifest_path = job_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["sequence_source"] = str(source_path.resolve())
    manifest["source_metadata"] = "imported_inputs/source_metadata.csv"
    manifest["duplicates_removed"] = deduplicate
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return returncode


def run_workflow_interactively(
    jobs_root: Path = Path("jobs"), config_path: Path = Path("enztra.config.json")
) -> int:
    """Select a prepared project and run or resume every downstream stage."""

    from .workflow import run_workflow, workflow_summary

    if not config_path.is_file():
        print(f"\nConfiguration not found: {config_path.resolve()}")
        return 1
    eligible = [job for job in list_design_jobs(jobs_root) if job.get("execution_ready")]
    if not eligible:
        print("\nNo RFdiffusion2-ready projects were found.")
        return 0
    print("\nRun or resume complete workflows")
    print("--------------------------------")
    for index, job in enumerate(eligible, 1):
        state_path = jobs_root / job["job_id"] / "workflow_state.json"
        workflow_status = "not started"
        if state_path.is_file():
            try:
                workflow_status = json.loads(
                    state_path.read_text(encoding="utf-8")
                ).get("workflow_status", "unknown")
            except (OSError, json.JSONDecodeError):
                workflow_status = "invalid state"
        print(f"{index}. {job['job_id']} — {workflow_status}")
    choice = _ask("Select a project", _integer(1, len(eligible)))
    job_dir = jobs_root / eligible[int(choice) - 1]["job_id"]

    print("\nExisting stage status")
    print("---------------------")
    try:
        rows = workflow_summary(job_dir)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print(f"Cannot read workflow state: {error}")
        return 1
    for name, status, attempts, message in rows:
        detail = f" — {message}" if message else ""
        print(f"{name:14} {status:11} attempts={attempts}{detail}")

    print("\nBoltz-2 settings for this workflow")
    print("1. Single-sequence mode (reliable locally, lower expected accuracy)")
    print("2. ColabFold MSA server (requires internet and may fail)")
    msa_choice = _ask("Select an MSA mode", _integer(1, 2), "1")
    samples = int(_ask("Boltz-2 structure samples per survivor", _integer(1, 25), "1"))
    print("\nCompleted stages with validated outputs will be skipped.")
    print("Execution stops at the first failure and can be resumed from this menu.")
    if not _confirm("Run or resume this complete workflow now?"):
        print("Cancelled. Existing outputs and workflow state were not changed.")
        return 0
    try:
        returncode = run_workflow(
            job_dir,
            config_path,
            msa_mode="single" if msa_choice == "1" else "server",
            diffusion_samples=samples,
            seed=42,
        )
    except KeyboardInterrupt:
        print("\nWorkflow interrupted. Progress was saved and can be resumed.")
        return 130
    if returncode:
        print(f"Workflow stopped with exit code {returncode}. Progress was saved.")
        return returncode
    print(f"Workflow completed: {job_dir}")
    return 0


def generate_publication_exports_interactively(
    jobs_root: Path = Path("jobs"),
) -> int:
    """Select a completed kinetic job and create publication exports."""

    from .publication import generate_publication_exports

    eligible = []
    if jobs_root.is_dir():
        for selection in jobs_root.rglob("selection.csv"):
            job_dir = selection.parent
            if (job_dir / "summary.json").is_file() and (job_dir / "manifest.json").is_file():
                eligible.append(job_dir)
    eligible.sort(key=lambda path: path.stat().st_mtime, reverse=True)
    if not eligible:
        print("\nNo completed kinetic jobs were found.")
        return 0
    print("\nPublication-ready projects")
    print("--------------------------")
    for index, job_dir in enumerate(eligible, 1):
        summary = json.loads((job_dir / "summary.json").read_text(encoding="utf-8"))
        print(
            f"{index}. {job_dir.relative_to(jobs_root)} — "
            f"{summary.get('survivor_count', 0)}/{summary.get('design_count', 0)} strict survivors"
        )
    choice = _ask("Select a project", _integer(1, len(eligible)))
    job_dir = eligible[int(choice) - 1]
    try:
        result = generate_publication_exports(job_dir)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print(f"Could not generate publication exports: {error}")
        return 1
    print(f"Publication exports created: {job_dir / 'publication'}")
    print(f"Candidates:             {result['candidate_count']}")
    print(f"Strict survivors:       {result['strict_survivor_count']}")
    print(f"Structurally ranked:    {result['structurally_ranked_count']}")
    counts = result["qc_counts"]
    print(f"Quality control:        {counts['PASS']} pass, {counts['WARN']} warning, {counts['FAIL']} fail")
    return 0


def run_controlled_pilot_interactively(
    jobs_root: Path = Path("jobs"), config_path: Path = Path("enztra.config.json")
) -> int:
    """Select and run/resume a prepared controlled-pilot project."""

    from .pilot import prepare_pilot, run_controlled_pilot

    if not config_path.is_file():
        print(f"\nConfiguration not found: {config_path.resolve()}")
        return 1
    eligible = [job for job in list_design_jobs(jobs_root) if job.get("execution_ready")]
    if not eligible:
        print("\nNo RFdiffusion2-ready projects were found.")
        return 0
    print("\nControlled native-pilot projects")
    print("--------------------------------")
    for index, job in enumerate(eligible, 1):
        print(f"{index}. {job['job_id']} — {job['total_designs']} planned designs")
    choice = _ask("Select a project", _integer(1, len(eligible)))
    job_dir = jobs_root / eligible[int(choice) - 1]["job_id"]
    try:
        preflight = prepare_pilot(job_dir, config_path)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print(f"Pilot preflight failed: {error}")
        return 1
    print("\nControlled pilot preflight")
    print("--------------------------")
    print(f"Backbones:             {preflight['planned_backbones']}")
    print(f"Sequences per backbone:{preflight['sequences_per_backbone']:>4}")
    print(f"Total designs:         {preflight['planned_sequences']:>4}")
    print(f"Free disk space:       {preflight['free_disk_gib']:.2f} GiB")
    print(f"Installations:         {preflight['installation_status'].upper()}")
    for warning in preflight["warnings"]:
        print(f"WARNING: {warning}")
    if not preflight["ready"]:
        print("Pilot is blocked until failed installation checks are corrected.")
        return 1
    samples = int(_ask("Boltz-2 samples per strict survivor", _integer(1, 25), "1"))
    if not _confirm("Run or resume this controlled pilot now?"):
        print("Pilot preflight was saved; model execution was not started.")
        return 0
    try:
        returncode = run_controlled_pilot(
            job_dir, config_path, msa_mode="single",
            diffusion_samples=samples, seed=42,
        )
    except KeyboardInterrupt:
        print("\nPilot interrupted. Workflow state was saved and can be resumed.")
        return 130
    report_path = job_dir / "pilot/integration_report.json"
    print(f"Pilot integration report: {report_path}")
    if returncode:
        print(f"Pilot stopped with exit code {returncode}.")
    return returncode


def run_release_demonstration_interactively(
    jobs_root: Path = Path("jobs"), config_path: Path = Path("enztra.config.json")
) -> int:
    """Run/resume the exact 50 × 10 release-qualification workflow."""

    from .demonstration import (
        prepare_release_demonstration,
        run_release_demonstration,
    )

    if not config_path.is_file():
        print(f"\nConfiguration not found: {config_path.resolve()}")
        return 1
    eligible = [
        job for job in list_design_jobs(jobs_root)
        if job.get("execution_ready")
        and job.get("backbone_count") == 50
        and job.get("sequences_per_backbone") == 10
        and job.get("total_designs") == 500
    ]
    if not eligible:
        print("\nNo exact 50-backbone × 10-sequence release project was found.")
        print("Prepare a project with 50 RFdiffusion2 backbones and 10 sequences each.")
        return 0
    print("\nFull 500-sequence release demonstrations")
    print("----------------------------------------")
    for index, job in enumerate(eligible, 1):
        print(f"{index}. {job['job_id']} — 50 backbones × 10 sequences")
    choice = _ask("Select a project", _integer(1, len(eligible)))
    job_dir = jobs_root / eligible[int(choice) - 1]["job_id"]
    calibration = jobs_root / "pilot-50"
    calibration_job = calibration if calibration.is_dir() else None
    try:
        preflight = prepare_release_demonstration(
            job_dir, config_path, calibration_job
        )
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print(f"Release preflight failed: {error}")
        return 1
    print("\nRelease demonstration preflight")
    print("-------------------------------")
    print("Plan:                 50 backbones × 10 sequences = 500 designs")
    print(f"Estimated output:     {preflight['estimated_output_gib']:.2f} GiB")
    print(f"Required free space:  {preflight['required_free_gib']:.2f} GiB")
    print(f"Available free space: {preflight['available_free_gib']:.2f} GiB")
    print(f"Installations:        {preflight['installation_status'].upper()}")
    for blocker in preflight["blockers"]:
        print(f"BLOCKER: {blocker}")
    if not preflight["ready"]:
        print("Release demonstration is blocked until preflight issues are corrected.")
        return 2
    samples = int(_ask("Boltz-2 samples per strict survivor", _integer(1, 25), "1"))
    print("Completed stages with validated outputs will be skipped when resumed.")
    if not _confirm("Run or resume the full release demonstration now?"):
        print("Release preflight was saved; model execution was not started.")
        return 0
    try:
        returncode = run_release_demonstration(
            job_dir,
            config_path,
            calibration_job=calibration_job,
            msa_mode="single",
            diffusion_samples=samples,
            seed=42,
        )
    except KeyboardInterrupt:
        print("\nRelease demonstration interrupted. Progress was saved and can be resumed.")
        return 130
    print(f"Release qualification report: {job_dir / 'release/qualification_report.json'}")
    return returncode


def run_menu(jobs_root: Path = Path("jobs"), config: Path = Path("enztra.config.json")) -> int:
    """Run the main ENZTRA menu until the user chooses Exit."""

    from .config import EnztraConfig
    from .doctor import run_doctor
    from .web import serve

    while True:
        print("\nENZTRA — Enzyme Redesign and Triage Automation")
        print("================================================")
        print("1. Prepare a new enzyme-redesign project")
        print("2. Run a prepared RFdiffusion2 project")
        print("3. Run LigandMPNN on completed backbones")
        print("4. Run kinetic triage on protein sequences")
        print("5. Run Boltz-2 on strict kinetic survivors")
        print("6. Generate or refresh Boltz-2 result reports")
        print("7. Run or resume the complete workflow")
        print("8. Run or resume a controlled native pilot")
        print("9. Run or resume the full 500-sequence release demonstration")
        print("10. Generate publication and quality-control exports")
        print("11. Export or verify a functional-site mapping")
        print("12. Inspect prepared projects")
        print("13. Open the results dashboard")
        print("14. Check local installations")
        print("15. Exit")
        choice = input("Select an option [1-15]: ").strip()
        if choice == "1":
            prepare_design_interactively(jobs_root)
        elif choice == "2":
            run_prepared_project_interactively(jobs_root, config)
        elif choice == "3":
            run_ligandmpnn_interactively(jobs_root, config)
        elif choice == "4":
            run_kinetics_interactively(jobs_root, config)
        elif choice == "5":
            run_boltz2_interactively(jobs_root, config)
        elif choice == "6":
            refresh_boltz2_reports_interactively(jobs_root)
        elif choice == "7":
            run_workflow_interactively(jobs_root, config)
        elif choice == "8":
            run_controlled_pilot_interactively(jobs_root, config)
        elif choice == "9":
            run_release_demonstration_interactively(jobs_root, config)
        elif choice == "10":
            generate_publication_exports_interactively(jobs_root)
        elif choice == "11":
            export_mapping_interactively(jobs_root, config)
        elif choice == "12":
            inspect_projects(jobs_root)
        elif choice == "13":
            print("Your browser address will be http://127.0.0.1:8000")
            serve(jobs_root)
        elif choice == "14":
            if not config.is_file():
                print(f"Configuration not found: {config.resolve()}")
                print("Copy enztra.config.example.json to enztra.config.json first.")
                continue
            report = run_doctor(EnztraConfig.from_json(config))
            for check in report["checks"]:
                print(f"[{'PASS' if check['ok'] else 'FAIL'}] {check['name']}: {check['detail']}")
            print(f"ENZTRA status: {report['status'].upper()}")
        elif choice == "15":
            print("Goodbye.")
            return 0
        else:
            print("Invalid choice. Please enter a number from 1 to 15.")
