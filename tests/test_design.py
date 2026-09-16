import json

import pytest

from enztra.design import create_design_job, parse_design_request


def request(**updates):
    data = {
        "job_id": "enzyme-test-01",
        "reference_id": "WT_A",
        "reference_format": "fasta",
        "reference_content": ">WT_A\nACDEFGHIKLMNPQRSTVWY",
        "substrate_name": "test substrate",
        "substrate_smiles": "CCO",
        "ligand_code": "LIG",
        "functional_site_residues": "A:C2, A:L10",
        "protein_chain": "A",
        "backbone_count": 10,
        "sequences_per_backbone": 4,
    }
    data.update(updates)
    return data


def test_parse_and_create_design_job(tmp_path):
    parsed = parse_design_request(request())
    status = create_design_job(tmp_path, parsed)
    specification = json.loads((tmp_path / "enzyme-test-01" / "design_request.json").read_text())
    assert status["stage"] == "prepared"
    assert specification["total_designs"] == 40
    assert specification["functional_site_residues"] == ["A:C2", "A:L10"]
    assert specification["design_length"] == "20"
    assert "reference_content" not in specification
    assert (tmp_path / "enzyme-test-01" / "inputs" / "reference.fasta").is_file()


@pytest.mark.parametrize(
    "updates, message",
    [
        ({"functional_site_residues": "58, 116"}, "must look like"),
        ({"functional_site_residues": "B:C2"}, "selected protein chain"),
        ({"reference_content": ">one\nAAAA\n>two\nCCCC"}, "exactly one"),
        ({"functional_site_residues": "A:21"}, "exceeds the FASTA"),
        ({"backbone_count": 1001}, "between 1 and 1000"),
        ({"design_length": "220-180"}, "lower to higher"),
        ({"design_length": "2000-2001"}, "may not exceed 2000"),
    ],
)
def test_invalid_requests_are_rejected(updates, message):
    with pytest.raises(ValueError, match=message):
        parse_design_request(request(**updates))


def test_existing_job_is_not_overwritten(tmp_path):
    parsed = parse_design_request(request())
    create_design_job(tmp_path, parsed)
    with pytest.raises(FileExistsError):
        create_design_job(tmp_path, parsed)


def test_exact_and_ranged_design_lengths_are_normalised():
    assert parse_design_request(request(design_length="200")).design_length == "200"
    assert parse_design_request(request(design_length=" 180 - 220 ")).design_length == "180-220"


def test_identity_mismatch_creates_inactive_mutant_warning():
    parsed = parse_design_request(request(functional_site_residues="A:H2"))
    assert len(parsed.validation_warnings) == 1
    assert "reference contains C" in parsed.validation_warnings[0]


def test_legacy_position_only_notation_remains_supported():
    parsed = parse_design_request(request(functional_site_residues="A:2, A:10"))
    assert parsed.functional_site_residues == ("A:2", "A:10")
    assert parsed.validation_warnings == ()


def test_pdb_checks_expected_identity_and_ligand_code():
    pdb = "\n".join(
        [
            "ATOM      1  CA  ALA A  17      10.000  10.000  10.000  1.00 20.00           C",
            "HETATM    2  P   IHP B 501      11.000  10.000  10.000  1.00 20.00           P",
        ]
    )
    parsed = parse_design_request(
        request(
            reference_format="pdb",
            reference_content=pdb,
            functional_site_residues="A:H17",
            ligand_code="IHP",
        )
    )
    assert len(parsed.validation_warnings) == 1
    assert "Expected H at A:17" in parsed.validation_warnings[0]


def test_pdb_warns_when_ligand_code_is_absent():
    pdb = "ATOM      1  CA  HIS A  17      10.000  10.000  10.000  1.00 20.00           C"
    parsed = parse_design_request(
        request(
            reference_format="pdb",
            reference_content=pdb,
            functional_site_residues="A:H17",
            ligand_code="IHP",
        )
    )
    assert len(parsed.validation_warnings) == 1
    assert "CCD code IHP was not found" in parsed.validation_warnings[0]


def test_boltz_generic_ligand_is_mapped_without_warning():
    pdb = "\n".join(
        [
            "ATOM      1  CA  HIS A  17      10.000  10.000  10.000  1.00 20.00           C",
            "HETATM    2  P   LIG B   1      11.000  10.000  10.000  1.00 20.00           P",
        ]
    )
    parsed = parse_design_request(
        request(
            reference_format="pdb",
            reference_content=pdb,
            functional_site_residues="A:H17",
            ligand_code="IHP",
        )
    )
    assert parsed.validation_warnings == ()
    assert parsed.structure_ligand_selector == "B:1"
    assert "generic ligand label LIG" in parsed.validation_notes[0]
    assert "mapped to IHP" in parsed.validation_notes[0]


def test_multiple_generic_ligands_require_a_valid_selector():
    pdb = "\n".join(
        [
            "ATOM      1  CA  HIS A  17      10.000  10.000  10.000  1.00 20.00           C",
            "HETATM    2  P   LIG B   1      11.000  10.000  10.000  1.00 20.00           P",
            "HETATM    3  P   LIG C   2      12.000  10.000  10.000  1.00 20.00           P",
        ]
    )
    data = request(
        reference_format="pdb",
        reference_content=pdb,
        functional_site_residues="A:H17",
        ligand_code="IHP",
    )
    with pytest.raises(ValueError, match="multiple LIG residues"):
        parse_design_request(data)
    parsed = parse_design_request({**data, "structure_ligand_selector": "C:2"})
    assert parsed.structure_ligand_selector == "C:2"
    assert parsed.validation_warnings == ()


def test_boltz_generic_label_cannot_replace_scientific_ccd_code():
    pdb = "\n".join(
        [
            "ATOM      1  CA  HIS A  17      10.000  10.000  10.000  1.00 20.00           C",
            "HETATM    2  P   LIG B   1      11.000  10.000  10.000  1.00 20.00           P",
        ]
    )
    with pytest.raises(ValueError, match="generic structure label"):
        parse_design_request(
            request(
                reference_format="pdb",
                reference_content=pdb,
                functional_site_residues="A:H17",
                ligand_code="LIG",
            )
        )


def test_pdb_atom_selection_creates_rfdiffusion2_preflight(tmp_path):
    pdb = "\n".join(
        [
            "ATOM      1  ND1 HIS A  17      10.000  10.000  10.000  1.00 20.00           N",
            "ATOM      2  NE2 HIS A  17      11.000  10.000  10.000  1.00 20.00           N",
            "HETATM    3  P   LIG B   1      12.000  10.000  10.000  1.00 20.00           P",
        ]
    )
    parsed = parse_design_request(
        request(
            reference_format="pdb",
            reference_content=pdb,
            functional_site_residues="A:H17",
            functional_site_atoms={"A:17": "ND1,NE2"},
            ori_mode="custom",
            ori_coordinates="1 2 3",
            ligand_code="IHP",
            design_length="150-180",
        )
    )
    status = create_design_job(tmp_path, parsed)
    spec = json.loads((tmp_path / "enzyme-test-01/rfdiffusion2_spec.json").read_text())
    assert status["execution_ready"] is True
    assert spec["enzyme-test-01"]["ligand"] == "LIG"
    assert spec["enzyme-test-01"]["length"] == "150-180"
    assert spec["enzyme-test-01"]["select_fixed_atoms"]["A17"] == "ND1,NE2"
    cleaned = (tmp_path / "enzyme-test-01/inputs/rfdiffusion2_input.pdb").read_text()
    assert " LIG B   1" in cleaned
    assert " ORI ORI Z   1" in cleaned


def test_functional_atom_mapping_preserves_residue_input_order():
    pdb = "\n".join(
        [
            "ATOM      1  CA  LYS A  24      10.000  10.000  10.000  1.00 20.00           C",
            "ATOM      2  NZ  LYS A  24      10.500  10.000  10.000  1.00 20.00           N",
            "ATOM      3  CA  HIS A  17      11.000  10.000  10.000  1.00 20.00           C",
            "ATOM      4  ND1 HIS A  17      11.500  10.000  10.000  1.00 20.00           N",
            "HETATM    5  P   IHP B   1      12.000  10.000  10.000  1.00 20.00           P",
        ]
    )
    parsed = parse_design_request(
        request(
            reference_format="pdb",
            reference_content=pdb,
            functional_site_residues="A:K24, A:H17",
            functional_site_atoms={"A:24": "NZ", "A:17": "ND1"},
            ligand_code="IHP",
        )
    )
    assert list(parsed.functional_site_atoms) == ["A:24", "A:17"]


def test_unknown_functional_atom_is_rejected():
    pdb = "\n".join(
        [
            "ATOM      1  ND1 HIS A  17      10.000  10.000  10.000  1.00 20.00           N",
            "HETATM    2  P   IHP B   1      12.000  10.000  10.000  1.00 20.00           P",
        ]
    )
    with pytest.raises(ValueError, match="unknown atom"):
        parse_design_request(
            request(
                reference_format="pdb",
                reference_content=pdb,
                functional_site_residues="A:H17",
                functional_site_atoms={"A:17": "NZ"},
                ligand_code="IHP",
            )
        )
