import json

from enztra.interactive import prepare_design_interactively, run_menu


def test_guided_preparation_saves_confirmed_project(tmp_path, monkeypatch, capsys):
    fasta = tmp_path / "reference.fasta"
    fasta.write_text(">WT_A\nACDEFGHIKLMNPQRSTVWY\n")
    answers = iter(
        [
            str(fasta),
            "general-enzyme-01",
            "",
            "",
            "test substrate",
            "LIG",
            "CCO",
            "A:C2, A:L10",
            "18-22",
            "2",
            "3",
            "test project",
            "y",
            "y",
        ]
    )
    monkeypatch.setattr("builtins.input", lambda prompt="": next(answers))
    assert prepare_design_interactively(tmp_path / "jobs") == 0
    assert (tmp_path / "jobs/general-enzyme-01/design_request.json").is_file()
    assert "Total planned designs:  6" in capsys.readouterr().out


def test_guided_preparation_writes_nothing_when_cancelled(tmp_path, monkeypatch):
    fasta = tmp_path / "reference.fasta"
    fasta.write_text(">WT_A\nACDEFGHIKLMNPQRSTVWY\n")
    answers = iter(
        [str(fasta), "cancelled-job", "", "", "substrate", "LIG", "CCO", "A:C2", "20", "1", "1", "", "n"]
    )
    monkeypatch.setattr("builtins.input", lambda prompt="": next(answers))
    assert prepare_design_interactively(tmp_path / "jobs") == 0
    assert not (tmp_path / "jobs").exists()


def test_pdb_quick_mode_preserves_all_heavy_atoms(tmp_path, monkeypatch, capsys):
    pdb = tmp_path / "reference.pdb"
    pdb.write_text(
        "\n".join(
            [
                "ATOM      1  N   HIS A  17      10.000  10.000  10.000  1.00 20.00           N",
                "ATOM      2  ND1 HIS A  17      11.000  10.000  10.000  1.00 20.00           N",
                "ATOM      3  H   HIS A  17      11.000  11.000  10.000  1.00 20.00           H",
                "ATOM      4  CA  ARG A  20      12.000  10.000  10.000  1.00 20.00           C",
                "ATOM      5  NH1 ARG A  20      13.000  10.000  10.000  1.00 20.00           N",
                "HETATM    6  P   IHP B   1      14.000  10.000  10.000  1.00 20.00           P",
            ]
        )
    )
    answers = iter(
        [
            str(pdb), "quick-atoms", "", "", "substrate", "IHP", "CCO",
            "A:H17, A:R20", "", "", "2", "1", "1", "", "y", "y",
        ]
    )
    monkeypatch.setattr("builtins.input", lambda prompt="": next(answers))
    assert prepare_design_interactively(tmp_path / "jobs") == 0
    specification = json.loads(
        (tmp_path / "jobs/quick-atoms/design_request.json").read_text()
    )
    assert specification["atom_preservation_mode"] == "all_heavy"
    assert specification["functional_site_atoms"] == {
        "A:17": ["N", "ND1"],
        "A:20": ["CA", "NH1"],
    }
    output = capsys.readouterr().out
    assert "Atom-preservation mode: All available heavy atoms" in output
    assert "Residues covered:       2/2" in output


def test_main_menu_reprompts_then_exits(tmp_path, monkeypatch, capsys):
    answers = iter(["wrong", "15"])
    monkeypatch.setattr("builtins.input", lambda prompt="": next(answers))
    assert run_menu(tmp_path / "jobs", tmp_path / "config.json") == 0
    output = capsys.readouterr().out
    assert "Invalid choice" in output
    assert "Goodbye" in output
