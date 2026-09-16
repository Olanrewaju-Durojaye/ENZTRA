import csv
from pathlib import Path

import pytest

from enztra.adapters import combine_predictions, prepare_inputs
from enztra.fasta import FastaRecord


RECORDS = [
    FastaRecord("reference", "reference", "ACDE"),
    FastaRecord("design_1", "design_1", "FGHI"),
]


def test_prepare_inputs_uses_same_pairs_for_both_tools(tmp_path: Path) -> None:
    dlkcat, catpred = prepare_inputs(RECORDS, "substrate", "CCO", tmp_path)
    with dlkcat.open(newline="") as handle:
        dlkcat_rows = list(csv.DictReader(handle, delimiter="\t"))
    with catpred.open(newline="") as handle:
        catpred_rows = list(csv.DictReader(handle))
    assert [row["Protein Sequence"] for row in dlkcat_rows] == ["ACDE", "FGHI"]
    assert [row["pdbpath"] for row in catpred_rows] == [
        "reference.pdb",
        "design_1.pdb",
    ]
    assert {row["Substrate SMILES"] for row in dlkcat_rows} == {"CCO"}
    assert {row["SMILES"] for row in catpred_rows} == {"CCO"}


def test_combine_predictions_validates_identity(tmp_path: Path) -> None:
    dlkcat = tmp_path / "dlkcat.tsv"
    dlkcat.write_text(
        "Substrate Name\tSubstrate SMILES\tProtein Sequence\tKcat value (1/s)\n"
        "s\tCCO\tACDE\t10\n"
        "s\tCCO\tFGHI\t20\n"
    )
    catpred = tmp_path / "catpred.csv"
    catpred.write_text(
        "pdbpath,Prediction_(mM),SD_total,SD_aleatoric,SD_epistemic\n"
        "reference.pdb,2,0.3,0.2,0.1\n"
        "design_1.pdb,1,0.4,0.3,0.2\n"
    )
    predictions = combine_predictions(RECORDS, dlkcat, catpred)
    assert predictions[0].kcat_s == 10
    assert predictions[1].km_mm == 1
    assert predictions[1].km_sd_total == pytest.approx(0.4)


def test_combine_predictions_rejects_row_loss(tmp_path: Path) -> None:
    dlkcat = tmp_path / "dlkcat.tsv"
    dlkcat.write_text(
        "Protein Sequence\tKcat value (1/s)\nACDE\t10\n"
    )
    catpred = tmp_path / "catpred.csv"
    catpred.write_text("pdbpath,Prediction_(mM)\nreference.pdb,2\n")
    with pytest.raises(ValueError, match="row-count mismatch"):
        combine_predictions(RECORDS, dlkcat, catpred)

