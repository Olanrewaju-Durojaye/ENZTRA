import csv
import json

from enztra.design import pdb_chain_sequence
from enztra.fasta import FastaRecord
from enztra.kinetic_import import (
    normalize_external_records,
    write_imported_kinetic_inputs,
)


def test_metadata_headers_are_cleaned_and_preserved(tmp_path):
    source = tmp_path / "designs.txt"
    source.write_text(
        ">run_Main1-bb-False, id=0, overall_confidence=0.3890\nACDE\n"
    )
    records, provenance = normalize_external_records(source)
    assert records[0].record_id == "run_Main1-bb-False"
    assert provenance[0]["original_header"].endswith("overall_confidence=0.3890")
    metadata = json.loads(provenance[0]["source_metadata"])
    assert metadata["overall_confidence"] == "0.3890"


def test_identical_sequences_can_be_deduplicated(tmp_path):
    source = tmp_path / "designs.fasta"
    source.write_text(">first\nACDE\n>second\nACDE\n>third\nFGHI\n")
    records, provenance = normalize_external_records(source, deduplicate=True)
    assert [record.record_id for record in records] == ["first", "third"]
    assert provenance[1]["included"] == "False"
    assert provenance[1]["duplicate_of"] == "first"


def test_identical_sequences_can_be_retained(tmp_path):
    source = tmp_path / "designs.fasta"
    source.write_text(">first\nACDE\n>second\nACDE\n")
    records, provenance = normalize_external_records(source, deduplicate=False)
    assert len(records) == 2
    assert all(row["included"] == "True" for row in provenance)


def test_normalized_identifier_collisions_receive_suffixes(tmp_path):
    source = tmp_path / "designs.fasta"
    source.write_text(">same/a\nACDE\n>same_a\nFGHI\n")
    records, _ = normalize_external_records(source)
    assert [record.record_id for record in records] == ["same_a", "same_a_2"]


def test_pdb_reference_chain_is_converted_to_sequence():
    pdb = "\n".join([
        "ATOM      1  CA  HIS A  17      10.000  10.000  10.000  1.00 20.00           C",
        "ATOM      2  CB  HIS A  17      11.000  10.000  10.000  1.00 20.00           C",
        "ATOM      3  CA  ARG A  20      12.000  10.000  10.000  1.00 20.00           C",
    ])
    assert pdb_chain_sequence(pdb, "A") == "HR"


def test_import_writes_normalized_fastas_and_provenance(tmp_path):
    source = tmp_path / "source.fasta"
    source.write_text(">design one\nACDE\n")
    designs, provenance = normalize_external_records(source)
    reference = FastaRecord("reference", "reference", "FGHI")
    reference_path, designs_path = write_imported_kinetic_inputs(
        tmp_path / "job", reference, designs, provenance, source, True
    )
    assert reference_path.read_text() == ">reference\nFGHI\n"
    assert designs_path.read_text() == ">design\nACDE\n"
    with (tmp_path / "job/imported_inputs/source_metadata.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["original_header"] == "design one"
