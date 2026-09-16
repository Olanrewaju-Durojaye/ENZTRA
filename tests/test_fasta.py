from pathlib import Path

import pytest

from enztra.fasta import read_fasta


def test_reads_multiline_fasta(tmp_path: Path) -> None:
    path = tmp_path / "records.fasta"
    path.write_text(">ref description\nACDE\nFGHI\n>design\nKLMN\n")
    records = read_fasta(path)
    assert [record.record_id for record in records] == ["ref", "design"]
    assert records[0].sequence == "ACDEFGHI"


def test_rejects_ambiguous_residue(tmp_path: Path) -> None:
    path = tmp_path / "bad.fasta"
    path.write_text(">bad\nACDX\n")
    with pytest.raises(ValueError, match="unsupported residues"):
        read_fasta(path)


def test_rejects_duplicate_ids(tmp_path: Path) -> None:
    path = tmp_path / "duplicates.fasta"
    path.write_text(">same first\nACDE\n>same second\nFGHI\n")
    with pytest.raises(ValueError, match="Duplicate"):
        read_fasta(path)

