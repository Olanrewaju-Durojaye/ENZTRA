"""Small, strict FASTA reader used at ENZTRA's workflow boundary."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True, slots=True)
class FastaRecord:
    record_id: str
    description: str
    sequence: str


def read_fasta(path: Path) -> list[FastaRecord]:
    """Parse FASTA records and reject ambiguous or malformed protein input."""

    records: list[FastaRecord] = []
    header: str | None = None
    chunks: list[str] = []

    def finish_record() -> None:
        if header is None:
            return
        sequence = "".join(chunks).replace(" ", "").upper()
        if not sequence:
            raise ValueError(f"FASTA record {header!r} has no sequence")
        invalid = sorted(set(sequence) - AMINO_ACIDS)
        if invalid:
            raise ValueError(
                f"FASTA record {header!r} contains unsupported residues: {invalid}"
            )
        record_id = header.split()[0]
        if not record_id:
            raise ValueError("FASTA record identifier must not be empty")
        records.append(FastaRecord(record_id, header, sequence))

    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish_record()
            header = line[1:].strip()
            chunks = []
            if not header:
                raise ValueError(f"Empty FASTA header at line {line_number}")
        else:
            if header is None:
                raise ValueError(f"Sequence found before FASTA header at line {line_number}")
            chunks.append(line)
    finish_record()

    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    identifiers = [record.record_id for record in records]
    duplicates = sorted({item for item in identifiers if identifiers.count(item) > 1})
    if duplicates:
        raise ValueError(f"Duplicate FASTA identifiers: {duplicates}")
    return records

