"""Source-neutral FASTA normalization for ENZTRA kinetic jobs."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

from .fasta import FastaRecord, read_fasta


def safe_identifier(header: str) -> str:
    preferred = header.split(",", 1)[0].split()[0]
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", preferred).strip("_.-")
    return value or "design"


def _metadata(header: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for field in header.split(",")[1:]:
        if "=" in field:
            key, value = field.split("=", 1)
            key = re.sub(r"[^A-Za-z0-9_.-]+", "_", key.strip()).strip("_")
            if key:
                values[key] = value.strip()
    return values


def normalize_external_records(
    source: Path,
    deduplicate: bool = True,
    reserved_ids: set[str] | None = None,
) -> tuple[list[FastaRecord], list[dict[str, str]]]:
    """Normalize external headers and optionally remove identical sequences."""

    source_records = read_fasta(source)
    used = set(reserved_ids or set())
    sequence_owner: dict[str, str] = {}
    imported: list[FastaRecord] = []
    provenance: list[dict[str, str]] = []
    for record in source_records:
        base = safe_identifier(record.description)
        design_id = base
        suffix = 2
        while design_id in used:
            design_id = f"{base}_{suffix}"
            suffix += 1
        used.add(design_id)
        duplicate_of = sequence_owner.get(record.sequence, "")
        included = not (deduplicate and duplicate_of)
        if not duplicate_of:
            sequence_owner[record.sequence] = design_id
        if included:
            imported.append(FastaRecord(design_id, record.description, record.sequence))
        provenance.append({
            "design_id": design_id,
            "original_header": record.description,
            "included": str(included),
            "duplicate_of": duplicate_of,
            "source_metadata": json.dumps(_metadata(record.description), sort_keys=True),
        })
    return imported, provenance


def write_imported_kinetic_inputs(
    job_dir: Path,
    reference: FastaRecord,
    designs: list[FastaRecord],
    provenance: list[dict[str, str]],
    source: Path,
    deduplicated: bool,
) -> tuple[Path, Path]:
    """Persist normalized FASTAs and complete source provenance."""

    inputs = job_dir / "imported_inputs"
    inputs.mkdir(parents=True, exist_ok=False)
    reference_path = inputs / "reference.fasta"
    designs_path = inputs / "designs.fasta"
    reference_path.write_text(
        f">{reference.record_id}\n{reference.sequence}\n", encoding="utf-8"
    )
    designs_path.write_text(
        "".join(f">{record.record_id}\n{record.sequence}\n" for record in designs),
        encoding="utf-8",
    )
    with (inputs / "source_metadata.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=(
            "design_id", "original_header", "included", "duplicate_of",
            "source_metadata",
        ))
        writer.writeheader()
        writer.writerows(provenance)
    (inputs / "import_summary.json").write_text(json.dumps({
        "source_file": str(source.expanduser().resolve()),
        "source_record_count": len(provenance),
        "imported_design_count": len(designs),
        "duplicate_count": sum(bool(row["duplicate_of"]) for row in provenance),
        "duplicates_removed": deduplicated,
        "reference_id": reference.record_id,
    }, indent=2) + "\n", encoding="utf-8")
    return reference_path, designs_path
