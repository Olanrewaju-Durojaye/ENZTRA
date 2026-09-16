"""Publication-friendly, dependency-free reporting for Boltz-2 results."""

from __future__ import annotations

import csv
import html
from pathlib import Path

from .plotting import display_design_label


BEST_MODEL_FIELDS = (
    "structural_rank",
    "design_id",
    "model_index",
    "confidence_score",
    "ligand_iptm",
    "complex_plddt",
    "iptm",
    "ptm",
    "protein_iptm",
    "complex_iplddt",
    "complex_pde",
    "complex_ipde",
    "kinetic_survivor_rank",
    "kcat_s",
    "km_mm",
    "catalytic_efficiency_s-1_mM-1",
    "structure_file",
    "confidence_file",
)


def write_best_models_csv(path: Path, rows: list[dict[str, object]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=BEST_MODEL_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return path


def write_boltz2_confidence_plot(
    path: Path,
    best_models: list[dict[str, object]],
    maximum_designs: int = 50,
    display_labels: dict[str, str] | None = None,
    rank_field: str = "structural_rank",
) -> Path:
    """Plot the best model for up to 50 structurally ranked survivors."""

    displayed = best_models[:maximum_designs]
    row_height = 30
    width = 1200
    top, bottom, left, right = 125, 90, 310, 70
    height = max(520, top + bottom + row_height * max(1, len(displayed)))
    plot_width = width - left - right
    colors = {
        "confidence_score": "#13795b",
        "ligand_iptm": "#7a4fa3",
        "complex_plddt": "#397b80",
    }
    labels = {
        "confidence_score": "Confidence score",
        "ligand_iptm": "Ligand ipTM",
        "complex_plddt": "Complex pLDDT",
    }
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">',
        '<title id="title">ENZTRA Boltz-2 best-model confidence</title>',
        '<desc id="desc">Best structural sample per strict kinetic survivor, ranked by Boltz-2 confidence score.</desc>',
        f'<rect width="{width}" height="{height}" fill="#f4f2ea"/>',
        f'<rect x="30" y="25" width="{width-60}" height="{height-50}" rx="18" fill="#fffdf7" stroke="#d9ddd5"/>',
        '<text x="65" y="65" font-family="Georgia,serif" font-size="28" font-weight="700" fill="#10231d">Boltz-2 best-model confidence</text>',
        '<text x="65" y="90" font-family="sans-serif" font-size="13" fill="#64736d">One automatically selected model per strict kinetic survivor; higher scores indicate greater predicted structural confidence.</text>',
    ]
    legend_x = 650
    for field in colors:
        parts.append(f'<circle cx="{legend_x}" cy="64" r="6" fill="{colors[field]}"/>')
        parts.append(f'<text x="{legend_x+11}" y="68" font-family="sans-serif" font-size="12" fill="#64736d">{labels[field]}</text>')
        legend_x += 170

    axis_bottom = height - bottom
    for tick in range(6):
        value = tick / 5
        x = left + value * plot_width
        parts.append(f'<line x1="{x:.2f}" y1="{top-12}" x2="{x:.2f}" y2="{axis_bottom}" stroke="#e4e7e1"/>')
        parts.append(f'<text x="{x:.2f}" y="{axis_bottom+27}" text-anchor="middle" font-family="sans-serif" font-size="12" fill="#64736d">{value:.1f}</text>')
    parts.append(f'<text x="{left+plot_width/2:.2f}" y="{height-40}" text-anchor="middle" font-family="sans-serif" font-size="14" font-weight="700" fill="#10231d">Boltz-2 confidence metric (0–1)</text>')

    offsets = {"confidence_score": -7, "ligand_iptm": 0, "complex_plddt": 7}
    for index, row in enumerate(displayed):
        cy = top + index * row_height
        design_id = str(row["design_id"])
        model_index = int(row["model_index"])
        rank = int(row[rank_field])
        full_title = html.escape(
            f"Rank {rank}: {design_id}, model {model_index}; "
            f"confidence={float(row['confidence_score']):.4f}, "
            f"ligand ipTM={float(row['ligand_iptm']):.4f}, "
            f"complex pLDDT={float(row['complex_plddt']):.4f}"
        )
        override = (display_labels or {}).get(design_id)
        visible = html.escape(override or display_design_label(design_id, 34))
        prefix = "" if override else f"#{rank} "
        parts.append(f'<text x="{left-16}" y="{cy+4}" text-anchor="end" font-family="sans-serif" font-size="12" fill="#10231d">{prefix}{visible}</text>')
        for field, color in colors.items():
            value = float(row[field])
            cx = left + value * plot_width
            point_y = cy + offsets[field]
            parts.append(f'<circle cx="{cx:.2f}" cy="{point_y:.2f}" r="5.5" fill="{color}" stroke="#fff" stroke-width="1.5"><title>{full_title}</title></circle>')

    if len(best_models) > maximum_designs:
        omitted = len(best_models) - maximum_designs
        parts.append(f'<text x="65" y="{height-40}" font-family="sans-serif" font-size="12" fill="#64736d">Showing the top {maximum_designs} designs; {omitted} additional designs remain in best_models.csv.</text>')
    if not displayed:
        parts.append(f'<text x="{width/2}" y="{height/2}" text-anchor="middle" font-family="sans-serif" font-size="16" fill="#64736d">No completed Boltz-2 models were found.</text>')
    parts.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")
    return path
