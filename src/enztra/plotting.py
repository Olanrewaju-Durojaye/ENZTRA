"""Dependency-free, permanent kinetic-selection plots."""

from __future__ import annotations

import html
import re
from pathlib import Path
from typing import Iterable

from .models import ReferenceThresholds, SelectionResult


def _scale(value: float, low: float, high: float, start: float, end: float) -> float:
    if high == low:
        return (start + end) / 2
    return start + ((value - low) / (high - low)) * (end - start)


def display_design_label(label: str, limit: int = 32) -> str:
    """Make a compact label while preserving backbone and sequence identity."""

    if len(label) <= limit:
        return label
    generated = re.search(
        r"_cond\d+_(\d+).*_(seq(?:uence)?_?\d+)$", label, re.IGNORECASE
    )
    if generated:
        backbone, sequence = generated.groups()
        return f"bb{int(backbone):03d} · {sequence.lower()}"
    left = (limit - 1) // 2
    right = limit - left - 1
    return f"{label[:left]}…{label[-right:]}"


def write_kinetic_selection_plot(
    path: Path,
    reference: ReferenceThresholds,
    results: Iterable[SelectionResult],
    display_labels: dict[str, str] | None = None,
) -> Path:
    """Write a standalone SVG decision map and return its path."""

    results = list(results)
    points = [(reference.km_mm, reference.kcat_s)] + [
        (result.candidate.km_mm, result.candidate.kcat_s) for result in results
    ]
    xs = [point[0] for point in points]
    ys = [point[1] for point in points]
    x_span = max(xs) - min(xs)
    y_span = max(ys) - min(ys)
    x_pad = x_span * 0.18 if x_span else max(abs(xs[0]) * 0.08, 0.01)
    y_pad = y_span * 0.18 if y_span else max(abs(ys[0]) * 0.08, 1.0)
    x_min, x_max = max(0.0, min(xs) - x_pad), max(xs) + x_pad
    y_min, y_max = max(0.0, min(ys) - y_pad), max(ys) + y_pad

    width, height = 1100, 720
    left, right, top, bottom = 110, 55, 100, 105
    plot_right, plot_bottom = width - right, height - bottom
    x = lambda value: _scale(value, x_min, x_max, left, plot_right)
    y = lambda value: _scale(value, y_min, y_max, plot_bottom, top)
    ref_x, ref_y = x(reference.km_mm), y(reference.kcat_s)

    elements = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="720" '
        'viewBox="0 0 1100 720" role="img" aria-labelledby="title desc">',
        '<title id="title">ENZTRA kinetic selection decision map</title>',
        '<desc id="desc">Kcat versus Km. Strict survivors have higher kcat and lower Km than the reference.</desc>',
        '<rect width="1100" height="720" fill="#f4f2ea"/>',
        '<rect x="35" y="30" width="1030" height="655" rx="18" fill="#fffdf7" stroke="#d9ddd5"/>',
        '<text x="70" y="67" font-family="Georgia,serif" font-size="28" font-weight="700" fill="#10231d">Kinetic selection decision map</text>',
        '<text x="70" y="90" font-family="sans-serif" font-size="13" fill="#64736d">Upper-left shaded region: kcat above reference and Km below reference (strict inequalities).</text>',
        f'<rect x="{left:.2f}" y="{top:.2f}" width="{max(0, ref_x-left):.2f}" height="{max(0, ref_y-top):.2f}" fill="#dff3df" opacity="0.85"/>',
        f'<line x1="{left}" y1="{plot_bottom}" x2="{plot_right}" y2="{plot_bottom}" stroke="#89958f"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{plot_bottom}" stroke="#89958f"/>',
        f'<line x1="{ref_x:.2f}" y1="{top}" x2="{ref_x:.2f}" y2="{plot_bottom}" stroke="#829a91" stroke-dasharray="6 6"/>',
        f'<line x1="{left}" y1="{ref_y:.2f}" x2="{plot_right}" y2="{ref_y:.2f}" stroke="#829a91" stroke-dasharray="6 6"/>',
    ]

    for tick in range(5):
        fraction = tick / 4
        x_value = x_min + fraction * (x_max - x_min)
        y_value = y_min + fraction * (y_max - y_min)
        x_pos = left + fraction * (plot_right - left)
        y_pos = plot_bottom - fraction * (plot_bottom - top)
        elements.extend([
            f'<text x="{x_pos:.2f}" y="{plot_bottom + 27}" text-anchor="middle" font-family="sans-serif" font-size="12" fill="#64736d">{x_value:.4g}</text>',
            f'<text x="{left - 13}" y="{y_pos + 4:.2f}" text-anchor="end" font-family="sans-serif" font-size="12" fill="#64736d">{y_value:.4g}</text>',
        ])

    elements.extend([
        f'<text x="{(left + plot_right) / 2:.2f}" y="{height - 50}" text-anchor="middle" font-family="sans-serif" font-size="15" font-weight="700" fill="#10231d">Km (mM) — lower is better</text>',
        f'<text x="40" y="{(top + plot_bottom) / 2:.2f}" transform="rotate(-90 40 {(top + plot_bottom) / 2:.2f})" text-anchor="middle" font-family="sans-serif" font-size="15" font-weight="700" fill="#10231d">kcat (s⁻¹) — higher is better</text>',
    ])

    def add_point(
        label: str,
        km: float,
        kcat: float,
        color: str,
        radius: int,
        show_label: bool,
    ) -> None:
        safe_label = html.escape(label)
        cx, cy = x(km), y(kcat)
        elements.append(
            f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{radius}" fill="{color}" stroke="#ffffff" stroke-width="3"><title>{safe_label}: kcat={kcat:.6g} s⁻¹, Km={km:.6g} mM</title></circle>'
        )
        if show_label:
            visible = (display_labels or {}).get(label, display_design_label(label))
            visible_label = html.escape(visible)
            elements.append(
                f'<text x="{cx + 11:.2f}" y="{cy - 10:.2f}" font-family="sans-serif" font-size="12" fill="#10231d">{visible_label}</text>'
            )

    add_point(
        reference.enzyme_id,
        reference.km_mm,
        reference.kcat_s,
        "#397b80",
        8,
        True,
    )
    for result in results:
        add_point(
            result.candidate.design_id,
            result.candidate.km_mm,
            result.candidate.kcat_s,
            "#13795b" if result.passes_strict_gate else "#c55a4b",
            7,
            result.passes_strict_gate,
        )

    legend_y = 660
    for legend_x, color, label in (
        (700, "#397b80", "Reference"),
        (820, "#13795b", "Strict survivor"),
        (965, "#c55a4b", "Rejected"),
    ):
        elements.append(f'<circle cx="{legend_x}" cy="{legend_y}" r="6" fill="{color}"/>')
        elements.append(f'<text x="{legend_x + 11}" y="{legend_y + 4}" font-family="sans-serif" font-size="12" fill="#64736d">{label}</text>')

    elements.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(elements) + "\n", encoding="utf-8")
    return path
