import csv
import json

from enztra.publication import generate_publication_exports


SELECTION_FIELDS = [
    "design_id", "sequence", "kcat_s", "km_mm",
    "catalytic_efficiency_s-1_mM-1", "passes_kcat", "passes_km",
    "passes_strict_gate", "decision", "kcat_fold_change", "km_fold_change",
    "efficiency_fold_change", "survivor_rank",
]


def _write_csv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _kinetic_job(tmp_path, survivor=True):
    job = tmp_path / "job"
    job.mkdir()
    passed = "True" if survivor else "False"
    _write_csv(job / "selection.csv", SELECTION_FIELDS, [{
        "design_id": "design_1", "sequence": "ACDEFG", "kcat_s": "20",
        "km_mm": "0.5" if survivor else "2", "catalytic_efficiency_s-1_mM-1": "40",
        "passes_kcat": "True", "passes_km": passed,
        "passes_strict_gate": passed, "decision": "ACCEPTED" if survivor else "REJECTED",
        "kcat_fold_change": "2", "km_fold_change": "0.5" if survivor else "2",
        "efficiency_fold_change": "4", "survivor_rank": "1" if survivor else "",
    }])
    _write_csv(job / "kinetics.csv", [
        "design_id", "role", "sequence", "kcat_s", "km_mm", "km_sd_total",
        "km_sd_aleatoric", "km_sd_epistemic",
    ], [{
        "design_id": "design_1", "role": "design", "sequence": "ACDEFG",
        "kcat_s": "20", "km_mm": "0.5", "km_sd_total": "0.2",
        "km_sd_aleatoric": "0.1", "km_sd_epistemic": "0.1",
    }])
    (job / "summary.json").write_text(json.dumps({
        "reference": {"enzyme_id": "WT", "kcat_s": 10, "km_mm": 1},
        "design_count": 1, "survivor_count": int(survivor),
        "boltz2_queue": ["design_1"] if survivor else [],
    }))
    (job / "manifest.json").write_text(json.dumps({
        "reference_id": "WT", "substrate_name": "S", "substrate_smiles": "CCO",
    }))
    (job / "results").mkdir()
    (job / "results/kinetic_selection.svg").write_text("<svg/>")
    return job


def test_zero_survivor_publication_export_is_complete(tmp_path):
    job = _kinetic_job(tmp_path, survivor=False)
    result = generate_publication_exports(job)
    assert result["strict_survivor_count"] == 0
    assert result["structurally_ranked_count"] == 0
    assert (job / "publication/tables/master_ranking.csv").is_file()
    assert (job / "publication/tables/master_ranking.json").is_file()
    assert (job / "publication/figures/candidate_outcomes.svg").is_file()
    assert (job / "publication/reproducibility_manifest.json").is_file()
    rows = list(csv.DictReader((job / "publication/tables/master_ranking.csv").open()))
    assert rows[0]["decision"] == "REJECTED"
    assert rows[0]["final_priority_rank"] == ""
    assert rows[0]["qc_status"] == "PASS"


def test_survivor_combines_kinetic_structure_provenance_and_qc(tmp_path):
    job = _kinetic_job(tmp_path, survivor=True)
    structure = job / "model.cif"
    structure.write_text("data_model")
    reports = job / "boltz2/results/reports"
    _write_csv(reports / "best_models.csv", [
        "structural_rank", "design_id", "model_index", "confidence_score",
        "ligand_iptm", "complex_plddt", "structure_file",
    ], [{
        "structural_rank": "1", "design_id": "design_1", "model_index": "0",
        "confidence_score": "0.8", "ligand_iptm": "0.9",
        "complex_plddt": "0.85", "structure_file": str(structure),
    }])
    (job / "boltz2/status.json").parent.mkdir(parents=True, exist_ok=True)
    (job / "boltz2/status.json").write_text(json.dumps({
        "stage": "boltz2_complete", "returncode": 0,
        "best_models_csv": str(reports / "best_models.csv"),
    }))
    metadata = job / "imported_inputs/source_metadata.csv"
    _write_csv(metadata, ["design_id", "original_header"], [{
        "design_id": "design_1", "original_header": "long original header",
    }])
    result = generate_publication_exports(job)
    assert result["structurally_ranked_count"] == 1
    rows = list(csv.DictReader((job / "publication/tables/master_ranking.csv").open()))
    assert rows[0]["original_header"] == "long original header"
    assert rows[0]["final_priority_rank"] == "1"
    assert rows[0]["publication_label"] == "Priority 1"
    assert float(rows[0]["combined_priority_score"]) > 0
    assert rows[0]["qc_status"] == "PASS"
    kinetic_plot = (job / "publication/figures/kinetic_selection.svg").read_text()
    boltz_plot = (job / "publication/figures/boltz2_confidence.svg").read_text()
    assert "Priority 1" in kinetic_plot
    assert "Priority 1" in boltz_plot


def test_missing_expected_structure_is_qc_failure(tmp_path):
    job = _kinetic_job(tmp_path, survivor=True)
    (job / "boltz2").mkdir()
    (job / "boltz2/status.json").write_text(json.dumps({
        "stage": "boltz2_complete", "returncode": 0,
    }))
    generate_publication_exports(job)
    rows = list(csv.DictReader((job / "publication/tables/master_ranking.csv").open()))
    assert rows[0]["qc_status"] == "FAIL"
    assert "missing_boltz2_best_model" in rows[0]["qc_errors"]


def test_low_structural_confidence_is_warning_not_gate_reversal(tmp_path):
    job = _kinetic_job(tmp_path, survivor=True)
    (job / "model.cif").write_text("data_model")
    reports = job / "boltz2/results/reports"
    _write_csv(reports / "best_models.csv", [
        "structural_rank", "design_id", "model_index", "confidence_score",
        "ligand_iptm", "complex_plddt", "structure_file",
    ], [{
        "structural_rank": "1", "design_id": "design_1", "model_index": "0",
        "confidence_score": "0.6", "ligand_iptm": "0.65",
        "complex_plddt": "0.68", "structure_file": "model.cif",
    }])
    generate_publication_exports(job)
    rows = list(csv.DictReader((job / "publication/tables/master_ranking.csv").open()))
    assert rows[0]["passes_strict_gate"] == "True"
    assert rows[0]["qc_status"] == "WARN"
    assert "low_confidence_score" in rows[0]["qc_warnings"]


def test_km_uncertainty_threshold_is_opt_in(tmp_path):
    job = _kinetic_job(tmp_path, survivor=False)
    generate_publication_exports(job)
    rows = list(csv.DictReader((job / "publication/tables/master_ranking.csv").open()))
    assert "km_uncertainty_above_review_threshold" not in rows[0]["qc_warnings"]

    generate_publication_exports(job, km_uncertainty_review_threshold=0.1)
    rows = list(csv.DictReader((job / "publication/tables/master_ranking.csv").open()))
    assert "km_uncertainty_above_review_threshold" in rows[0]["qc_warnings"]
