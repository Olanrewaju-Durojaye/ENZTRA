import csv
import json
import threading
from http.server import ThreadingHTTPServer
from urllib.request import urlopen

from enztra.web import list_jobs, load_job, make_handler


def _job(tmp_path):
    job = tmp_path / "kinetics_run"
    job.mkdir()
    (job / "summary.json").write_text(
        json.dumps({"reference": {"enzyme_id": "ref", "substrate_id": "S"}, "design_count": 1, "survivor_count": 1})
    )
    with (job / "selection.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["design_id", "kcat_s", "passes_strict_gate"])
        writer.writeheader()
        writer.writerow({"design_id": "d1", "kcat_s": 12.5, "passes_strict_gate": True})
    return job


def test_load_job_converts_numbers_and_booleans(tmp_path):
    job = _job(tmp_path)
    data = load_job(job)
    assert data["selection"][0] == {"design_id": "d1", "kcat_s": 12.5, "passes_strict_gate": True}


def test_list_jobs_ignores_incomplete_directories(tmp_path):
    job = _job(tmp_path)
    (tmp_path / "unfinished").mkdir()
    (job / "boltz2").mkdir()
    (job / "boltz2/summary.json").write_text(json.dumps({"survivor_count": 1}))
    jobs = list_jobs(tmp_path)
    assert [job["job_id"] for job in jobs] == ["kinetics_run"]


def test_load_job_includes_boltz2_best_models(tmp_path):
    job = _job(tmp_path)
    reports = job / "boltz2/results/reports"
    reports.mkdir(parents=True)
    (job / "boltz2/summary.json").write_text(json.dumps({
        "survivor_count": 1, "diffusion_samples": 1, "msa_mode": "single"
    }))
    with (reports / "best_models.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=(
            "structural_rank", "design_id", "confidence_score"
        ))
        writer.writeheader()
        writer.writerow({
            "structural_rank": 1, "design_id": "d1", "confidence_score": 0.8
        })

    data = load_job(job)
    assert data["boltz2"]["summary"]["msa_mode"] == "single"
    assert data["boltz2"]["best_models"][0]["confidence_score"] == 0.8


def test_http_api_serves_job_data(tmp_path):
    _job(tmp_path)
    server = ThreadingHTTPServer(("127.0.0.1", 0), make_handler(tmp_path))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        port = server.server_address[1]
        with urlopen(f"http://127.0.0.1:{port}/api/jobs") as response:
            jobs = json.load(response)
        assert jobs[0]["job_id"] == "kinetics_run"
    finally:
        server.shutdown()
        server.server_close()
