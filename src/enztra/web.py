"""Local HTTP service for ENZTRA result exploration."""

from __future__ import annotations

import csv
import json
import mimetypes
from importlib.resources import files
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import unquote, urlparse

def _read_csv(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        rows: list[dict[str, Any]] = []
        for row in csv.DictReader(handle):
            converted: dict[str, Any] = {}
            for key, value in row.items():
                if value in (None, ""):
                    converted[key] = None
                elif value == "True":
                    converted[key] = True
                elif value == "False":
                    converted[key] = False
                else:
                    try:
                        converted[key] = float(value)
                    except ValueError:
                        converted[key] = value
            rows.append(converted)
        return rows


def load_job(job_dir: Path) -> dict[str, Any]:
    """Load one completed ENZTRA job into a JSON-ready document."""

    summary_path = job_dir / "summary.json"
    if not summary_path.is_file():
        raise FileNotFoundError(f"No summary.json in {job_dir}")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    boltz_summary_path = job_dir / "boltz2/summary.json"
    boltz_summary = (
        json.loads(boltz_summary_path.read_text(encoding="utf-8"))
        if boltz_summary_path.is_file()
        else None
    )
    return {
        "job_id": job_dir.name,
        "summary": summary,
        "selection": _read_csv(job_dir / "selection.csv"),
        "kinetics": _read_csv(job_dir / "kinetics.csv"),
        "boltz2": {
            "summary": boltz_summary,
            "models": _read_csv(job_dir / "boltz2/confidence.csv"),
            "best_models": _read_csv(
                job_dir / "boltz2/results/reports/best_models.csv"
            ),
        } if boltz_summary is not None else None,
    }


def list_jobs(jobs_root: Path) -> list[dict[str, Any]]:
    """List completed jobs below a root, newest first."""

    if not jobs_root.is_dir():
        return []
    jobs = []
    for summary_path in jobs_root.rglob("summary.json"):
        job_dir = summary_path.parent
        if not (job_dir / "selection.csv").is_file():
            continue
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        jobs.append(
            {
                "job_id": job_dir.relative_to(jobs_root).as_posix(),
                "design_count": summary.get("design_count", 0),
                "survivor_count": summary.get("survivor_count", 0),
                "substrate_id": summary.get("reference", {}).get("substrate_id", ""),
                "modified": summary_path.stat().st_mtime,
            }
        )
    return sorted(jobs, key=lambda job: job["modified"], reverse=True)


def make_handler(jobs_root: Path, static_dir: Path | None = None):
    """Build an HTTP handler bound to one local jobs directory."""

    root = jobs_root.resolve()
    assets = static_dir or Path(str(files("enztra").joinpath("static")))

    class EnztraHandler(BaseHTTPRequestHandler):
        def _send(self, status: int, body: bytes, content_type: str) -> None:
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)

        def _json(self, value: Any, status: int = 200) -> None:
            self._send(status, json.dumps(value).encode(), "application/json; charset=utf-8")

        def do_GET(self) -> None:  # noqa: N802 - required by BaseHTTPRequestHandler
            path = unquote(urlparse(self.path).path)
            if path == "/api/health":
                self._json({"status": "ready", "jobs_root": str(root)})
                return
            if path == "/api/jobs":
                self._json(list_jobs(root))
                return
            if path.startswith("/api/jobs/"):
                candidate = (root / path.removeprefix("/api/jobs/")).resolve()
                if candidate != root and root not in candidate.parents:
                    self._json({"detail": "Job not found"}, 404)
                    return
                try:
                    self._json(load_job(candidate))
                except (FileNotFoundError, json.JSONDecodeError):
                    self._json({"detail": "Job not found"}, 404)
                return
            relative = "index.html" if path == "/" else path.removeprefix("/static/")
            target = (assets / relative).resolve()
            if assets.resolve() not in target.parents or not target.is_file():
                self._json({"detail": "Not found"}, 404)
                return
            content_type = mimetypes.guess_type(target.name)[0] or "application/octet-stream"
            self._send(200, target.read_bytes(), content_type)

        def log_message(self, format: str, *args: Any) -> None:
            return

    return EnztraHandler


def serve(jobs_root: Path, host: str = "127.0.0.1", port: int = 8000) -> int:
    """Serve the dashboard until interrupted by the user."""

    server = ThreadingHTTPServer((host, port), make_handler(jobs_root))
    print(f"ENZTRA dashboard: http://{host}:{port}")
    print(f"Reading completed jobs from: {jobs_root.resolve()}")
    print("Press Ctrl+C to stop.")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nENZTRA dashboard stopped.")
    finally:
        server.server_close()
    return 0
