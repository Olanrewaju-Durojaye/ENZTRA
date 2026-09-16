import json

import pytest

from enztra.workflow import load_workflow_state, run_workflow


def _prepared_job(tmp_path):
    job = tmp_path / "jobs" / "resume-test"
    job.mkdir(parents=True)
    (job / "design_request.json").write_text(json.dumps({
        "job_id": "resume-test",
        "reference_id": "WT",
        "reference_format": "pdb",
        "reference_file": "inputs/reference.pdb",
        "protein_chain": "A",
        "substrate_name": "substrate",
        "substrate_smiles": "CCO",
    }))
    (job / "status.json").write_text(json.dumps({
        "job_id": "resume-test", "stage": "prepared", "execution_ready": True,
    }))
    config = tmp_path / "config.json"
    config.write_text(json.dumps({
        "rfdiffusion2_root": str(tmp_path / "rf"),
        "dlkcat_root": str(tmp_path / "dlkcat"),
        "catpred_root": str(tmp_path / "catpred"),
        "catpred_data_root": str(tmp_path / "catpred-data"),
    }))
    return job, config


def test_failure_stops_downstream_and_resume_skips_completed(tmp_path):
    job, config = _prepared_job(tmp_path)
    calls = []

    def first(stage):
        calls.append(stage)
        return 7 if stage == "ligandmpnn" else 0

    assert run_workflow(job, config, stage_runner=first) == 7
    state = load_workflow_state(job)
    assert calls == ["rfdiffusion2", "ligandmpnn"]
    assert state["stages"]["rfdiffusion2"]["status"] == "completed"
    assert state["stages"]["ligandmpnn"]["status"] == "failed"
    assert state["stages"]["kinetics"]["status"] == "pending"

    # Create the durable evidence required by resume reconciliation.
    (job / "status.json").write_text(json.dumps({
        "rfdiffusion2_returncode": 0, "rfdiffusion2_validated": True,
    }))
    output = job / "outputs/rfdiffusion2"
    output.mkdir(parents=True)
    (output / "design.pdb").write_text("MODEL\n")
    (output / "design.trb").write_bytes(b"trb")
    (job / "outputs/functional_site_mapping.csv").write_text("backbone\n")
    second_calls = []

    def second(stage):
        second_calls.append(stage)
        if stage == "kinetics":
            (job / "summary.json").write_text(json.dumps({"boltz2_queue": ["d1"]}))
        return 0

    assert run_workflow(job, config, stage_runner=second) == 0
    assert second_calls == ["ligandmpnn", "kinetics", "boltz2", "publication"]
    state = load_workflow_state(job)
    assert state["workflow_status"] == "completed"
    assert state["stages"]["ligandmpnn"]["attempts"] == 2


def test_interruption_is_saved_for_retry(tmp_path):
    job, config = _prepared_job(tmp_path)

    def interrupted(stage):
        if stage == "rfdiffusion2":
            raise KeyboardInterrupt
        return 0

    with pytest.raises(KeyboardInterrupt):
        run_workflow(job, config, stage_runner=interrupted)
    state = load_workflow_state(job)
    assert state["workflow_status"] == "interrupted"
    assert state["stages"]["rfdiffusion2"]["status"] == "interrupted"
    assert state["stages"]["rfdiffusion2"]["attempts"] == 1


def test_zero_survivors_skip_boltz2_successfully(tmp_path):
    job, config = _prepared_job(tmp_path)
    calls = []

    def runner(stage):
        calls.append(stage)
        if stage == "kinetics":
            (job / "summary.json").write_text(json.dumps({"boltz2_queue": []}))
        return 0

    assert run_workflow(job, config, stage_runner=runner) == 0
    assert calls == ["rfdiffusion2", "ligandmpnn", "kinetics", "publication"]
    state = load_workflow_state(job)
    assert state["workflow_status"] == "completed"
    assert state["stages"]["boltz2"]["status"] == "skipped"


def test_running_state_becomes_interrupted_on_load(tmp_path):
    job, _ = _prepared_job(tmp_path)
    state = load_workflow_state(job)
    state["stages"]["rfdiffusion2"]["status"] = "running"
    (job / "workflow_state.json").write_text(json.dumps(state))
    from enztra.workflow import reconcile_workflow_state

    reconciled = reconcile_workflow_state(job)
    assert reconciled["stages"]["rfdiffusion2"]["status"] == "interrupted"
