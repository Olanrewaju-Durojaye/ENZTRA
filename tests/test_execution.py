import json
from types import SimpleNamespace

from enztra.config import EnztraConfig
from enztra.design import create_design_job, parse_design_request
from enztra.execution import (
    export_existing_mapping,
    prepare_rfdiffusion2_execution,
    run_rfdiffusion2,
)


def _prepared_job(tmp_path, design_length="150-180"):
    pdb = "\n".join(
        [
            "ATOM      1  CA  HIS A  17      10.000  20.000  30.000  1.00 20.00           C",
            "ATOM      2  ND1 HIS A  17      11.000  20.000  30.000  1.00 20.00           N",
            "HETATM    3  P   LIG B   1      12.000  22.000  32.000  1.00 20.00           P",
        ]
    )
    request = parse_design_request(
        {
            "job_id": "gpu-test", "reference_id": "active",
            "reference_format": "pdb", "reference_content": pdb,
            "substrate_name": "phytate", "substrate_smiles": "CCO",
            "ligand_code": "IHP", "structure_ligand_selector": "B:1",
            "functional_site_residues": "A:H17",
            "functional_site_atoms": {"A:17": "CA,ND1"},
            "atom_preservation_mode": "all_heavy", "protein_chain": "A",
            "design_length": design_length, "backbone_count": 3,
            "sequences_per_backbone": 4,
        }
    )
    jobs = tmp_path / "jobs"
    create_design_job(jobs, request)
    rfroot = tmp_path / "RFdiffusion2"
    (rfroot / "rf_diffusion/exec").mkdir(parents=True)
    (rfroot / "rf_diffusion/benchmark").mkdir(parents=True)
    (rfroot / "rf_diffusion/exec/bakerlab_rf_diffusion_aa.sif").touch()
    (rfroot / "rf_diffusion/benchmark/pipeline.py").touch()
    (rfroot / "rf_diffusion/run_inference.py").touch()
    config = EnztraConfig("conda", "apptainer", rfroot, tmp_path, tmp_path,
                          tmp_path, tmp_path, tmp_path)
    return jobs / "gpu-test", config


def test_execution_plan_uses_validated_workload(tmp_path):
    job, config = _prepared_job(tmp_path)
    plan = prepare_rfdiffusion2_execution(job, config)
    assert plan["backbone_count"] == 3
    assert "inference.num_designs=3" in plan["command"]
    assert "inference.ligand=LIG" in plan["command"]
    assert 'contigmap.contigs=["150-180,A17-17"]' in plan["command"]
    assert "contigmap.length=151-181" in plan["command"]
    assert 'contigmap.contig_atoms={A17:"CA,ND1"}' in plan["command"]
    assert "MKL_THREADING_LAYER=GNU" in plan["command"]


def test_execution_plan_preserves_exact_length_for_rfdiffusion2(tmp_path):
    job, config = _prepared_job(tmp_path, design_length="410")
    plan = prepare_rfdiffusion2_execution(job, config)
    assert 'contigmap.contigs=["410-410,A17-17"]' in plan["command"]
    assert 'contigmap.length="411"' in plan["command"]


def test_successful_execution_updates_status(tmp_path):
    job, config = _prepared_job(tmp_path)
    calls = []

    def runner(command, **kwargs):
        calls.append((command, kwargs))
        return SimpleNamespace(returncode=0)

    assert run_rfdiffusion2(job, config, runner) == 0
    status = json.loads((job / "status.json").read_text())
    assert status["stage"] == "rfdiffusion2_complete"
    assert status["rfdiffusion2_validated"] is True
    assert len(calls) == 2
    assert calls[0][1]["check"] is False


def test_missing_protein_guideposts_fail_the_stage(tmp_path):
    job, config = _prepared_job(tmp_path)
    returncodes = iter((0, 2))

    def runner(command, **kwargs):
        return SimpleNamespace(returncode=next(returncodes))

    assert run_rfdiffusion2(job, config, runner) == 2
    status = json.loads((job / "status.json").read_text())
    assert status["stage"] == "rfdiffusion2_failed"
    assert status["rfdiffusion2_validated"] is False
    assert status["rfdiffusion2_validation_returncode"] == 2


def test_existing_trb_mapping_can_be_exported_without_diffusion(tmp_path):
    job, config = _prepared_job(tmp_path)

    def runner(command, **kwargs):
        mapping = job / "outputs/functional_site_mapping.csv"
        mapping.parent.mkdir(parents=True, exist_ok=True)
        mapping.write_text(
            "backbone,reference_residue,scaffold_residue,fixed_atoms\n"
            'backbone_0,A:H17,A:H8,"CA,ND1"\n'
        )
        return SimpleNamespace(returncode=0)

    mapping = export_existing_mapping(job, config, runner)
    assert mapping.name == "functional_site_mapping.csv"
    status = json.loads((job / "status.json").read_text())
    assert status["functional_site_mapping"] == str(mapping)
