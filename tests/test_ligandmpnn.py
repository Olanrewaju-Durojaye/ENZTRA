import json
from types import SimpleNamespace

from enztra.ligandmpnn import prepare_ligandmpnn_execution, run_ligandmpnn

from test_execution import _prepared_job


def _completed_backbone_job(tmp_path):
    job, config = _prepared_job(tmp_path)
    output = job / "outputs/rfdiffusion2"
    output.mkdir(parents=True)
    (output / "backbone_0.pdb").write_text("ATOM\nHETATM\n")
    (output / "backbone_0.trb").write_bytes(b"test")
    mapping = job / "outputs/functional_site_mapping.csv"
    mapping.write_text(
        "backbone,reference_residue,scaffold_residue,fixed_atoms\n"
        'backbone_0,A:H17,A:H8,"CA,ND1"\n'
    )
    status_path = job / "status.json"
    status = json.loads(status_path.read_text())
    status.update({
        "stage": "rfdiffusion2_complete",
        "rfdiffusion2_returncode": 0,
        "rfdiffusion2_validated": True,
    })
    status_path.write_text(json.dumps(status))
    return job, config


def test_ligandmpnn_plan_uses_stored_sequence_count(tmp_path):
    job, config = _completed_backbone_job(tmp_path)
    plan = prepare_ligandmpnn_execution(job, config)
    assert plan["backbone_count"] == 1
    assert plan["sequences_per_backbone"] == 4
    assert plan["total_sequences"] == 4
    assert plan["backbone_ids"] == ["backbone_0"]
    assert "start_step=mpnn" in plan["command"]
    assert "stop_step=thread_mpnn" in plan["command"]
    assert "use_ligand=True" in plan["command"]
    assert "mpnn.num_seq_per_target=4" in plan["command"]
    assert "mpnn.slurm.in_proc=True" in plan["command"]


def test_successful_ligandmpnn_collects_clean_fasta(tmp_path):
    job, config = _completed_backbone_job(tmp_path)

    def runner(command, **kwargs):
        ligmpnn = job / "outputs/rfdiffusion2/ligmpnn"
        seqs = ligmpnn / "seqs"
        seqs.mkdir(parents=True)
        ligmpnn.joinpath("pdbs_position_fixed_0.jsonl").write_text(
            json.dumps({str(job / "outputs/rfdiffusion2/backbone_0.pdb"): ["A8"]})
            + "\n"
        )
        seqs.joinpath("backbone_0.fa").write_text(
            ">native\nACDE\n"
            ">sample=1\nACDF\n"
            ">sample=2\nACDG\n"
            ">sample=3\nACDH\n"
            ">sample=4\nACDI\n"
        )
        return SimpleNamespace(returncode=0)

    assert run_ligandmpnn(job, config, runner) == 0
    status = json.loads((job / "status.json").read_text())
    assert status["stage"] == "ligandmpnn_complete"
    assert status["ligandmpnn_sequence_count"] == 4
    fasta = (job / "outputs/ligandmpnn/designs.fasta").read_text()
    assert fasta.count(">") == 4
    assert "ACDE" not in fasta


def test_ligandmpnn_rejects_wrong_fixed_positions(tmp_path):
    job, config = _completed_backbone_job(tmp_path)

    def runner(command, **kwargs):
        ligmpnn = job / "outputs/rfdiffusion2/ligmpnn"
        seqs = ligmpnn / "seqs"
        seqs.mkdir(parents=True)
        ligmpnn.joinpath("pdbs_position_fixed_0.jsonl").write_text(
            json.dumps({str(job / "outputs/rfdiffusion2/backbone_0.pdb"): []}) + "\n"
        )
        seqs.joinpath("backbone_0.fa").write_text(">native\nACDE\n>sample=1\nACDF\n")
        return SimpleNamespace(returncode=0)

    assert run_ligandmpnn(job, config, runner) == 2
    status = json.loads((job / "status.json").read_text())
    assert status["stage"] == "ligandmpnn_failed"
    assert status["ligandmpnn_fixed_positions_validated"] is False
    assert "fixed-position mismatch" in status["ligandmpnn_validation_error"]


def test_ligandmpnn_ignores_stale_fasta_from_another_backbone(tmp_path):
    job, config = _completed_backbone_job(tmp_path)

    def runner(command, **kwargs):
        ligmpnn = job / "outputs/rfdiffusion2/ligmpnn"
        seqs = ligmpnn / "seqs"
        seqs.mkdir(parents=True)
        ligmpnn.joinpath("pdbs_position_fixed_0.jsonl").write_text(
            json.dumps({str(job / "outputs/rfdiffusion2/backbone_0.pdb"): ["A8"]})
            + "\n"
        )
        seqs.joinpath("backbone_0.fa").write_text(
            ">native\nACDE\n>sample=1\nACDF\n>sample=2\nACDG\n"
            ">sample=3\nACDH\n>sample=4\nACDI\n"
        )
        seqs.joinpath("old_backbone.fa").write_text(
            ">native\nACDE\n>old\nAAAA\n"
        )
        return SimpleNamespace(returncode=0)

    assert run_ligandmpnn(job, config, runner) == 0
    fasta = (job / "outputs/ligandmpnn/designs.fasta").read_text()
    assert fasta.count(">") == 4
    assert "old_backbone" not in fasta
