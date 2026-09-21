# ENZTRA

<img width="1230" height="1278" alt="ENZTRA_Logo" src="https://github.com/user-attachments/assets/811250b1-574c-4afe-913b-07e438e28534" />

**ENZTRA** (ENzyme redesign and TRiage Automation) is a reference-guided
workflow for enzyme redesign, kinetic triage, and structural analysis.

ENZTRA is being developed to connect four independent scientific tools without
redistributing their source code or model weights:

1. RFdiffusion2 and LigandMPNN generate candidate backbones and sequences.
2. DLKcat predicts `kcat` for the reference and every design.
3. CatPred predicts `Km` and uncertainty for the same enzyme–substrate pairs.
4. A strict kinetic gate retains only designs with both higher `kcat` and lower
   `Km` than the reference.
5. Boltz-2 models only the survivors and supports structural reliability and
   catalytic-geometry analysis. Affinity prediction is not required.

> ENZTRA predictions are computational hypotheses. They do not establish
> catalytic activity and must be validated experimentally.

## Step-by-step tutorial

A complete beginner-friendly installation and reproduction guide is available in
[TUTORIAL.md](TUTORIAL.md). It includes a preliminary smoke test and the exact
input conditions used for the `release-500` and `release-500b` demonstrations.

## Strict selection rule

For reference values `kcat_ref` and `Km_ref`, a design survives only when:

```text
kcat_design > kcat_ref AND Km_design < Km_ref
```

Equality fails. Catalytic efficiency (`kcat/Km`) ranks strict survivors but
never rescues a design that fails either threshold.

On the interactive plot, `Km` is the x-axis and `kcat` is the y-axis.
The reference point defines vertical and horizontal boundaries; survivors
occupy the upper-left region.

## Version 1.0

The package currently implements the model-independent scientific core and
the first local-tool integration layer:

- validated reference and candidate records;
- strict two-condition selection;
- fold-change and catalytic-efficiency calculations;
- survivor-only ranking and Boltz-2 queue generation;
- a command-line interface;
- automated tests, including boundary and high-efficiency failure cases;
- an illustrative 1DKP–IHP phytase example;
- strict multi-record FASTA validation;
- equivalent DLKcat and CatPred batch-input generation;
- isolated Conda execution adapters;
- row-count and identity checks before predictions are combined;
- a read-only installation doctor;
- a prepare-only mode for validating jobs without launching model inference;
- a local, read-only browser dashboard for completed kinetic jobs;
- an interactive Km-versus-kcat decision map, strict-gate shading, tooltips,
  linear/log axes, searchable results, and the Boltz-2 survivor queue;
- a guided terminal menu for validated design-job preparation;
- FASTA/PDB, substrate, ligand, functional-site-residue, and library-size validation;
- identity-aware functional-site constraints and inactive-mutation warnings;
- reproducible local `design_request.json` and reference-input storage;
- guarded RFdiffusion2 backbone generation;
- ligand-aware LigandMPNN sequence generation from completed backbones;
- consolidated design FASTA output for subsequent kinetic prediction;
- survivor-only Boltz-2 protein–ligand structure prediction;
- validated confidence collection and best-model summaries;
- source-neutral kinetic intake for external or ENZTRA-generated sequences.
- resumable end-to-end execution with output validation and failure recovery.
- consolidated rankings, QC flags, and publication/reproducibility exports.
- controlled-pilot preflight and planned-versus-observed integration reports.

Version `1.0.1` is the current stable ENZTRA release. It retains the native
v1.0.0 qualification and corrects RFdiffusion2 length accounting when protein
guideposts are present. The original native qualification
run completed the exact 50-backbone × 10-sequence plan: all 500 sequences
reached kinetic evaluation, 14 strict survivors reached Boltz-2 structural
evaluation, all publication exports completed, and no integration check
failed. A concise verification record is provided under `examples/release_500/`.

## Full 500-sequence release qualification

Prepare a new project with exactly 50 RFdiffusion2 backbones and 10 LigandMPNN
sequences per backbone. Before starting expensive inference, run the preflight:

```bash
enztra demonstrate \
  --config enztra.config.json \
  --job-dir jobs/release-500 \
  --calibration-job jobs/pilot-50 \
  --prepare-only
```

The successful 50-design pilot calibrates the disk estimate. Review
`jobs/release-500/release/preflight.json`, then remove `--prepare-only` to run
or resume the workflow. Single-sequence Boltz-2 mode and one structural sample
per strict survivor are the reproducible defaults:

```bash
enztra demonstrate \
  --config enztra.config.json \
  --job-dir jobs/release-500 \
  --calibration-job jobs/pilot-50 \
  --msa-mode single \
  --diffusion-samples 1
```

You can also choose **Run or resume the full 500-sequence release
demonstration** in the terminal menu. A stopped run resumes from its first
incomplete stage. Final evidence is written to
`release/qualification_report.json` and `.md`. A qualifying run records
`github_release_ready: true`.

## Developer quick start

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[test]"
pytest
```

Run the illustrative kinetic gate:

```bash
enztra select \
  --reference examples/phytase/reference.json \
  --designs examples/phytase/designs.csv \
  --output outputs/phytase/results.csv \
  --summary outputs/phytase/summary.json
```

Expected result: two of four designs survive, ranked by predicted catalytic
efficiency and placed in the Boltz-2 queue.

## Local tool configuration

Copy the example configuration and adjust paths only if your installations are
elsewhere:

```bash
cp enztra.config.example.json enztra.config.json
enztra doctor --config enztra.config.json
```

Prepare a kinetic job without running either model:

```bash
enztra kinetics \
  --config enztra.config.json \
  --reference-fasta examples/phytase/reference.fasta \
  --designs-fasta examples/phytase/designs.fasta \
  --substrate-name IHP \
  --substrate-smiles 'O=P(O)(O)OC1C(OP(=O)(O)O)C(OP(=O)(O)O)C(OP(=O)(O)O)C(OP(=O)(O)O)C1OP(=O)(O)O' \
  --job-dir jobs/phytase_prepare_test \
  --prepare-only
```

Remove `--prepare-only` to run DLKcat and GPU CatPred, combine their outputs,
apply the strict kinetic gate, and create `summary.json`, `selection.csv`, the
Boltz-2 queue, and a permanent `results/kinetic_selection.svg` decision plot in
the job directory. The SVG opens directly in a browser or figure editor; the
interactive dashboard remains optional.

To keep large screens readable, permanent plots label only the reference and
strict survivors. Rejected candidates remain as red points. Long generated
identifiers are shortened only on the visible plot label; their complete IDs
remain in the SVG hover text, CSV/JSON outputs, and downstream Boltz-2 queue.

## External and ENZTRA-generated sequence intake

Choose **Run kinetic triage on protein sequences** from the terminal menu to
use either a completed ENZTRA LigandMPNN FASTA or an external `.fasta`, `.fa`,
`.faa`, or FASTA-formatted `.txt` file. The reference may be a one-record FASTA
or a PDB; for PDB input, ENZTRA asks for the active protein chain and extracts
its sequence.

Metadata-rich headers from tools such as BOLTRA and LigandMPNN are normalized
to safe unique design IDs. Original headers and parsed key-value metadata are
retained in `imported_inputs/source_metadata.csv`. Identical sequences are
reported before execution; the default removes redundant copies while recording
which retained design they duplicate, and an alternative keeps every record.
The normalized inputs feed the same DLKcat/CatPred and strict-gate workflow as
native ENZTRA sequences.

## Guided terminal workflow

Run ENZTRA without arguments:

```bash
enztra
```

The numbered menu prepares a new enzyme-redesign project, runs RFdiffusion2,
runs LigandMPNN, exports or verifies functional-site mappings, inspects projects,
opens the results dashboard, or checks local installations. During
project preparation, ENZTRA asks one question at a time, re-prompts invalid
answers, prints a final review, and writes files only after confirmation.
The designed-protein length prompt accepts either an exact residue count such
as `200` or an inclusive range such as `180-220`; pressing Enter retains the
selected reference chain's length.

For PDB references, ENZTRA offers two atom-preservation modes. The default
one-step mode preserves all available heavy atoms for every selected
functional-site residue. Advanced mode lists the atoms in each residue and
lets the user choose a custom subset. Successful preflight writes the complete
mapping to `rfdiffusion2_spec.json`; it does not start GPU inference
automatically.

## Resumable complete workflow

Choose **Run or resume the complete workflow** after preparing a PDB-based
project, or use the non-interactive command:

```bash
enztra workflow \
  --config enztra.config.json \
  --job-dir jobs/<project> \
  --msa-mode single \
  --diffusion-samples 1
```

ENZTRA runs RFdiffusion2, LigandMPNN, kinetic prediction and strict selection,
then Boltz-2 for strict survivors. Progress is written after every transition
to `jobs/<project>/workflow_state.json`. Each completed stage is skipped only
after its required outputs validate on disk. If execution fails or is stopped,
run the same menu action or command again: the failed or interrupted stage is
retried while validated earlier stages are preserved.

If no sequence beats both reference thresholds, the workflow completes
normally and records Boltz-2 as skipped. This is a valid result rather than an
execution error. In v0.12.0, it then continues to the publication-export stage.

## Publication and quality-control exports

Choose **Generate publication and quality-control exports** from the terminal
menu, or run:

```bash
enztra publication-export --job-dir jobs/<project>
```

The command creates `jobs/<project>/publication/` containing:

```text
tables/master_ranking.csv
tables/master_ranking.json
tables/quality_control.csv
figures/candidate_outcomes.svg
figures/kinetic_selection.svg
figures/boltz2_confidence.svg        # when Boltz-2 results exist
publication_summary.json
reproducibility_manifest.json
publication_report.md
```

The strict kinetic decision remains authoritative. Structural rank and the
combined priority rank apply only to strict survivors with valid Boltz-2
results. The combined score is the dimensionless geometric mean of kcat fold
improvement, inverse Km fold, confidence score, ligand ipTM, and complex
pLDDT; it is a transparent prioritization aid, not evidence of activity.

QC failures identify invalid or missing required data. CatPred uncertainty is
always exported, but an uncertainty review threshold is opt-in because the
initial `0.5` value has not been established as a universal scientific cutoff.
Structural metrics below the documented 0.70 review threshold remain warnings
without reversing the strict kinetic decision. SHA-256
values in the reproducibility manifest simply identify whether files changed;
users do not need to calculate them manually.

To request an explicit CatPred uncertainty review threshold:

```bash
enztra publication-export \
  --job-dir jobs/<project> \
  --km-uncertainty-review-threshold 0.5
```

## Controlled native pilot

Prepare a project with 10 RFdiffusion2 backbones and 5 LigandMPNN sequences per
backbone, then choose **Run or resume a controlled native pilot** from the menu.
The equivalent command is:

```bash
enztra pilot \
  --config enztra.config.json \
  --job-dir jobs/pilot-50 \
  --msa-mode single \
  --diffusion-samples 1
```

Use `--prepare-only` to run preflight without starting a model. The pilot writes
`pilot/preflight.json`, `pilot/integration_report.json`, and
`pilot/integration_report.md`. The report verifies every planned-versus-observed
handoff and is the primary diagnostic artifact for finalizing v0.13.0.

### Validated 50-design pilot

The qualifying `pilot-50` run used 10 RFdiffusion2 backbones and 5 LigandMPNN
sequences per backbone. Every planned-versus-observed check passed:

| Measure | Planned | Observed |
| --- | ---: | ---: |
| Completed backbones | 10 | 10 |
| Protein sequences | 50 | 50 |
| Kinetic candidates | 50 | 50 |
| Strict kinetic survivors | — | 2 |
| Structurally ranked survivors | 2 | 2 |
| Publication candidates | 50 | 50 |

Observed stage times on an NVIDIA RTX 4000 Ada system were approximately 77.6
minutes for RFdiffusion2, 22.5 seconds for LigandMPNN, 15.3 seconds for kinetic
triage, 39.2 seconds for Boltz-2, and 0.6 seconds for publication export. These
times describe one local run and are not hardware-independent benchmarks.
The machine-readable and Markdown qualification records are included under
`examples/pilot_50/`.

## Guarded RFdiffusion2 backbone generation

Choose **Run a prepared RFdiffusion2 project** from the terminal menu after a
PDB project passes preflight. ENZTRA creates a direct `run_inference.py` plan,
shows the exact Apptainer command and workload, and starts no GPU process until
the user confirms. Cancelling leaves a reusable command plan in
`rfdiffusion2_execution.json`. This stage generates backbones only.
ENZTRA isolates RFdiffusion2 from host Intel-MKL threading settings and passes
each atom-level Hydra override as one process argument. After generation, each
TRB must contain the expected number of mapped protein guideposts; otherwise the
stage is marked failed and the output is blocked from LigandMPNN.
The validated original-to-scaffold residue correspondence and fixed atom names
are exported automatically to `outputs/functional_site_mapping.csv`.

## Guarded LigandMPNN sequence generation

After backbone generation succeeds, choose **Run LigandMPNN on completed
backbones**. ENZTRA detects complete PDB/TRB backbone pairs, uses RFdiffusion2's
ligand-aware MPNN v2 pipeline, and retains the functional-site positions fixed
during sequence generation. It shows the exact command and workload before
starting and uses the sequence count selected when the project was prepared.

Generated sequences are consolidated into:

```text
jobs/<project>/outputs/ligandmpnn/designs.fasta
```

Before accepting this FASTA, ENZTRA requires LigandMPNN's fixed-position file
to match the TRB-derived scaffold mapping exactly. It also validates sequence
count separately for every planned backbone and rejects stale output belonging
to other backbones. Amino-acid symbols are checked before the project is marked
complete, making the FASTA the handoff to kinetic prediction.

## Guarded Boltz-2 structural validation

After kinetic prediction creates `summary.json`, `selection.csv`, and
`manifest.json`, choose **Run Boltz-2 on strict kinetic survivors**. ENZTRA
refuses to model rejected candidates and cross-checks every queued identifier,
strict-gate decision, sequence, and substrate SMILES before creating inputs.

The default single-sequence mode writes `msa: empty`. It is the most reliable
local option but may reduce structural accuracy. MSA-server mode omits the
`msa` field and adds `--use_msa_server`; it can improve evolutionary context
but requires network access and depends on the external server. The execution
review shows the exact command before the GPU starts.

ENZTRA does not request Boltz-2's optional affinity property. It collects
`confidence_score`, pTM, ipTM, ligand ipTM, complex pLDDT, interface pLDDT,
PDE, and interface PDE into `boltz2/confidence.csv`, with best models summarized
in `boltz2/summary.json`. Each sample is ranked within its design; the best
sample from every survivor is then ranked globally by confidence score, ligand
ipTM, and complex pLDDT. Publication-friendly reports are written to:

```text
jobs/<kinetic-job>/boltz2/results/reports/best_models.csv
jobs/<kinetic-job>/boltz2/results/reports/boltz2_confidence.svg
```

The best-model table also carries each design's kinetic rank, kcat, Km, and
catalytic efficiency. Complete structure and confidence-file paths remain
traceable. These are structural-confidence estimates, not proof of catalytic
activity or experimental binding affinity.

For a completed Boltz-2 run created by an earlier ENZTRA release, choose
**Generate or refresh Boltz-2 result reports** from the menu. This reads the
existing confidence and structure files and does not rerun prediction or use
the GPU. The equivalent direct command is:

```bash
enztra boltz2-report --job-dir jobs/completed_kinetics_job
```

The same stage can be prepared or launched directly:

```bash
enztra boltz2 \
  --config enztra.config.json \
  --job-dir jobs/completed_kinetics_job \
  --msa-mode single \
  --diffusion-samples 1 \
  --prepare-only
```

Remove `--prepare-only` to start Boltz-2.

ENZTRA also prepares the RFdiffusion2 ORI token. If the reference contains an
ORI token, it can be reused. Otherwise, the default places ORI at the selected
protein chain's C-alpha centroid; advanced users may enter custom X, Y, and Z
coordinates. The generated `rfdiffusion2_input.pdb` retains the selected
protein and intended substrate but removes unrelated HETATM records such as
waters, Hg ions, and crystallization additives. The untouched reference PDB is
kept alongside it for provenance.

### Reference structures and functional-site residues

Experimental PDB structures may contain inactive mutations, engineered
substitutions, unresolved residues, or missing loops. Confirm that the supplied
sequence represents the intended active enzyme. When possible, specify each
catalytic or substrate-interacting constraint with its expected identity:

```text
A:R16, A:H17, A:R20, A:D88, A:R92, A:H303, A:D304
```

ENZTRA compares these identities with the reference and reports mismatches such
as an expected catalytic histidine being represented by alanine. Position-only
notation such as `A:17` remains accepted, but it cannot detect an identity
mismatch.

Use the full descriptive compound name for **substrate name**, the PDB Chemical
Component Dictionary identifier for **ligand CCD code**, and a stereochemically
defined structure for **substrate SMILES**. For example:

```text
Substrate name: myo-inositol hexakisphosphate
Ligand CCD code: IHP
```

Boltz-generated structures commonly label every small molecule as `LIG`. ENZTRA
automatically maps a single `LIG` residue to the substrate code entered above
and uses the stereochemical SMILES as the authoritative chemical identity. If
several `LIG` residues are present, the terminal workflow asks which
chain/residue identifies the intended substrate. Do not enter `LIG` as the
substrate CCD code: retain the original scientific code, such as `IHP`.

For direct access to the questionnaire, use:

```bash
enztra prepare-design --jobs-root jobs
```

## Local results dashboard

Launch the dependency-free interface from the ENZTRA repository:

```bash
enztra serve --jobs-root jobs
```

Open <http://127.0.0.1:8000>. The results-only dashboard automatically discovers
completed kinetic jobs, displays kinetic selection and any available Boltz-2
best-model ranking, and never modifies them. Stop the server with `Ctrl+C` in
its terminal.

## Planned deployment model

The public repository will provide one installer and one browser interface,
while keeping incompatible model stacks isolated. Large third-party weights
and containers will be downloaded from their official sources and will not be
committed to GitHub.

## License

ENZTRA's original orchestration code is released under the MIT License.
RFdiffusion2, DLKcat, CatPred, Boltz-2, their weights, and their dependencies
retain their respective licenses and citation requirements.
