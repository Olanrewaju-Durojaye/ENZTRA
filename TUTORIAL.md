# ENZTRA step-by-step tutorial

This tutorial is intended for first-time ENZTRA users, manuscript reviewers, and researchers who want to reproduce the two phytase demonstrations described with ENZTRA v1.0.1.

ENZTRA (**Enzyme Redesign and Triage Automation**) connects separately installed tools for:

1. reference-enzyme preparation;
2. RFdiffusion2 backbone generation;
3. LigandMPNN sequence generation;
4. DLKcat and CatPred kinetic prediction;
5. strict kinetic selection;
6. Boltz-2 structure prediction for strict kinetic survivors; and
7. quality control, traceability, and publication-oriented export.

> **Important:** ENZTRA produces computational hypotheses. Predicted kinetic improvements, structural models, confidence scores, and rankings do not establish enzyme activity. Selected candidates require experimental validation.

---

## Contents

- [1. What this tutorial reproduces](#1-what-this-tutorial-reproduces)
- [2. Before installation](#2-before-installation)
- [3. Install ENZTRA](#3-install-enztra)
- [4. Install and configure external tools](#4-install-and-configure-external-tools)
- [5. Check the installation](#5-check-the-installation)
- [6. Prepare the phytase reference](#6-prepare-the-phytase-reference)
- [7. Recommended smoke test](#7-recommended-smoke-test)
- [8. Reproduce `release-500`](#8-reproduce-release-500)
- [9. Reproduce `release-500b`](#9-reproduce-release-500b)
- [10. Run and resume a complete workflow](#10-run-and-resume-a-complete-workflow)
- [11. Monitor a long run](#11-monitor-a-long-run)
- [12. Understand the outputs](#12-understand-the-outputs)
- [13. Verify a completed demonstration](#13-verify-a-completed-demonstration)
- [14. Common problems](#14-common-problems)
- [15. Reproducibility checklist](#15-reproducibility-checklist)

---

## 1. What this tutorial reproduces

The two demonstration projects use the same enzyme, substrate, functional-site constraints, number of backbones, and number of sequences per backbone. Their principal difference is the requested scaffold length.

| Setting | `release-500` | `release-500b` |
|---|---:|---:|
| Reference | active *E. coli* phytase, PDB 1DKP | same |
| Substrate | myo-inositol hexakisphosphate (IHP; phytate) | same |
| RFdiffusion2 backbones | 50 | 50 |
| LigandMPNN sequences per backbone | 10 | 10 |
| Total planned sequences | 500 | 500 |
| Requested design length | 150–180 residues | exactly 410 residues |
| Observed strict kinetic survivors in the reported run | 14 | 0 |
| Boltz-2 behavior in the reported run | modeled 14 survivors | skipped because no candidate passed the strict gate |

The ENZTRA strict kinetic gate is:

```text
predicted kcat > reference kcat
AND
predicted Km < reference Km
```

The gate uses strict inequalities. Passing quality control is not the same as passing this kinetic gate.

---

## 2. Before installation

### 2.1 Recommended platform

The demonstrated workflow was run on Linux with an NVIDIA GPU. A recent Ubuntu release is recommended.

You will need:

- Git;
- Miniconda, Anaconda, or Mambaforge;
- Python 3.11 or 3.12 for ENZTRA;
- an NVIDIA GPU with a working driver for the GPU stages;
- Apptainer for the demonstrated RFdiffusion2 installation;
- sufficient disk space for external repositories, model weights, containers, trajectories, structures, and results.

The exact GPU-memory and runtime requirements depend on sequence length, external-tool versions, prediction options, and hardware. Do not begin with 500 sequences. Complete the smoke test in [Section 7](#7-recommended-smoke-test) first.

### 2.2 Check basic software

```bash
git --version
conda --version
python --version
apptainer --version
nvidia-smi
```

If `nvidia-smi` fails, repair the NVIDIA driver before attempting RFdiffusion2 or Boltz-2 GPU inference.

### 2.3 Suggested directory layout

The following layout keeps ENZTRA and its external tools separate:

```text
~/Applications/
├── ENZTRA/
├── RFdiffusion2/
├── DLKcat/
├── catpred_pipeline/
└── boltz/                 # or a separately installed Boltz environment
```

The names are examples. ENZTRA should be configured with the actual absolute paths on your machine.

---

## 3. Install ENZTRA

### 3.1 Clone the repository

```bash
mkdir -p "$HOME/Applications"
cd "$HOME/Applications"
git clone https://github.com/Olanrewaju-Durojaye/ENZTRA.git
cd ENZTRA
```

If you downloaded a release ZIP instead, extract it and enter the extracted directory:

```bash
cd "$HOME/Downloads"
unzip ENZTRA-v1.0.1.zip -d "$HOME/Applications"
cd "$HOME/Applications/ENZTRA-v1.0.1"
```

Use either the cloned repository or the extracted release directory, not both for the same run.

### 3.2 Create a clean environment

```bash
conda create -n enztra python=3.12 -y
conda activate enztra
```

### 3.3 Install ENZTRA

From the ENZTRA repository root:

```bash
python -m pip install --upgrade pip
python -m pip install -e ".[test]"
```

The editable installation is convenient for a local clone. For a fixed release snapshot, a regular installation can be used instead:

```bash
python -m pip install ".[test]"
```

### 3.4 Run the test suite

```bash
pytest
```

For ENZTRA v1.0.1, the complete distributed test suite should pass before scientific runs are started. The release used for the manuscript passed 91 tests.

### 3.5 Confirm the installed version and command

```bash
python -c "import enztra; print(enztra.__version__)"
command -v enztra
enztra
```

The version command should report `1.0.1` when reproducing these demonstrations.

---

## 4. Install and configure external tools

ENZTRA orchestrates external software; it does not redistribute those programs, their containers, or their model weights. Install each dependency according to its own license and official instructions.

### 4.1 RFdiffusion2 and LigandMPNN

Follow the official RFdiffusion2 installation instructions:

- repository: <https://github.com/RosettaCommons/RFdiffusion2>
- documentation: <https://rosettacommons.github.io/RFdiffusion2/>

The demonstrated installation used Apptainer and the RFdiffusion2 image:

```text
rf_diffusion/exec/bakerlab_rf_diffusion_aa.sif
```

A typical setup is:

```bash
cd "$HOME/Applications"
git clone https://github.com/RosettaCommons/RFdiffusion2.git
cd RFdiffusion2
python setup.py
```

The RFdiffusion2 setup downloads large files. Follow the upstream recovery instructions if a download is interrupted.

LigandMPNN is invoked through the RFdiffusion2 pipeline used by the demonstrated workflow. Confirm that the RFdiffusion2 setup has the expected LigandMPNN capability before running ENZTRA.

### 4.2 DLKcat

Install DLKcat from its official repository and verify its example prediction independently:

<https://github.com/SysBioChalmers/DLKcat>

ENZTRA supplies protein sequences and the substrate representation to the configured DLKcat predictor and imports the resulting predicted `kcat` values.

### 4.3 CatPred

Install CatPred and its pretrained models by following:

<https://github.com/maranasgroup/CatPred>

The upstream installation uses a dedicated Conda environment. Verify the upstream batch example before connecting it to ENZTRA. ENZTRA uses CatPred for predicted `Km` and its reported uncertainty components.

### 4.4 Boltz-2

Install Boltz in a separate environment using the current upstream instructions:

<https://github.com/jwohlwend/boltz>

For a CUDA installation, the upstream project currently recommends:

```bash
pip install "boltz[cuda]" -U
```

The demonstrated ENZTRA workflow runs Boltz-2 only for candidates that pass the strict kinetic gate.

### 4.5 Configure ENZTRA paths

From the ENZTRA root, make a local configuration from the distributed example:

```bash
cd "$HOME/Applications/ENZTRA"
cp enztra.config.example.json enztra.config.json
```

Open `enztra.config.json` in a text editor and replace the example locations with the absolute locations of the external installations and executables on your machine. Do not commit machine-specific paths, access tokens, or credentials to GitHub.

```bash
nano enztra.config.json
```

If the example file in your release uses a different documented configuration filename, preserve that filename and edit the existing keys rather than inventing new keys.

> **Why absolute paths?** Long GPU jobs may be launched from different working directories. Absolute paths avoid accidentally resolving an executable, container, or checkpoint relative to the wrong directory.

---

## 5. Check the installation

Activate ENZTRA and use the installation check before preparing a project:

```bash
conda activate enztra
cd "$HOME/Applications/ENZTRA"
enztra
```

At the main menu, choose:

```text
14. Check local installations
```

Resolve every missing required dependency before a complete run. A command being present is not enough: run the dependency's own small example to confirm that its model weights, environment, CUDA access, and output permissions work.

Useful independent checks include:

```bash
apptainer exec --nv /absolute/path/to/bakerlab_rf_diffusion_aa.sif nvidia-smi
boltz --help
```

DLKcat and CatPred should be tested using the examples distributed by their maintainers.

---

## 6. Prepare the phytase reference

### 6.1 Reference file

The demonstrations use an active representation derived from the phytase–phytate structure PDB 1DKP.

Example local path:

```text
~/Documents/ENZTRA_Test/Test1/1DKP_Active.pdb
```

Replace `~` with your actual home directory if a program requires an absolute path. Confirm the file exists:

```bash
test -s "$HOME/Documents/ENZTRA_Test/Test1/1DKP_Active.pdb" \
  && echo "Reference found" \
  || echo "Reference missing"
```

The reference must represent the intended active enzyme. Experimental PDB files can contain engineered mutations, missing residues, alternative conformers, unresolved loops, nonstandard ligand labels, or crystallization additives. ENZTRA cannot determine biological intent on the user's behalf.

### 6.2 Shared substrate settings

Use the following values for both demonstrations:

```text
Substrate name: Myo-inositol hexakisphosphate
Ligand CCD code: IHP
```

Stereochemical SMILES:

```text
[C@@H]1([C@@H]([C@@H]([C@@H]([C@H]([C@@H]1OP(=O)(O)O)OP(=O)(O)O)OP(=O)(O)O)OP(=O)(O)O)OP(=O)(O)O)OP(=O)(O)O
```

Keep the stereochemistry symbols (`@` and `@@`) unchanged. For the archived demonstration, the structure's generic ligand label `LIG` at `B:1` was mapped to IHP, while this stereochemical SMILES remained the authoritative substrate identity.

### 6.3 Functional-site residues

Use this explicit list:

```text
A:R16, A:H17, A:R20, A:K24, A:D88, A:R92, A:S212, A:S215, A:M216, A:Q253, A:R267, A:H303, A:D304, A:T305
```

`A:S215` is written explicitly here. An early saved request displayed the shortened token `A:215`; the reference residue and retained `OG` atom identify this position as serine. New projects should use the unambiguous `A:S215` form.

### 6.4 Shared interactive choices

For both demonstrations:

- protein chain: `A`;
- atom-preservation mode: `1`, preserve all available heavy atoms;
- ORI placement mode: `1`, automatic protein C-alpha centroid;
- RFdiffusion2 backbones: `50`;
- LigandMPNN sequences per backbone: `10`;
- total planned designs: `500`.

The automatic ORI coordinates recorded for the demonstrated reference were approximately:

```text
-0.330, -1.790, -0.421
```

These coordinates should be calculated from the exact input structure. Do not force the recorded values onto a different PDB preparation.

---

## 7. Recommended smoke test

Before either 500-sequence demonstration, prepare a one-backbone, one-sequence project. This verifies the complete interface while limiting wasted GPU time.

Launch ENZTRA:

```bash
conda activate enztra
cd "$HOME/Applications/ENZTRA"
enztra
```

Choose:

```text
1. Prepare a new enzyme-redesign project
```

Enter the shared reference, substrate, chain, functional residues, atom-preservation mode, and ORI choices from Section 6. Use:

```text
Project name: exact-410-smoke
Designed protein length: 410
RFdiffusion2 backbones: 1
LigandMPNN sequences per backbone: 1
```

Confirm the active reference and save the project. Then choose:

```text
7. Run or resume the complete workflow
```

Select `exact-410-smoke` and review the displayed commands before approving execution.

The first 410-residue RFdiffusion2 design can appear quiet for several minutes. Progress should eventually show flow-matching steps. Use the monitoring commands in Section 11 rather than assuming the process has frozen.

Do not proceed to 500 sequences unless:

- RFdiffusion2 completes and its outputs validate;
- LigandMPNN produces the requested sequence;
- kinetic inputs and outputs are created;
- the selection decision is recorded;
- Boltz-2 either runs for a strict survivor or is explicitly skipped when none exists; and
- publication and qualification outputs are generated.

---

## 8. Reproduce `release-500`

### 8.1 Prepare the project

Start ENZTRA and select option 1:

```text
1. Prepare a new enzyme-redesign project
```

Use the following responses. Text in angle brackets is explanatory and should not be typed.

| Prompt | Entry |
|---|---|
| Reference PDB or FASTA path | `~/Documents/ENZTRA_Test/Test1/1DKP_Active.pdb` |
| Project name | `release-500` |
| Reference enzyme ID | `1DKP_Active` |
| Protein chain | `A` |
| Descriptive substrate name | `Myo-inositol hexakisphosphate` |
| Substrate CCD code | `IHP` |
| Substrate stereochemical SMILES | use the exact SMILES in Section 6.2 |
| Functional-site residues | use the explicit 14-residue list in Section 6.3 |
| Atom-preservation mode | `1` |
| ORI-placement mode | `1` |
| Designed protein length or range | `150-180` |
| Number of RFdiffusion2 backbones | `50` |
| LigandMPNN sequences per backbone | `10` |
| Optional notes | `Short-scaffold 150-180-residue demonstration` |
| Confirm intended active reference | `y` |
| Save prepared project | `y` |

Read the project review carefully. It should report 50 backbones, 10 sequences per backbone, and 500 planned designs.

### 8.2 Run the release demonstration

From the main menu, choose:

```text
9. Run or resume the full 500-sequence release demonstration
```

Select `release-500`, review the preflight summary, and confirm only when:

- the output directory is correct;
- installation status is `READY`;
- available disk space exceeds the reported requirement; and
- no other job is using the required GPU.

Option 9 runs or resumes the release workflow and skips stages only when their expected outputs validate. Alternatively, option 7 runs or resumes the same prepared project's general complete workflow.

### 8.3 Expected reported outcome

The manuscript demonstration produced:

```text
Planned backbones:        50
Completed backbones:      50
Generated sequences:      500
Kinetic candidates:       500
Strict kinetic survivors: 14
Structurally ranked:       14
```

The values above are a reproducibility reference, not a universal guarantee. Differences in software versions, model weights, seeds, hardware behavior, or reference preparation can affect results.

---

## 9. Reproduce `release-500b`

`release-500b` uses the same conditions as `release-500`, except that its requested protein length is exactly 410 residues, matching the reference length.

### 9.1 Prepare the project

Choose option 1 and use:

| Prompt | Entry |
|---|---|
| Reference PDB or FASTA path | `~/Documents/ENZTRA_Test/Test1/1DKP_Active.pdb` |
| Project name | `release-500b` |
| Reference enzyme ID | `1DKP_Active` |
| Protein chain | `A` |
| Descriptive substrate name | `Myo-inositol hexakisphosphate` |
| Substrate CCD code | `IHP` |
| Substrate stereochemical SMILES | use the exact SMILES in Section 6.2 |
| Functional-site residues | use the explicit 14-residue list in Section 6.3 |
| Atom-preservation mode | `1` |
| ORI-placement mode | `1` |
| Designed protein length | `410` |
| Number of RFdiffusion2 backbones | `50` |
| LigandMPNN sequences per backbone | `10` |
| Optional notes | `Full-length 410-residue comparison run; all other settings match release-500.` |
| Confirm intended active reference | `y` |
| Save prepared project | `y` |

Use ENZTRA v1.0.1 or later for this fixed-length guidepost case. Version 1.0.1 corrected RFdiffusion2 guidepost-length compatibility handling while retaining the requested 410-residue scaffold length.

### 9.2 Run the project

Choose:

```text
9. Run or resume the full 500-sequence release demonstration
```

Select `release-500b`, approve the validated preflight, and leave the ENZTRA terminal open while the job is active.

### 9.3 Expected reported outcome

The manuscript demonstration produced:

```text
Planned backbones:        50
Completed backbones:      50
Generated sequences:      500
Kinetic candidates:       500
Strict kinetic survivors: 0
Structurally ranked:       0
QC PASS:                  500
QC WARN:                  0
QC FAIL:                  0
```

All candidate predicted `kcat` values were below the reference prediction, so no candidate passed the strict kinetic gate and Boltz-2 was not run. This is a scientifically valid negative selection outcome, not a failed workflow.

---

## 10. Run and resume a complete workflow

The main menu separates individual stages from complete workflows:

```text
1.  Prepare a new enzyme-redesign project
2.  Run a prepared RFdiffusion2 project
3.  Run LigandMPNN on completed backbones
4.  Run kinetic triage on protein sequences
5.  Run Boltz-2 on strict kinetic survivors
6.  Generate or refresh Boltz-2 result reports
7.  Run or resume the complete workflow
8.  Run or resume a controlled native pilot
9.  Run or resume the full 500-sequence release demonstration
10. Generate publication and quality-control exports
11. Export or verify a functional-site mapping
12. Inspect prepared projects
13. Open the results dashboard
14. Check local installations
15. Exit
```

Use option 7 when you want ENZTRA to continue a normal prepared project through its remaining stages. Use option 9 for the exact 50 × 10 release-demonstration workflow and its release checks.

If a run is interrupted:

1. return to the ENZTRA repository;
2. activate the `enztra` environment;
3. launch `enztra`;
4. select option 7 or 9, as appropriate;
5. select the same project; and
6. inspect the resume summary before confirming.

ENZTRA validates completed-stage artifacts and resumes at the first incomplete stage. It does not checkpoint the middle of an external GPU inference operation. If interruption occurs during one backbone, that external operation may need to restart.

---

## 11. Monitor a long run

### 11.1 Watch the ENZTRA terminal

RFdiffusion2 prints a line such as:

```text
Making design 0 of 0:50: .../run_release-500b_cond0_0
```

During sampling, progress advances through 100 flow-matching steps. A long interval before the first visible step can reflect model loading, input preparation, or the first large inference calculation.

### 11.2 Check the process

In a second terminal:

```bash
ps -eo pid,etime,%cpu,%mem,state,cmd \
  | grep -E 'run_inference|boltz|apptainer' \
  | grep -v grep
```

An active Python process with a running state and substantial CPU usage is evidence that the process is alive, even if the main terminal is temporarily quiet.

### 11.3 Check the GPU

```bash
watch -n 2 nvidia-smi
```

Press `Ctrl+C` to leave `watch`; this does not stop ENZTRA when ENZTRA is running in another terminal.

GPU memory allocation with briefly low utilization does not by itself prove a hang. Combine GPU information with process state, elapsed time, terminal output, and output-file timestamps.

### 11.4 Count generated RFdiffusion2 structures

Replace the project name as needed:

```bash
find "$HOME/Applications/ENZTRA/jobs/release-500b/outputs/rfdiffusion2" \
  -maxdepth 1 -type f -name '*-atomized-bb-False.pdb' \
  | wc -l
```

The output directory may not exist until the external tool creates its first output. A missing directory early in the first design is not, by itself, proof of failure.

### 11.5 Inspect the newest files

```bash
find "$HOME/Applications/ENZTRA/jobs/release-500b" \
  -type f -printf '%T@ %TY-%Tm-%Td %TH:%TM:%TS %p\n' \
  | sort -nr \
  | head -20
```

### 11.6 Hardware-specific timing example

On the workstation used for the reported demonstrations, the recorded RFdiffusion2 stage durations were approximately:

- `release-500` (150–180 residues): 23,057 seconds, about 6.4 hours;
- `release-500b` (410 residues): 73,335 seconds, about 20.4 hours.

These are observations from one machine and software setup, not promised runtimes. The 410-residue designs took roughly 23 minutes per backbone on that system.

---

## 12. Understand the outputs

Each project is saved under:

```text
jobs/<project-name>/
```

A completed release project can contain:

```text
jobs/release-500/
├── design_request.json
├── rfdiffusion2_spec.json
├── rfdiffusion2_execution.json
├── ligandmpnn_execution.json
├── inputs/
├── outputs/
├── raw/
├── kinetics.csv
├── selection.csv
├── boltz2/
├── publication/
├── release/
├── manifest.json
├── summary.json
├── status.json
└── workflow_state.json
```

Important files include:

- `design_request.json`: validated biological and design inputs;
- `rfdiffusion2_spec.json`: generated RFdiffusion2 specification;
- `rfdiffusion2_execution.json`: exact RFdiffusion2 command and validation context;
- `ligandmpnn_execution.json`: sequence-generation execution record;
- `outputs/functional_site_mapping.csv`: reference-to-design residue mapping;
- `kinetics.csv`: consolidated kinetic predictions;
- `selection.csv`: strict-gate decisions and fold changes;
- `boltz2/inputs/`: Boltz-2 input YAML files for kinetic survivors;
- `boltz2/results/`: predicted structures and confidence records;
- `publication/tables/master_ranking.csv`: combined candidate table;
- `publication/tables/quality_control.csv`: data-quality classifications;
- `publication/figures/`: generated summary figures;
- `publication/reproducibility_manifest.json`: file and provenance record;
- `release/qualification_report.json`: release-level integration checks;
- `workflow_state.json`: stage-completion and resume state.

Do not delete state or execution records merely to make the folder look smaller. They are central to reproducibility.

---

## 13. Verify a completed demonstration

### 13.1 Read the qualification report

```bash
cd "$HOME/Applications/ENZTRA"
python -m json.tool jobs/release-500/release/qualification_report.json
python -m json.tool jobs/release-500b/release/qualification_report.json
```

Look for:

```text
"qualification_status": "passed"
"exact_500_sequence_plan": true
"failed_checks": []
```

### 13.2 Read the publication summary

```bash
python -m json.tool jobs/release-500/publication/publication_summary.json
python -m json.tool jobs/release-500b/publication/publication_summary.json
```

### 13.3 Count sequences in the master table

```bash
python - <<'PY'
import csv
from pathlib import Path

for project in ("release-500", "release-500b"):
    path = Path("jobs") / project / "publication" / "tables" / "master_ranking.csv"
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    print(project, "candidates:", len(rows))
PY
```

### 13.4 Inspect selection ranges

```bash
python - <<'PY'
import csv
from pathlib import Path

for project in ("release-500", "release-500b"):
    path = Path("jobs") / project / "selection.csv"
    rows = list(csv.DictReader(path.open()))
    print(f"\n{project}")
    for field in ("kcat_s", "km_mm", "kcat_fold_change", "km_fold_change"):
        values = [float(row[field]) for row in rows if row.get(field)]
        print(f"  {field}: min={min(values):.6g}, max={max(values):.6g}")
PY
```

### 13.5 Generate or refresh publication outputs

From the ENZTRA menu, use:

```text
10. Generate publication and quality-control exports
```

This is useful after a completed stage or after report-generation code has been updated. It should not alter the underlying scientific predictions.

---

## 14. Common problems

### `enztra: command not found`

Activate the environment and reinstall from the repository root:

```bash
conda activate enztra
cd "$HOME/Applications/ENZTRA"
python -m pip install -e ".[test]"
command -v enztra
```

### ENZTRA imports the wrong version

```bash
python -c "import enztra, pathlib; print(enztra.__version__); print(pathlib.Path(enztra.__file__).resolve())"
```

The printed source path should point to the intended repository or environment.

### Apptainer cannot access the GPU

Test it outside ENZTRA:

```bash
apptainer exec --nv /absolute/path/to/bakerlab_rf_diffusion_aa.sif nvidia-smi
```

If this fails, fix Apptainer/NVIDIA integration before retrying ENZTRA.

### RFdiffusion2 appears frozen at `Making design ...`

Use all of the following before deciding that it is hung:

- `ps` process state and elapsed time;
- CPU utilization;
- `nvidia-smi` process and memory information;
- newest output timestamps;
- whether flow-matching steps begin within the normal hardware-specific interval.

A 410-residue first design can take much longer than a 150–180-residue design.

### No `outputs/rfdiffusion2` directory exists yet

Confirm that you used the full absolute project path. If the process is active and still on its first design, the directory or final PDB may not yet have been written.

### A completed run has zero survivors

This can be a correct result. Check that:

- all planned candidates were evaluated;
- QC failures are absent or explained;
- reference predictions are present;
- strict inequalities were applied; and
- the qualification report passed.

Boltz-2 should be skipped when there are no strict kinetic survivors.

### A residue token is rejected

Use the full format:

```text
CHAIN:ONE_LETTER_AMINO_ACID_POSITION
```

For example:

```text
A:S215
```

Confirm the residue identity and numbering directly in the exact PDB file being used.

### A project was interrupted

Do not create a second project with the same scientific intent merely to continue. Relaunch ENZTRA, select option 7 or 9, select the original project, and review the resume report.

---

## 15. Reproducibility checklist

Before reporting or depositing a run, record:

- [ ] ENZTRA version and Git commit or release tag;
- [ ] operating system;
- [ ] Python and Conda versions;
- [ ] GPU model, driver, and reported CUDA compatibility;
- [ ] Apptainer version;
- [ ] RFdiffusion2 version/commit, container, and checkpoint identity;
- [ ] LigandMPNN version/commit;
- [ ] DLKcat version/commit and model files;
- [ ] CatPred version/commit and checkpoint set;
- [ ] Boltz version/commit and model version;
- [ ] exact reference PDB and its checksum;
- [ ] substrate CCD code and stereochemical SMILES;
- [ ] functional residues and preserved atoms;
- [ ] ORI mode and resulting coordinates;
- [ ] design length, backbone count, and sequences per backbone;
- [ ] random/deterministic settings recorded by the execution files;
- [ ] complete `design_request.json` and execution records;
- [ ] `workflow_state.json`, qualification report, and reproducibility manifest;
- [ ] complete candidate denominator, including rejected candidates;
- [ ] a clear statement that all predictions require experimental validation.

For reviewer-friendly archiving, exclude credentials, caches, local Conda environments, model weights, containers, third-party installation directories, and files whose redistribution is prohibited. Keep the validated inputs, outputs, mappings, prediction tables, structures permitted for redistribution, manifests, logs, and checksums.

---

## Getting help

When opening a GitHub issue, include:

- ENZTRA version;
- operating system and GPU;
- the selected menu option;
- the project review or preflight summary;
- the relevant error message;
- `status.json` and the relevant stage execution record; and
- the smallest reproducible example that does not expose confidential data.

Do not publish credentials, private paths containing sensitive information, access tokens, or licensed model weights.

---

## For questions
please contact the corresponding authors:

Corresponding authors:
    Olanrewaju Ayodeji Durojaye;
    Rachid Daoud

Institution:
    Chemical and Biochemical Sciences, Green Process Engineering,
    University Mohammed VI Polytechnic,
    43150 Ben Guerir, Morocco

Email:
    olanrewaju.ayodeji-durojaye-ext@um6p.ma;
    rachid.daoud@um6p.ma
