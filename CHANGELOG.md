# Changelog

## 1.0.1 — 2026-09-16

- Fix RFdiffusion2 length construction in protein-guidepost mode by retaining
  the requested scaffold span in the contig and offsetting ContigMap's length
  compatibility constraint by the number of functional-residue guideposts.
- Quote exact Hydra length scalars and add regression coverage for both exact
  lengths and ranges, preventing incompatible-contig sampling loops.

## 1.0.0 — 2026-09-16

- Promote ENZTRA to its first stable release after the exact native
  50-backbone × 10-sequence demonstration passed all integration checks.
- Verify 500/500 sequences reached kinetic evaluation and all 14 strict
  kinetic survivors reached Boltz-2 structural evaluation and ranking.
- Confirm complete publication and quality-control exports with no failed
  release-qualification checks and `github_release_ready: true`.
- Retain the release-candidate scientific calculations unchanged while using
  traceable `Priority N` labels consistently in final publication figures and
  `master_ranking.csv`.
- Include citation, contribution, security, release-demonstration, and GitHub
  release documentation suitable for a public source repository.

## 1.0.0rc3 — 2026-09-16

- Label publication kinetic and Boltz-2 figures by the authoritative combined
  `final_priority_rank`, using `Priority 1`, `Priority 2`, and so on.
- Add `publication_label` to `master_ranking.csv` so every plotted survivor can
  be traced directly to its complete design ID, sequence, selected structure,
  QC status, and deposited Boltz-2 folder.
- Continue retaining full design identifiers in SVG tooltips and data exports.

## 1.0.0rc2 — 2026-09-16

- Fix generated-design labels in kinetic and Boltz-2 figures so compact labels
  retain both backbone and sequence provenance, for example
  `bb011 · seq0007`, instead of the ambiguous shared sequence suffix alone.
- Preserve the complete original identifier in SVG hover text and every
  CSV/JSON export.

## 1.0.0rc1 — 2026-09-16

- Add an exact 50-backbone × 10-sequence release-demonstration profile for the
  final 500-design qualification run.
- Add installation and calibrated disk-capacity preflight checks, including a
  protected free-space reserve and explicit blockers before GPU execution.
- Add `enztra demonstrate` and a guided terminal action for preparing,
  running, interrupting, and resuming the release demonstration.
- Generate JSON and Markdown qualification reports that certify planned and
  observed handoffs and expose whether the tested code is ready for v1.0.0.
- Add GitHub-facing citation, contribution, security, and release-demonstration
  documentation while retaining the validated v0.13.0 scientific workflow.
- Mark this package as a release candidate: the stable version is promoted
  only after the complete 500-sequence run passes every integration check.

## 0.13.0 — 2026-09-16

- Promote the controlled-pilot release candidate after a complete native
  10-backbone × 5-sequence run passed every integration check.
- Confirm exact recovery of all 10 requested RFdiffusion2 backbones, all 50
  LigandMPNN sequences, and all 50 kinetic candidates.
- Confirm two strict kinetic survivors were both processed and ranked by
  Boltz-2 and all 50 candidates reached publication export.
- Record the observed stage timings and verified pilot report under
  `examples/pilot_50/` for transparent release qualification.
- Retain resumable execution, controlled-pilot diagnostics, opt-in CatPred
  uncertainty review thresholds, and non-overriding structural QC warnings.

## 0.13.0rc1 — 2026-09-16

- Add a controlled native-pilot runner with installation, disk-space, pilot
  size, and backbone-diversity preflight checks.
- Recommend 50–100 pilot designs and record the exact plan in
  `pilot/preflight.json` without blocking intentionally smaller diagnostic runs.
- Generate machine-readable and Markdown integration reports comparing planned
  versus observed backbone, sequence, kinetic, structural, and publication counts.
- Record available stage durations and expose incomplete or mismatched handoffs
  instead of relying only on external-tool exit codes.
- Add `enztra pilot`, a prepare-only mode, and a guided terminal action for
  starting or resuming the controlled native pilot.
- Remove the uncalibrated CatPred Km SD `0.5` warning as a universal default;
  uncertainty remains in every export and threshold warnings are now opt-in.
- Retain structural review warnings as prioritization aids that never reverse
  the strict kinetic decision.

## 0.12.0 — 2026-09-15

- Add a consolidated master ranking that joins sequence provenance, predicted
  kinetics, strict-gate decisions, CatPred uncertainty, and the selected
  Boltz-2 model for every candidate.
- Keep strict kinetic acceptance, structural rank, and final prioritization as
  separate fields; the optional combined score cannot rescue a rejected design.
- Define the combined priority score transparently as the geometric mean of
  kcat fold improvement, inverse Km fold, confidence score, ligand ipTM, and
  complex pLDDT.
- Add per-candidate QC results for invalid sequences or kinetics, missing
  expected structures, missing/high Km uncertainty, and low structural metrics.
- Generate publication-ready CSV, JSON, Markdown, and standalone SVG exports,
  including complete zero-survivor reports.
- Add a reproducibility manifest containing ENZTRA version, requests, workflow
  state, commands, input/output inventory, sizes, and SHA-256 change identifiers.
- Add `enztra publication-export`, a guided terminal action, and an automatic
  publication stage at the end of resumable workflows.

## 0.11.0 — 2026-09-15

- Add a resumable end-to-end runner connecting RFdiffusion2, LigandMPNN,
  DLKcat/CatPred strict kinetic triage, and survivor-only Boltz-2 validation.
- Persist an authoritative `workflow_state.json` with per-stage status,
  attempts, timestamps, return codes, and failure messages.
- Validate durable output files before skipping a completed stage, preventing
  an incomplete or corrupt result from being treated as successful.
- Stop immediately at the first failed stage while leaving downstream stages
  pending; resuming retries the failed stage without repeating validated work.
- Convert abandoned `running` states to `interrupted` after a stopped process
  and support safe continuation from the terminal menu or CLI.
- Treat an empty strict-survivor queue as a successful scientific outcome and
  mark Boltz-2 as skipped rather than failed.
- Add `enztra workflow` and the guided **Run or resume the complete workflow**
  menu action, including Boltz-2 MSA and sampling controls.

## 0.10.0 — 2026-09-15

- Rank every Boltz-2 sample within its design and automatically select the best
  model using confidence score, ligand ipTM, and complex pLDDT.
- Rank the best models across strict kinetic survivors and combine structural
  confidence with kinetic values and kinetic survivor rank.
- Write `boltz2/results/reports/best_models.csv` and a permanent
  `boltz2_confidence.svg` plot while retaining the complete all-sample table.
- Add a GPU-free terminal action and `enztra boltz2-report` command to rebuild
  reports from completed Boltz-2 outputs without rerunning prediction.
- Display Boltz-2 summary cards and the structurally ranked best-model table in
  the local results dashboard.
- Prevent nested Boltz-2 summary files from appearing as false dashboard jobs.

## 0.9.1 — 2026-09-15

- Label only the reference and strict kinetic survivors in permanent plots;
  rejected candidates remain visible as unlabelled red points.
- Preserve every complete identifier in SVG hover text and all tabular outputs.
- Shorten long RFdiffusion2/LigandMPNN survivor labels with a middle ellipsis
  while retaining both distinguishing ends of the identifier.
- Prevent Km and kcat plot axes from extending into impossible negative values.

## 0.9.0 — 2026-09-15

- Generate a permanent, standalone kinetic-selection decision plot after every
  completed kinetic triage run.
- Save the SVG under `jobs/<job-name>/results/kinetic_selection.svg`, so the
  result remains available without starting the browser dashboard.
- Plot the reference, strict survivors, rejected designs, strict decision
  boundaries, and the qualifying upper-left region using no new dependency.
- Print the exact plot path when kinetic processing finishes while retaining
  the optional interactive dashboard.

## 0.8.1 — 2026-09-14

- Add guided kinetic intake for external FASTA or FASTA-formatted text files,
  as well as completed ENZTRA LigandMPNN projects.
- Accept either a one-record reference FASTA or an active reference PDB with
  guided protein-chain selection and sequence extraction.
- Normalize metadata-rich external headers into safe, unique design IDs while
  preserving every original header and key-value field for provenance.
- Detect identical protein sequences and offer either provenance-preserving
  deduplication or explicit retention of every input record.
- Connect imported sequences to the existing DLKcat/CatPred preparation and
  execution pipeline, producing the standard strict-gate outputs required by
  the Boltz-2 survivor stage.

## 0.8.0 — 2026-09-14

- Add a guarded Boltz-2 structural-validation stage restricted to designs in
  the strict kinetic survivor queue.
- Generate one schema-valid protein–ligand YAML input per survivor using the
  kinetic manifest's stereochemical substrate SMILES.
- Offer reliable single-sequence mode or optional ColabFold MSA-server mode.
- Show the exact Conda/Boltz GPU command and workload before confirmation.
- Deliberately omit affinity prediction and distinguish structural confidence
  from catalytic activity and experimental affinity.
- Require confidence JSON and mmCIF output for every requested sample, then
  consolidate confidence metrics and the best model per survivor.

## 0.7.4 — 2026-09-14

- Restrict LigandMPNN collection to the exact completed backbones recorded in
  the current execution plan, so unrelated or stale FASTA files cannot alter
  the result of a repeated run.
- Validate the generated sequence count independently for every backbone and
  report the precise missing or malformed output.
- Name consolidated sequences with a per-backbone counter and reject invalid
  amino-acid symbols instead of silently dropping affected records.
- Limit fixed-position comparison to the planned backbone set while requiring
  every planned backbone to be present in the TRB-derived mapping.

## 0.7.3 — 2026-09-14

- Export each validated original-to-scaffold functional-site correspondence to
  `outputs/functional_site_mapping.csv`, including fixed atom names.
- Preserve residue identities in the exported mapping when they are supplied
  by the user or available from the reference PDB.
- Compare LigandMPNN's generated fixed-position JSON with the TRB-derived
  scaffold mapping and reject sequence output if they differ.
- Add a terminal menu action that backfills or verifies the mapping from an
  existing completed TRB without rerunning RFdiffusion2.

## 0.7.2 — 2026-09-14

- Preserve the RFdiffusion2 contig as one comma-delimited string inside
  Hydra's list instead of allowing Hydra to split it into separate list items.
- Ensure `ContigMap` receives the requested length and all functional-site
  residue segments, rather than silently reading only the length element.
- Retain mandatory post-run TRB validation before LigandMPNN eligibility.

## 0.7.1 — 2026-09-14

- Call `run_inference.py` directly so Hydra receives the functional-atom map as
  one intact argument instead of passing it through a lossy nested shell.
- Validate every completed TRB inside the RFdiffusion2 container and require the
  expected number of mapped protein guideposts.
- Mark missing-guidepost runs as failed and prevent their backbones from being
  offered to LigandMPNN.

## 0.7.0 — 2026-09-14

- Add a separate guided LigandMPNN execution stage for completed RFdiffusion2
  backbone jobs.
- Reuse RFdiffusion2's ligand-aware MPNN v2 pipeline and preserve contig motif
  positions during sequence generation.
- Show the detected backbone count, sequences per backbone, total workload,
  output directory, and exact command before requesting confirmation.
- Consolidate generated sequences into a validated, uniquely named
  `outputs/ligandmpnn/designs.fasta` file for downstream kinetic prediction.
- Track LigandMPNN running, completed, failed, and interrupted states.

## 0.6.2 — 2026-09-14

- Generate guidepost contigs with one scaffold segment matching the requested
  designed-protein length, followed by all preserved functional-site residues.
- Prevent RFdiffusion2 from collapsing repeated `0-40` gaps into an impossible
  0–40-residue scaffold when the requested design length is larger.
- Add regression coverage for ranged design lengths in RFdiffusion2 execution.

## 0.6.1 — 2026-09-14

- Preserve escaped atom-dictionary quotes through RFdiffusion2's nested shell
  command so Hydra can parse multi-residue `contigmap.contig_atoms` mappings.
- Set `MKL_THREADING_LAYER=GNU` inside Apptainer to prevent Intel MKL and
  `libgomp` incompatibility on mixed Conda installations.
- Allow failed, preflight-ready projects to be selected and retried without
  recreating their scientific inputs.

## 0.6.0 — 2026-09-14

- Add guarded local RFdiffusion2 backbone execution for prepared PDB projects.
- Generate an official custom benchmark JSON from the validated motif, atom
  mapping, length range, cleaned PDB, ligand, and ORI preflight.
- Show the exact Apptainer command, requested backbone count, and output path
  before asking for explicit GPU confirmation.
- Save a reusable `rfdiffusion2_execution.json` command plan even when launch
  is cancelled.
- Track running, completed, failed, and interrupted RFdiffusion2 job states.
- Keep backbone generation separate from LigandMPNN sequence generation.

## 0.5.2 — 2026-09-14

- Add guided RFdiffusion2 ORI-token placement with existing-token,
  automatic C-alpha-centroid, and custom-coordinate modes.
- Default new PDB projects to automatic centroid placement when no ORI token
  is present.
- Write a cleaned `rfdiffusion2_input.pdb` containing only the selected protein
  chain, intended substrate, and one validated ORI token.
- Exclude unrelated HETATM records such as Hg ions, waters, and crystallization
  additives from the RFdiffusion2 execution input.
- Keep the original reference PDB unchanged for provenance.

## 0.5.1 — 2026-09-14

- Add a one-step default mode that preserves every available heavy atom for
  all selected functional-site residues.
- Retain advanced residue-by-residue atom selection for users who need custom
  constraints.
- Show a compact atom-preservation summary for the default mode while keeping
  the complete mapping in the saved project specifications.
- Preserve the user's functional-site residue order in generated RFdiffusion2
  specifications.

## 0.5.0 — 2026-09-14

- Add guided, residue-by-residue selection of the exact heavy atoms that
  RFdiffusion2 must preserve in the functional motif.
- Read available atom names directly from the selected PDB and reject unknown
  or missing selections rather than guessing chemically important atoms.
- Support `ALL` as an explicit choice when the complete residue must be fixed.
- Generate `rfdiffusion2_spec.json` containing the reference PDB, structural
  ligand label, unindexed motif residues, design-length constraint, and fixed
  atom mapping.
- Mark prepared PDB projects as RFdiffusion2 execution-ready after atom-level
  preflight succeeds. GPU inference remains a separately confirmed operation.

## 0.4.5 — 2026-09-14

- Add a guided **Designed protein length or range** prompt supporting exact
  values such as `200` and ranges such as `180-220`.
- Default the prompt to the selected reference chain length while allowing the
  user to replace it explicitly.
- Validate range order and reject lengths above 2,000 residues.
- Store the normalized length constraint in `design_request.json` and display
  it during final project review.

## 0.4.4 — 2026-09-14

- Detect Boltz's generic `LIG` label before requesting substrate metadata.
- Rename the prompt to **Substrate CCD code** and explicitly tell users to enter
  the original scientific code, such as `IHP`, rather than `LIG`.
- Reject `LIG` as the scientific substrate code when it is clearly being used
  as a generic Boltz structure label.
- Display the resolved mapping in project review, for example
  `LIG B:1 → CCD IHP`.

## 0.4.3 — 2026-09-14

- Recognize Boltz-generated PDB structures whose ligand residue name is the
  generic `LIG` rather than the substrate's original CCD code.
- Automatically map a single `LIG` residue to the user-supplied substrate code
  while retaining stereochemical SMILES as the authoritative chemical identity.
- Request a chain/residue selector only when multiple generic `LIG` residues
  make the structural ligand ambiguous.
- Display successful generic-ligand mappings as reference notes instead of
  false missing-CCD warnings.

## 0.4.2 — 2026-09-14

- Rename catalytic-residue input to **functional-site residues to preserve** so
  catalytic, substrate-interacting, and other constrained pocket residues are
  represented without implying that every selected residue is catalytic.
- Add identity-aware residue notation such as `A:H17` while retaining legacy
  position-only entries such as `A:17`.
- Compare expected identities with FASTA sequences and PDB ATOM records and warn
  about possible inactive or engineered mutations.
- Check that functional-site positions exist and reject duplicate positions.
- Check whether the supplied ligand CCD code occurs in PDB HETATM records.
- Warn that experimental structures may contain mutations, missing residues, or
  unresolved loops, and require active-reference confirmation before saving.
- Clarify the different roles of descriptive substrate names, CCD codes, and
  stereochemical SMILES.

## 0.4.1 — 2026-09-14

- Replace the browser design form with a guided, BOLTRA-style terminal workflow.
- Make `enztra` with no arguments open a numbered main menu.
- Add one-question-at-a-time project preparation with immediate re-prompting for
  invalid paths, identifiers, chains, ligand codes, residues, and design counts.
- Detect and display chains from reference PDB files.
- Show a complete project review and require confirmation before writing files.
- Add prepared-project inspection and direct dashboard/doctor access from the menu.
- Restore the browser application to results visualization only.

## 0.4.0 — 2026-09-14

- Add a general browser interface for preparing enzyme-design requests.
- Accept a one-protein FASTA or PDB reference, protein chain, substrate name,
  stereochemical SMILES, ligand code, and catalytic-residue constraints.
- Add configurable RFdiffusion2 backbone and LigandMPNN sequence counts.
- Validate identifiers, reference contents, chains, residues, and job size before
  writing a reproducible `design_request.json`.
- Save prepared jobs locally without launching models prematurely.
- Add prepared-job history and separate New design and Results views.

This browser-input approach was superseded by the terminal workflow in 0.4.1.

## 0.3.0 — 2026-09-13

- Add `enztra serve` for a local, read-only results dashboard.
- Discover completed jobs recursively beneath a configurable jobs directory.
- Add an interactive Km-versus-kcat decision map with reference boundaries,
  strict-survivor shading, point tooltips, and linear/log axes.
- Add summary cards, searchable selection details, and the Boltz-2 queue.
- Serve the dashboard using only Python's standard library, with no new runtime
  dependencies.

## 0.2.1 — 2026-09-14

- Collect CatPred's processed results from its configured results directory
  instead of the raw prediction written beside the input CSV.
- Give CatPred inputs job-specific names to prevent cross-job collisions.
- Add `--reuse-raw` recovery so completed predictions can be combined without
  rerunning either model.

## 0.2.0 — 2026-09-13

- Add configurable discovery of local RFdiffusion2, DLKcat, CatPred, and
  Boltz-2 installations.
- Add `enztra doctor` with read-only path, environment, and GPU checks.
- Add strict multi-record protein FASTA parsing.
- Add equivalent DLKcat and CatPred batch-input generation.
- Add isolated Conda execution adapters for DLKcat and GPU CatPred.
- Add prediction row-count, sequence, and identity validation.
- Add `enztra kinetics` and `--prepare-only` workflows.
- Add unified kinetic results, strict selection output, summary, and Boltz-2
  survivor queue.

## 0.1.0 — 2026-09-13

- Add validated kinetic records and strict reference-relative selection.
- Rank only strict survivors by predicted catalytic efficiency.
- Add the illustrative 1DKP–IHP reference case and automated tests.
