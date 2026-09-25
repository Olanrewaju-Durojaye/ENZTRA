# ENZTRA

<img width="1919" height="820" alt="enztra-logo" src="https://github.com/user-attachments/assets/7c5789cb-a778-4c6b-beea-41ffbc0d9430" />

**ENZTRA** (**EN**zyme redesign and **TR**iage **A**utomation) is a
reference-guided computational workflow for enzyme redesign, kinetic triage,
structural assessment, and reproducible reporting.

ENZTRA connects separately installed scientific tools without redistributing
their source code, model weights, or containers:

1. **RFdiffusion2** generates motif-constrained protein backbones.
2. **LigandMPNN** designs ligand-aware sequences while retaining mapped
   functional positions.
3. **DLKcat** predicts turnover number (`kcat`).
4. **CatPred** predicts Michaelis constant (`Km`) and associated uncertainty.
5. **Boltz-2** models candidates that pass the strict kinetic gate.

> ENZTRA predictions are computational hypotheses. They do not establish
> catalytic activity, binding, or experimental improvement. Selected candidates
> require experimental validation.

## Workflow

```text
Reference enzyme and substrate
            |
            v
Functional-site mapping and validated project specification
            |
            v
RFdiffusion2 backbones -> LigandMPNN sequences
            |
            v
DLKcat kcat + CatPred Km predictions
            |
            v
Strict kinetic gate
            |
            +-- rejected designs retained with provenance
            |
            +-- strict survivors -> Boltz-2 structural assessment
                                      |
                                      v
                         QC, ranking, and publication export
```

Every candidate remains linked to its reference constraints, parent backbone,
sequence, kinetic record, structural result, and final classification.

## Strict kinetic gate

A design advances to structural assessment only when both conditions are met:

```text
kcat_design > kcat_reference AND Km_design < Km_reference
```

The inequalities are strict: equality does not pass. Catalytic efficiency
(`kcat/Km`) can rank candidates that already pass, but it cannot rescue a design
that fails either threshold. Quality-control status is reported separately from
the kinetic decision.

## Core capabilities

- guided preparation of PDB- or FASTA-based enzyme-redesign projects;
- identity-aware functional-residue and preserved-atom constraints;
- validated RFdiffusion2-to-LigandMPNN residue mapping;
- native and externally generated sequence intake;
- coordinated DLKcat and CatPred prediction with identity and row-count checks;
- survivor-only Boltz-2 input generation and confidence reporting;
- resumable execution with validation of completed stages;
- stable candidate identifiers and retained intermediate records;
- consolidated rankings, QC classifications, figures, reports, and manifests;
- read-only installation checks and results dashboard; and
- controlled pilot and release-qualification reporting.

## Documentation

The complete beginner-friendly guide is provided in
**[TUTORIAL.md](TUTORIAL.md)**. It contains:

- platform requirements and ENZTRA installation;
- installation and configuration of external tools;
- dependency and GPU checks;
- preparation of the 1DKP–IHP phytase reference;
- the recommended smoke test;
- exact input conditions for `release-500` and `release-500b`;
- complete workflow execution and resumption;
- long-run monitoring;
- output interpretation and qualification checks;
- common-problem diagnosis; and
- a reproducibility checklist.

For release history, see **[CHANGELOG.md](CHANGELOG.md)**. For instructions on
contributing, see **[CONTRIBUTING.md](CONTRIBUTING.md)**.

## Current release

**ENZTRA v1.0.1** is the current stable release. It retains the native v1.0.0
qualification and corrects RFdiffusion2 length accounting when protein
guideposts are present.

The qualified `release-500` demonstration completed the planned 50-backbone x
10-sequence design: 500 sequences reached kinetic evaluation, 14 strict
survivors proceeded to Boltz-2, publication exports completed, and no
integration check failed. The matched full-length `release-500b` demonstration
also evaluated 500 sequences; none passed the strict kinetic gate, so structural
assessment was correctly skipped. These are reported computational outcomes,
not experimental validation.

Concise verification records are provided under `examples/` where distributed.
The accompanying curated demonstration data are available from Zenodo:

<https://doi.org/10.5281/zenodo.22836909>

## Getting started

Start with **[TUTORIAL.md](TUTORIAL.md)** and complete its installation checks
and one-design smoke test before launching a pilot or 500-sequence workflow.
Running ENZTRA without arguments opens the guided terminal interface:

```bash
enztra
```

Command-specific help is available with:

```bash
enztra --help
```

## Outputs and reproducibility

Depending on the completed stages, ENZTRA retains validated project
specifications, execution records, functional-site mappings, generated
backbones and sequences, kinetic predictions, strict-gate decisions, structural
confidence summaries, workflow state, QC tables, figures, publication reports,
and reproducibility manifests.

The workflow records negative outcomes as valid results. If no candidate passes
both kinetic thresholds, the project completes normally and records Boltz-2 as
skipped rather than treating the absence of survivors as an execution error.

## Software and data availability

ENZTRA v1.0.1 is available from:

<https://github.com/Olanrewaju-Durojaye/ENZTRA>

The curated demonstration data are available from:

<https://doi.org/10.5281/zenodo.22836909>

Third-party programs must be installed separately from their official sources
and cited according to their respective requirements.

## Citation

If ENZTRA contributes to published work, cite the ENZTRA software or associated
article when available and cite every external scientific tool used in the
analysis. Repository citation metadata are provided in
**[CITATION.cff](CITATION.cff)**.

## License

ENZTRA's original orchestration code is released under the
**[MIT License](LICENSE)**. RFdiffusion2, LigandMPNN, DLKcat, CatPred, Boltz-2,
their model weights, and their dependencies retain their respective licenses
and citation requirements.
