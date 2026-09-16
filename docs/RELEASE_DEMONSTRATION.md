# ENZTRA v1 release demonstration

The stable v1 release is qualified by one native, resumable run of exactly 50
RFdiffusion2 backbones and 10 LigandMPNN sequences per backbone (500 protein
sequences total).

## Acceptance criteria

- Preflight reports no installation, plan, or disk-space blockers.
- RFdiffusion2 produces 50 validated backbone/TRB pairs.
- LigandMPNN produces 500 valid protein sequences.
- Kinetic prediction and strict selection process all 500 candidates.
- Every strict survivor reaches Boltz-2 and structural ranking.
- Publication exports contain all 500 candidate records.
- `release/qualification_report.json` records `github_release_ready: true`.

The number of strict survivors is a scientific outcome, not an acceptance
target. A zero-survivor run is valid if all candidates were processed and
Boltz-2 was correctly skipped.

## Reproducible command

```bash
enztra demonstrate \
  --config enztra.config.json \
  --job-dir jobs/release-500 \
  --calibration-job jobs/pilot-50 \
  --msa-mode single \
  --diffusion-samples 1 \
  --seed 42
```

Rerun the same command after interruption or recoverable failure. Validated
completed stages are preserved. Archive the qualification report,
`workflow_state.json`, publication reproducibility manifest, tables, plots,
and the exact source revision used for the run.
