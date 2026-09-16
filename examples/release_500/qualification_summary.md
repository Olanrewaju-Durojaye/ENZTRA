# ENZTRA v1.0.0 release qualification

- Qualification status: passed
- Exact release plan: 50 RFdiffusion2 backbones × 10 LigandMPNN sequences
- Protein sequences evaluated kinetically: 500
- Strict kinetic survivors: 14
- Survivors evaluated and ranked structurally: 14
- Failed integration checks: 0
- GitHub release ready: yes

The expensive model execution was completed with `1.0.0rc1`. Release
candidates `rc2` and `rc3` changed only reporting labels: first preserving
backbone provenance, then adopting the authoritative `Priority N` labels from
the final-priority table. They did not alter model inputs, predictions,
selection decisions, structural scores, or ranking calculations.

Computational predictions are hypotheses and require experimental validation.
