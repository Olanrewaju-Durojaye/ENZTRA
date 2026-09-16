# Contributing to ENZTRA

Thank you for helping improve ENZTRA. Please open an issue before a large
change so its scientific assumptions, scope, and validation plan are explicit.

## Development setup

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[test]"
pytest
```

Pull requests should include focused tests and update the README or changelog
when user-visible behavior changes. Do not commit model weights, Conda
environments, private structures, generated `jobs/` data, or local
`enztra.config.json` files.

## Scientific changes

Changes to thresholds, ranking equations, residue mappings, model inputs, or
acceptance rules must document the rationale and add boundary tests. Structural
confidence must not silently replace the strict kinetic gate. Computational
predictions must be described as hypotheses requiring experimental validation.

## External tools

ENZTRA integrates independent local installations but does not redistribute
their source code or weights. Contributions must preserve those boundaries and
the upstream projects' licenses and citation requirements.
