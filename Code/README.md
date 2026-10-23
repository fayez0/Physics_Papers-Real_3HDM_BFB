# Weinberg 3HDM BFB: publication code and data

This folder is a curated subset of the supplied research archive. It contains the exact classifiers, independent original-field objective, numerical validation inputs and per-start records, selected exact algebraic certificates, and the source comparison driver. It does not claim independent verification of the theorem.

## Contents
- `code/`: exact branch solvers, decision and field certification, selected analytical auditors.
- `validation_scripts/`: reproducible input generation, BFGS validation, dispatch and result summarization.
- `validation/validation_scripts/`: archived exact inputs (398), numerical records, summaries, environment and external comparison records.
- `certificates/`: exact SOS and separating-functional data for analytical examples.
- `examples/`: reproducibility input fixtures.

## Commands (from this Code directory)
```bash
python validation_scripts/exact_dispatch.py examples/charged_flat.json
python validation_scripts/audit_revision.py
python validation_scripts/summarize_validation.py
```
A full fresh numerical run is more expensive and requires the dependencies listed in `requirements-tested.txt` and matching documented runtime choices; it can be started with `python validation_scripts/run_validation.py --workers 4`. Check the script help and output paths before running, as it may overwrite archived result files.

## Provenance and release status
Files originate in the author-provided research archive; only publication-relevant components were selected, and script references to the old validation subdirectory name were changed to `validation_scripts`. Existing validation outputs were not regenerated. External StableWein software is NOT included. Some supporting historical certificates are excluded.

**License pending author decision.** This folder has not been granted a public redistribution license by its preparation. Before public release, author must approve license, audit included assets, check exact reproduction claims and finalize DOI.
