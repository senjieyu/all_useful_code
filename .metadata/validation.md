# Catalog validation

Validation completed on 2026-10-10.

- 493 Python files were syntax-checked with `python -m py_compile`.
- 489 files passed syntax compilation.
- 4 files are preserved as legacy/reference material but are not executable Python in their current form:
  - `atomistic/simulation_workflows/test_validate_train.py` — unfinished triple-quoted note.
  - `machine_learning/ovr_mvm_classification.py` — prose saved with a `.py` extension.
  - `python_learning/atmosk_merge.py` — Atomsk command text saved with a `.py` extension.
  - `python_learning/test__005.py` — incomplete learning exercise.

The catalog intentionally preserves these four files for provenance. They should be repaired or reclassified before use as executable scripts.
