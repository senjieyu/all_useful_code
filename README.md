# All Useful Code

Curated, de-duplicated collection of reusable scripts from five personal code folders.  The previous repository layout mixed source code with simulation outputs, generated structures, binary archives, and IDE caches; this revision keeps the reusable source code in a compact catalog.

## Layout

- `atomistic/format_conversion` — converters between CIF, POSCAR/VASP, XYZ, LAMMPS, CFG, and related formats.
- `atomistic/grain_boundaries` — grain-boundary construction and energy workflows.
- `atomistic/grain_boundaries/arbitrary_csl_grain_boundary_builder` — configuration-driven builder for arbitrary CSL grain-boundary parameters; it includes a validated Σ17{223}[110] Cu Atomsk/LAMMPS example and usage notes.
- `atomistic/gsfe` — generalized stacking-fault (GSFE) construction and analysis.
- `atomistic/material_properties` — EOS and elastic-property workflows.
- `atomistic/simulation_workflows` — VASP, LAMMPS, EAM/SUS2, relaxation, training, and sampling helpers.
- `atomistic/twins_and_stacking` — twinning, stacking, and supercell construction.
- `machine_learning` — regression, classification, feature selection, Bayesian optimization, and visualization utilities.
- `python_learning` — standalone Python learning exercises and small applications.
- `utilities_and_legacy` — miscellaneous utilities retained for reference.
- `visualization` — standalone plotting helpers.

## Naming and traceability

All catalog paths use lowercase English names separated by underscores. Exact duplicate files are stored once.  When two different scripts would have the same normalized name, a deterministic suffix such as `__002` is added.

`.metadata/source_manifest.csv` records the original file path and filename, curated path, SHA-256, and byte size for every retained source.  Treat this repository as a curated script catalog: scripts that originally belonged to a larger project may need their original configuration files, data, or environment before they can be run.

## Deliberately excluded

Generated calculation outputs, structure datasets, archives, binaries, IDE settings/cache files, and duplicate copies were omitted.  They are inappropriate for a source-code repository and remain recoverable from the prior Git history or the original desktop folders.

## Validation

489 of 493 Python files passed a syntax-compilation check. The four retained legacy/reference files that are not executable Python are listed in `.metadata/validation.md`.

## Basic checks

The Python scripts are intended to be checked with:

```powershell
Get-ChildItem -Recurse -Filter *.py | ForEach-Object { python -m py_compile $_.FullName }
```

This checks syntax only; individual scripts have external scientific-computing dependencies and may require a configured VASP/LAMMPS/ML environment.
