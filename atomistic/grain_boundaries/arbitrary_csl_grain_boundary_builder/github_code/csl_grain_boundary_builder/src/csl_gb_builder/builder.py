"""aimsgb-backed structure generation."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict

from .config import GBConfig


def _load_backend() -> Any:
    try:
        from aimsgb import Grain, GrainBoundary  # type: ignore
        from pymatgen.core import Lattice, Structure  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "aimsgb and pymatgen are required. Install dependencies with: "
            "python -m pip install -r requirements.txt"
        ) from exc
    return Grain, GrainBoundary, Lattice, Structure


def build_grain_boundary(config: GBConfig) -> Dict[str, Any]:
    """Construct a periodic FCC CSL grain boundary and write it to disk."""
    Grain, GrainBoundary, Lattice, Structure = _load_backend()
    lattice = Lattice.cubic(config.lattice_constant)
    fcc_sites = ([0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0])
    conventional = Structure(lattice, [config.element] * 4, fcc_sites)
    initial_grain = Grain.from_sites(conventional.sites)

    gb = GrainBoundary(
        list(config.rotation_axis), config.sigma, list(config.plane), initial_grain,
        config.uc_a, config.uc_b,
    )
    structure = Grain.stack_grains(
        gb.grain_a,
        gb.grain_b,
        vacuum=config.vacuum,
        gap=config.gap,
        direction=gb.direction,
        delete_layer=config.delete_layer,
        tol=config.layer_tolerance,
        to_primitive=config.to_primitive,
    )
    config.output_path.parent.mkdir(parents=True, exist_ok=True)
    structure.to(filename=str(config.output_path), fmt=config.output_format)

    metadata = {
        "element": config.element,
        "sigma": gb.sigma,
        "rotation_axis": list(config.rotation_axis),
        "plane": list(config.plane),
        "theta_degrees": gb.theta,
        "csl_matrix": gb.csl.tolist(),
        "stacking_direction": gb.direction,
        "atom_count": len(structure),
        "output_path": str(config.output_path),
        "output_format": config.output_format,
    }
    metadata_path = config.output_path.with_name(config.output_path.name + ".metadata.json")
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata
