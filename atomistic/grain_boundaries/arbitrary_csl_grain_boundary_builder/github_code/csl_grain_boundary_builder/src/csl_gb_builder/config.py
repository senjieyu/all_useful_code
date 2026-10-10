"""JSON configuration parsing and validation for CSL grain-boundary jobs."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Tuple


class ConfigurationError(ValueError):
    """Raised when a grain-boundary configuration is incomplete or invalid."""


def _vector(value: Any, field: str) -> Tuple[int, int, int]:
    if not isinstance(value, list) or len(value) != 3:
        raise ConfigurationError(f"{field} must be a three-integer list.")
    try:
        vector = tuple(int(component) for component in value)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"{field} must contain integers.") from exc
    if vector == (0, 0, 0):
        raise ConfigurationError(f"{field} cannot be [0, 0, 0].")
    return vector


def _section(data: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    value = data.get(name)
    if not isinstance(value, Mapping):
        raise ConfigurationError(f"Missing object section: {name}")
    return value


@dataclass(frozen=True)
class GBConfig:
    """All parameters required for one aimsgb grain-boundary construction."""

    element: str
    lattice_constant: float
    rotation_axis: Tuple[int, int, int]
    sigma: int
    plane: Tuple[int, int, int]
    uc_a: int
    uc_b: int
    vacuum: float
    gap: float
    delete_layer: str
    layer_tolerance: float
    to_primitive: bool
    output_path: Path
    output_format: str

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any], base_dir: Path) -> "GBConfig":
        material = _section(data, "material")
        boundary = _section(data, "grain_boundary")
        output = _section(data, "output")

        element = material.get("element")
        if not isinstance(element, str) or not element.strip():
            raise ConfigurationError("material.element must be a non-empty element symbol.")
        lattice_constant = float(material.get("lattice_constant_angstrom", 0))
        if lattice_constant <= 0:
            raise ConfigurationError("material.lattice_constant_angstrom must be positive.")

        sigma = int(boundary.get("sigma", 0))
        if sigma < 3:
            raise ConfigurationError("grain_boundary.sigma must be at least 3.")
        uc_a = int(boundary.get("uc_a", 1))
        uc_b = int(boundary.get("uc_b", 1))
        if min(uc_a, uc_b) < 1:
            raise ConfigurationError("grain_boundary.uc_a and uc_b must be positive.")

        output_path = output.get("path")
        if not isinstance(output_path, str) or not output_path.strip():
            raise ConfigurationError("output.path must be a non-empty path.")
        output_format = str(output.get("format", "poscar")).lower()
        if output_format not in {"poscar", "cif", "json", "cssr"}:
            raise ConfigurationError("output.format must be poscar, cif, json, or cssr.")

        to_primitive = boundary.get("to_primitive", False)
        if not isinstance(to_primitive, bool):
            raise ConfigurationError("grain_boundary.to_primitive must be true or false.")
        vacuum = float(boundary.get("vacuum_angstrom", 0.0))
        gap = float(boundary.get("gap_angstrom", 0.0))
        layer_tolerance = float(boundary.get("layer_tolerance_angstrom", 0.25))
        if vacuum < 0 or gap < 0 or layer_tolerance <= 0:
            raise ConfigurationError("vacuum, gap, and layer tolerance must be non-negative (tolerance > 0).")

        return cls(
            element=element.strip(),
            lattice_constant=lattice_constant,
            rotation_axis=_vector(boundary.get("rotation_axis"), "grain_boundary.rotation_axis"),
            sigma=sigma,
            plane=_vector(boundary.get("plane"), "grain_boundary.plane"),
            uc_a=uc_a,
            uc_b=uc_b,
            vacuum=vacuum,
            gap=gap,
            delete_layer=str(boundary.get("delete_layer", "0b0t0b0t")),
            layer_tolerance=layer_tolerance,
            to_primitive=to_primitive,
            output_path=(base_dir / output_path).resolve(),
            output_format=output_format,
        )


def load_config(path: str | Path) -> GBConfig:
    """Load one JSON configuration file and return a validated configuration."""
    config_path = Path(path).resolve()
    try:
        data = json.loads(config_path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise ConfigurationError(f"Configuration file does not exist: {config_path}") from exc
    except json.JSONDecodeError as exc:
        raise ConfigurationError(f"Invalid JSON in {config_path}: {exc.msg}") from exc
    if not isinstance(data, Mapping):
        raise ConfigurationError("Top-level JSON value must be an object.")
    return GBConfig.from_mapping(data, config_path.parent)
