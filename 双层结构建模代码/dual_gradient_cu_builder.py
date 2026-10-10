#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import ast
import json
import math
import os
import shutil
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable
try:
    import tomllib  # Python 3.11+
except ModuleNotFoundError:  # Python 3.10：使用下方内置简化 TOML 读取器
    tomllib = None


DEFAULT_ATOMSK = "auto"
DEFAULT_OUTPUT_ROOT = str(Path(__file__).with_name("output"))

LAYER_THICKNESS_NM = {
    "AB_A": 24.0,
    "AB_B": 20.0,
    "ABAB_A": 12.0,
    "ABAB_B": 10.0,
    "PURE": 44.0,
}

TWIN_N111 = {
    "T4": 4,
    "T6": 6,
    "T10": 10,
    "T15": 15,
    "T20": 20,
    "T24": 24,
    "T30": 30,
}

GRAIN_TEMPLATES = {
    "single": {"nx": 1, "ny": 1, "label": "single crystal"},
    "G2x2": {"nx": 2, "ny": 2, "label": "4 grains"},
    "G3x2": {"nx": 3, "ny": 2, "label": "6 grains"},
    "G4x2": {"nx": 4, "ny": 2, "label": "8 grains"},
    "G4x3": {"nx": 4, "ny": 3, "label": "12 grains"},
    "G5x3": {"nx": 5, "ny": 3, "label": "15 grains"},
}

PLAN_MODELS = [
    {
        "name": "REF",
        "family": "REF",
        "phase": "P0",
        "grain_a": "G3x2",
        "grain_b": "G3x2",
        "twin_a": "T15",
        "twin_b": "T15",
    },
    {
        "name": "TTG-1",
        "family": "TTG",
        "phase": "P0",
        "grain_a": "G3x2",
        "grain_b": "G3x2",
        "twin_a": "T10",
        "twin_b": "T15",
    },
    {
        "name": "TTG-2",
        "family": "TTG",
        "phase": "P0",
        "grain_a": "G3x2",
        "grain_b": "G3x2",
        "twin_a": "T10",
        "twin_b": "T20",
    },
    {
        "name": "TTG-3",
        "family": "TTG",
        "phase": "P2",
        "grain_a": "G3x2",
        "grain_b": "G3x2",
        "twin_a": "T10",
        "twin_b": "T24",
    },
    {
        "name": "TTG-6",
        "family": "TTG",
        "phase": "P2",
        "grain_a": "G3x2",
        "grain_b": "G3x2",
        "twin_a": "T6",
        "twin_b": "T20",
    },
    {
        "name": "TTG-7",
        "family": "TTG",
        "phase": "P0",
        "grain_a": "G3x2",
        "grain_b": "G3x2",
        "twin_a": "T6",
        "twin_b": "T24",
    },
    {
        "name": "TTG-8",
        "family": "TTG",
        "phase": "P2",
        "grain_a": "G3x2",
        "grain_b": "G3x2",
        "twin_a": "T6",
        "twin_b": "T30",
    },
    {
        "name": "TTG-12",
        "family": "TTG",
        "phase": "P2",
        "grain_a": "G3x2",
        "grain_b": "G3x2",
        "twin_a": "T4",
        "twin_b": "T30",
    },
    {
        "name": "GSG-1",
        "family": "GSG",
        "phase": "P0",
        "grain_a": "G4x2",
        "grain_b": "G3x2",
        "twin_a": "T15",
        "twin_b": "T15",
    },
    {
        "name": "GSG-2",
        "family": "GSG",
        "phase": "P0",
        "grain_a": "G4x2",
        "grain_b": "G2x2",
        "twin_a": "T15",
        "twin_b": "T15",
    },
    {
        "name": "GSG-3",
        "family": "GSG",
        "phase": "P2",
        "grain_a": "G4x3",
        "grain_b": "G3x2",
        "twin_a": "T15",
        "twin_b": "T15",
    },
    {
        "name": "GSG-4",
        "family": "GSG",
        "phase": "P0",
        "grain_a": "G4x3",
        "grain_b": "G2x2",
        "twin_a": "T15",
        "twin_b": "T15",
    },
    {
        "name": "GSG-6",
        "family": "GSG",
        "phase": "P2",
        "grain_a": "G5x3",
        "grain_b": "G2x2",
        "twin_a": "T15",
        "twin_b": "T15",
    },
    {
        "name": "DG-A",
        "family": "DG",
        "phase": "P0",
        "grain_a": "G4x2",
        "grain_b": "G2x2",
        "twin_a": "T10",
        "twin_b": "T20",
    },
    {
        "name": "DG-B",
        "family": "DG",
        "phase": "P2",
        "grain_a": "G4x2",
        "grain_b": "G2x2",
        "twin_a": "T6",
        "twin_b": "T24",
    },
    {
        "name": "DG-D",
        "family": "DG",
        "phase": "P2",
        "grain_a": "G4x3",
        "grain_b": "G2x2",
        "twin_a": "T10",
        "twin_b": "T20",
    },
    {
        "name": "DG-E",
        "family": "DG",
        "phase": "P0",
        "grain_a": "G4x3",
        "grain_b": "G2x2",
        "twin_a": "T6",
        "twin_b": "T24",
    },
    {
        "name": "DG-F",
        "family": "DG",
        "phase": "P2",
        "grain_a": "G4x3",
        "grain_b": "G2x2",
        "twin_a": "T4",
        "twin_b": "T30",
    },
    {
        "name": "DG-G",
        "family": "DG",
        "phase": "P2",
        "grain_a": "G5x3",
        "grain_b": "G2x2",
        "twin_a": "T10",
        "twin_b": "T20",
    },
    {
        "name": "DG-H",
        "family": "DG",
        "phase": "P2",
        "grain_a": "G5x3",
        "grain_b": "G2x2",
        "twin_a": "T6",
        "twin_b": "T24",
    },
    {
        "name": "DG-J",
        "family": "DG",
        "phase": "P0",
        "grain_a": "G4x2",
        "grain_b": "G3x2",
        "twin_a": "T10",
        "twin_b": "T20",
    },
]


@dataclass(frozen=True)
class BuildContext:
    atomsk: str
    root_dir: Path
    element: str = "Cu"
    a0: float = 3.633
    box_xy_nm: float = 18.0
    packing_scale: float = 1.0
    shiftx: float = 0.25
    shifty: float = 0.35
    twin_nx: int = 8
    twin_ny: int = 8


@dataclass(frozen=True)
class BlockSpec:
    grain_code: str
    twin_code: str | None
    thickness_nm: float
    box_xy_nm: float
    packing_scale: float
    shiftx: float
    shifty: float
    pattern_rotation_deg: float = 0.0
    grain_angle_offset_deg: float = 0.0
    start_label: int = 0
    start_direction: int = 1

    def block_id(self) -> str:
        twin_tag = self.twin_code or "FCC"
        dir_tag = "p" if self.start_direction >= 0 else "m"
        return (
            f"grain-{self.grain_code}_"
            f"twin-{twin_tag}_"
            f"thick-{format_float(self.thickness_nm)}nm_"
            f"xy-{format_float(self.box_xy_nm)}nm_"
            f"ps-{format_float(self.packing_scale)}_"
            f"phase-{self.start_label}{dir_tag}_"
            f"prot-{format_float(self.pattern_rotation_deg)}_"
            f"aoff-{format_float(self.grain_angle_offset_deg)}"
        )


def format_float(value: float) -> str:
    return f"{value:.2f}".replace(".", "p")


def nm_to_angstrom(value_nm: float) -> float:
    return value_nm * 10.0


def rotate_xy(x: float, y: float, cx: float, cy: float, degrees: float) -> tuple[float, float]:
    if degrees % 360 == 0:
        return x, y

    radians = math.radians(degrees)
    dx = x - cx
    dy = y - cy
    xr = dx * math.cos(radians) - dy * math.sin(radians)
    yr = dx * math.sin(radians) + dy * math.cos(radians)
    return cx + xr, cy + yr


def choose_orthogonal_single_twin_repeats(a0: float, box_xy_nm: float) -> tuple[int, int]:
    target = nm_to_angstrom(box_xy_nm)
    x_cell = a0 * math.sqrt(3.0 / 2.0)
    y_cell = a0 / math.sqrt(2.0)
    nx = max(1, int(round(target / x_cell)))
    ny = max(1, int(round(target / y_cell)))
    return nx, ny


def normalize_direction(direction: int) -> int:
    return 1 if direction >= 0 else -1


def choose_compatible_total_layers(target_layers: int, layers_per_twin: int) -> int:
    if target_layers < 1:
        raise ValueError("target_layers must be positive")
    if layers_per_twin < 1:
        raise ValueError("layers_per_twin must be positive")

    lower = (target_layers // layers_per_twin) * layers_per_twin
    upper = lower + layers_per_twin
    candidates = [value for value in (lower, upper) if value >= layers_per_twin]
    if not candidates:
        return layers_per_twin

    return min(candidates, key=lambda value: (abs(value - target_layers), -value))


def choose_centered_margins(total_layers: int, layers_per_twin: int) -> tuple[int, int]:
    if total_layers <= 0:
        raise ValueError("total_layers must be positive")
    if layers_per_twin < 2:
        raise ValueError("layers_per_twin must be at least 2")

    desired = layers_per_twin / 2.0
    best: tuple[tuple[float, int, int, int], int, int] | None = None

    for bottom_margin in range(1, layers_per_twin):
        for top_margin in range(1, layers_per_twin):
            interior_layers = total_layers - bottom_margin - top_margin
            if interior_layers < 0:
                continue
            if interior_layers % layers_per_twin != 0:
                continue
            score = (
                abs(bottom_margin - desired) + abs(top_margin - desired),
                abs(bottom_margin - top_margin),
                bottom_margin + top_margin,
                bottom_margin,
            )
            if best is None or score < best[0]:
                best = (score, bottom_margin, top_margin)

    if best is None:
        raise ValueError(
            f"Could not find nonzero edge margins for total_layers={total_layers}, "
            f"layers_per_twin={layers_per_twin}"
        )

    return best[1], best[2]


def plan_nanotwin_stack(
    total_layers: int,
    layers_per_twin: int,
    start_label: int = 0,
    start_direction: int = 1,
) -> dict[str, object]:
    if total_layers <= 0:
        raise ValueError("total_layers must be positive")
    if layers_per_twin < 2:
        raise ValueError("layers_per_twin must be at least 2")

    full_segments = total_layers // layers_per_twin
    top_margin = total_layers % layers_per_twin
    segment_lengths = ([layers_per_twin] * full_segments) + ([top_margin] if top_margin else [])
    if not segment_lengths:
        segment_lengths = [total_layers]

    current_label = start_label % 3
    current_direction = normalize_direction(start_direction)
    labels: list[int] = []
    boundary_positions: list[int] = [0]

    segment_index = 0
    layers_in_segment = 0
    current_segment_length = segment_lengths[segment_index]

    for layer_index in range(total_layers):
        labels.append(current_label)
        layers_in_segment += 1

        at_boundary = (
            segment_index < len(segment_lengths) - 1
            and layers_in_segment == current_segment_length
        )
        if at_boundary:
            boundary_positions.append(layer_index + 1)
            current_direction *= -1
            segment_index += 1
            layers_in_segment = 0
            current_segment_length = segment_lengths[segment_index]

        current_label = (current_label + current_direction) % 3

    boundary_spacings = [
        boundary_positions[idx + 1] - boundary_positions[idx]
        for idx in range(len(boundary_positions) - 1)
    ]

    return {
        "labels": labels,
        "boundary_positions": boundary_positions,
        "boundary_spacings": boundary_spacings,
        "bottom_margin": 0,
        "top_margin": top_margin,
        "segment_lengths": segment_lengths,
        "start_label": start_label % 3,
        "start_direction": normalize_direction(start_direction),
        "next_label": current_label,
        "next_direction": current_direction,
    }


def hcp_positions_from_labels(labels: list[int]) -> list[int]:
    count = len(labels)
    return [
        index
        for index in range(count)
        if labels[(index - 1) % count] == labels[(index + 1) % count]
    ]


def interface_hcp_gap_layers(
    upper_labels: list[int],
    lower_labels: list[int],
) -> int | None:
    combined = upper_labels + lower_labels
    interface_layer = len(upper_labels)
    hcp_positions = hcp_positions_from_labels(combined)
    previous = [value for value in hcp_positions if value < interface_layer]
    following = [value for value in hcp_positions if value >= interface_layer]
    if not previous or not following:
        return None
    return following[0] - previous[-1]


def select_ttg_ab_interface_phase(
    upper_spec: BlockSpec,
    lower_spec: BlockSpec,
    a0: float,
) -> tuple[int, int]:
    if upper_spec.twin_code is None or lower_spec.twin_code is None:
        return lower_spec.start_label, lower_spec.start_direction

    d111 = a0 / math.sqrt(3.0)

    upper_total_layers = choose_compatible_total_layers(
        max(2, int(round(nm_to_angstrom(upper_spec.thickness_nm) / d111))),
        TWIN_N111[upper_spec.twin_code],
    )
    upper_plan = plan_nanotwin_stack(
        total_layers=upper_total_layers,
        layers_per_twin=TWIN_N111[upper_spec.twin_code],
        start_label=upper_spec.start_label,
        start_direction=upper_spec.start_direction,
    )
    upper_labels = list(upper_plan["labels"])

    lower_total_layers = choose_compatible_total_layers(
        max(2, int(round(nm_to_angstrom(lower_spec.thickness_nm) / d111))),
        TWIN_N111[lower_spec.twin_code],
    )

    allowed_gaps = {
        TWIN_N111[upper_spec.twin_code],
        TWIN_N111[lower_spec.twin_code],
    }
    candidates: list[tuple[tuple[int, int, int, int], int, int]] = []

    for start_label in (0, 1, 2):
        for start_direction in (1, -1):
            lower_plan = plan_nanotwin_stack(
                total_layers=lower_total_layers,
                layers_per_twin=TWIN_N111[lower_spec.twin_code],
                start_label=start_label,
                start_direction=start_direction,
            )
            gap = interface_hcp_gap_layers(upper_labels, list(lower_plan["labels"]))
            if gap is None:
                continue
            if gap not in allowed_gaps:
                continue
            score = (
                0 if gap == TWIN_N111[lower_spec.twin_code] else 1,
                abs(gap - TWIN_N111[lower_spec.twin_code]),
                start_label,
                0 if start_direction == -1 else 1,
            )
            candidates.append((score, start_label, start_direction))

    if not candidates:
        return lower_spec.start_label, lower_spec.start_direction

    candidates.sort(key=lambda item: item[0])
    _, start_label, start_direction = candidates[0]
    return start_label, start_direction


def run(cmd: list[str], cwd: Path) -> None:
    print(">>", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=str(cwd))


def ensure_atomsk(atomsk_value: str) -> str:
    if atomsk_value == "auto":
        resolved = shutil.which("atomsk")
    else:
        resolved = atomsk_value
    if not resolved or not Path(resolved).exists():
        raise FileNotFoundError(
            "Atomsk executable was not found. Use --atomsk to provide a valid path."
        )
    return str(Path(resolved).resolve())


def get_honeycomb_angle(row: int, col: int) -> float:
    if row % 2 == 0:
        return 0.0 if col % 2 == 0 else 30.0
    return 90.0 if col % 2 == 0 else 60.0


def generate_honeycomb_nodes(
    nx_cols: int,
    ny_rows: int,
    lx: float,
    ly: float,
    lz: float,
    x_shift_factor: float,
    y_shift_factor: float,
    packing_scale: float,
    pattern_rotation_deg: float = 0.0,
    grain_angle_offset_deg: float = 0.0,
) -> list[tuple[float, float, float, float]]:
    if nx_cols < 1 or ny_rows < 1:
        raise ValueError("nx_cols and ny_rows must be positive")

    dx = lx / nx_cols
    dy = ly / ny_rows
    x_shift = x_shift_factor * dx
    y_shift = y_shift_factor * dy
    cx = 0.5 * lx
    cy = 0.5 * ly

    nodes: list[tuple[float, float, float, float]] = []
    for row in range(ny_rows):
        y = row * dy + y_shift
        x_offset = 0.5 * dx if row % 2 == 1 else 0.0
        for col in range(nx_cols):
            x = col * dx + x_offset + x_shift
            x = (cx + packing_scale * (x - cx)) % lx
            y_scaled = (cy + packing_scale * (y - cy)) % ly
            x_rot, y_rot = rotate_xy(x, y_scaled, cx, cy, pattern_rotation_deg)
            angle = (get_honeycomb_angle(row, col) + grain_angle_offset_deg) % 360.0
            nodes.append((x_rot % lx, y_rot % ly, 0.5 * lz, angle))
    return nodes


def write_poly_file(poly_file: Path, lx: float, ly: float, lz: float, nodes: list[tuple[float, float, float, float]]) -> None:
    lines = [f"box {lx:.6f} {ly:.6f} {lz:.6f}", ""]
    lines.append("# Honeycomb orientation map")
    lines.append("# even rows: 0 30 0 30 ...")
    lines.append("# odd rows: 90 60 90 60 ...")
    lines.append("")
    for x, y, z, angle in nodes:
        lines.append(f"node {x:.6f} {y:.6f} {z:.6f}   0   0   {angle:.6f}")
    lines.append("")
    poly_file.write_text("\n".join(lines), encoding="utf-8")


def write_fcc_seed(output_path: Path, atomsk: str, a0: float, element: str) -> None:
    run(
        [
            atomsk,
            "--create",
            "fcc",
            f"{a0}",
            element,
            "orient",
            "[100]",
            "[010]",
            "[001]",
            output_path.name,
        ],
        cwd=output_path.parent,
    )


def write_nanotwin_seed(
    output_path: Path,
    a0: float,
    element: str,
    total_thickness_angstrom: float,
    twin_code: str,
    nx_inplane: int,
    ny_inplane: int,
    start_label: int = 0,
    start_direction: int = 1,
) -> dict[str, object]:
    if twin_code not in TWIN_N111:
        raise KeyError(f"Unsupported twin code: {twin_code}")

    if nx_inplane < 1 or ny_inplane < 1:
        raise ValueError("nx_inplane and ny_inplane must be positive")

    d111 = a0 / math.sqrt(3.0)
    a2d = a0 / math.sqrt(2.0)
    layers_per_twin = TWIN_N111[twin_code]
    requested_layers = max(2, int(round(total_thickness_angstrom / d111)))
    total_layers = choose_compatible_total_layers(requested_layers, layers_per_twin)
    actual_thickness = total_layers * d111
    stack_plan = plan_nanotwin_stack(
        total_layers=total_layers,
        layers_per_twin=layers_per_twin,
        start_label=start_label,
        start_direction=start_direction,
    )

    v1 = (a2d, 0.0, 0.0)
    v2 = (0.5 * a2d, 0.5 * math.sqrt(3.0) * a2d, 0.0)
    v3 = (0.0, 0.0, actual_thickness)
    offsets = [
        (0.0, 0.0, 0.0),
        ((v1[0] + v2[0]) / 3.0, (v1[1] + v2[1]) / 3.0, 0.0),
        (2.0 * (v1[0] + v2[0]) / 3.0, 2.0 * (v1[1] + v2[1]) / 3.0, 0.0),
    ]

    atoms: list[tuple[str, float, float, float]] = []

    for layer, offset_index in enumerate(stack_plan["labels"]):
        ox, oy, _ = offsets[offset_index]
        z = layer * d111
        for i in range(nx_inplane):
            for j in range(ny_inplane):
                x = i * v1[0] + j * v2[0] + ox
                y = i * v1[1] + j * v2[1] + oy
                atoms.append((element, x, y, z))

    cell_a = (nx_inplane * v1[0], nx_inplane * v1[1], 0.0)
    cell_b = (ny_inplane * v2[0], ny_inplane * v2[1], 0.0)
    cell_c = v3

    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("CRYSTAL\n")
        handle.write("PRIMVEC\n")
        handle.write(f"{cell_a[0]:.12f} {cell_a[1]:.12f} {cell_a[2]:.12f}\n")
        handle.write(f"{cell_b[0]:.12f} {cell_b[1]:.12f} {cell_b[2]:.12f}\n")
        handle.write(f"{cell_c[0]:.12f} {cell_c[1]:.12f} {cell_c[2]:.12f}\n")
        handle.write("PRIMCOORD\n")
        handle.write(f"{len(atoms)} 1\n")
        for symbol, x, y, z in atoms:
            handle.write(f"{symbol} {x:.12f} {y:.12f} {z:.12f}\n")

    return {
        "d111_angstrom": d111,
        "layers_per_twin": layers_per_twin,
        "requested_layers": requested_layers,
        "total_layers": total_layers,
        "requested_thickness_angstrom": total_thickness_angstrom,
        "actual_thickness_angstrom": actual_thickness,
        "bottom_margin_layers": stack_plan["bottom_margin"],
        "top_margin_layers": stack_plan["top_margin"],
        "boundary_positions_layers": stack_plan["boundary_positions"],
        "boundary_spacings_layers": stack_plan["boundary_spacings"],
        "start_label": stack_plan["start_label"],
        "start_direction": stack_plan["start_direction"],
        "next_label": stack_plan["next_label"],
        "next_direction": stack_plan["next_direction"],
    }


def write_orthogonal_single_twin_seed(
    output_path: Path,
    a0: float,
    element: str,
    total_thickness_angstrom: float,
    twin_code: str,
    box_xy_nm: float,
    start_label: int = 0,
    start_direction: int = 1,
) -> dict[str, object]:
    if twin_code not in TWIN_N111:
        raise KeyError(f"Unsupported twin code: {twin_code}")

    d111 = a0 / math.sqrt(3.0)
    x_cell = a0 * math.sqrt(3.0 / 2.0)
    y_cell = a0 / math.sqrt(2.0)
    requested_layers = max(2, int(round(total_thickness_angstrom / d111)))
    total_layers = choose_compatible_total_layers(requested_layers, TWIN_N111[twin_code])
    actual_thickness = total_layers * d111
    nx, ny = choose_orthogonal_single_twin_repeats(a0, box_xy_nm)
    stack_plan = plan_nanotwin_stack(
        total_layers=total_layers,
        layers_per_twin=TWIN_N111[twin_code],
        start_label=start_label,
        start_direction=start_direction,
    )

    layer_offsets = {
        0: [(0.0, 0.0), (0.5, 0.5)],
        1: [(1.0 / 3.0, 0.0), (5.0 / 6.0, 0.5)],
        2: [(2.0 / 3.0, 0.0), (1.0 / 6.0, 0.5)],
    }

    atoms: list[tuple[str, float, float, float]] = []
    for layer, offset_index in enumerate(stack_plan["labels"]):
        z = layer * d111
        for i in range(nx):
            for j in range(ny):
                for fx, fy in layer_offsets[offset_index]:
                    x = (i + fx) * x_cell
                    y = (j + fy) * y_cell
                    atoms.append((element, x, y, z))

    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("CRYSTAL\n")
        handle.write("PRIMVEC\n")
        handle.write(f"{nx * x_cell:.12f} 0.000000000000 0.000000000000\n")
        handle.write(f"0.000000000000 {ny * y_cell:.12f} 0.000000000000\n")
        handle.write(f"0.000000000000 0.000000000000 {actual_thickness:.12f}\n")
        handle.write("PRIMCOORD\n")
        handle.write(f"{len(atoms)} 1\n")
        for sym, x, y, z in atoms:
            handle.write(f"{sym} {x:.12f} {y:.12f} {z:.12f}\n")

    return {
        "d111_angstrom": d111,
        "layers_per_twin": TWIN_N111[twin_code],
        "requested_layers": requested_layers,
        "total_layers": total_layers,
        "requested_thickness_angstrom": total_thickness_angstrom,
        "actual_thickness_angstrom": actual_thickness,
        "bottom_margin_layers": stack_plan["bottom_margin"],
        "top_margin_layers": stack_plan["top_margin"],
        "boundary_positions_layers": stack_plan["boundary_positions"],
        "boundary_spacings_layers": stack_plan["boundary_spacings"],
        "start_label": stack_plan["start_label"],
        "start_direction": stack_plan["start_direction"],
        "next_label": stack_plan["next_label"],
        "next_direction": stack_plan["next_direction"],
        "orthogonal_single_twin": True,
        "nx_repeat": nx,
        "ny_repeat": ny,
        "x_cell_angstrom": x_cell,
        "y_cell_angstrom": y_cell,
    }


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def parse_lammps_header(lmp_path: Path) -> dict[str, object]:
    header: dict[str, object] = {}
    with lmp_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for _ in range(20):
            line = handle.readline()
            if not line:
                break
            stripped = line.strip()
            if stripped.endswith("atoms"):
                header["atoms"] = int(stripped.split()[0])
            elif stripped.endswith("xlo xhi"):
                parts = stripped.split()
                header["xlo_xhi"] = [float(parts[0]), float(parts[1])]
            elif stripped.endswith("ylo yhi"):
                parts = stripped.split()
                header["ylo_yhi"] = [float(parts[0]), float(parts[1])]
            elif stripped.endswith("zlo zhi"):
                parts = stripped.split()
                header["zlo_zhi"] = [float(parts[0]), float(parts[1])]
    return header


def relative_to_root(path: Path, root_dir: Path) -> str:
    return str(path.resolve().relative_to(root_dir.resolve()))


def should_use_direct_block_build(spec: BlockSpec) -> bool:
    return spec.grain_code == "single"


def should_copy_single_block_model(block_count: int) -> bool:
    return block_count == 1


def should_unskew_for_lammps_export(spec: BlockSpec) -> bool:
    return should_use_direct_block_build(spec) and spec.twin_code is not None


def read_json_if_exists(path: Path) -> dict[str, object] | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def normalize_rel_path(text: str) -> str:
    return text.replace("/", "\\")


def block_cache_is_current(block_dir: Path, spec: BlockSpec) -> bool:
    final_lmp = block_dir / "wrapped.lmp"
    if not final_lmp.exists():
        return False

    if not (should_use_direct_block_build(spec) and spec.twin_code is not None):
        return True

    metadata = read_json_if_exists(block_dir / "metadata.json")
    if not metadata:
        return False

    seed_info = metadata.get("seed")
    if not isinstance(seed_info, dict):
        return False

    return (
        seed_info.get("mode") == "orthogonal_nanotwin"
        and seed_info.get("orthogonal_single_twin") is True
        and seed_info.get("start_label") == (spec.start_label % 3)
        and seed_info.get("start_direction") == normalize_direction(spec.start_direction)
    )


def model_output_is_current(
    model_dir: Path,
    effective_specs: list[BlockSpec],
    family: str,
    phase: str,
    variant: str,
    root_dir: Path,
) -> bool:
    final_lmp = model_dir / "wrapped.lmp"
    if not final_lmp.exists():
        return False

    metadata = read_json_if_exists(model_dir / "metadata.json")
    if not metadata:
        return False

    if metadata.get("family") != family:
        return False
    if metadata.get("phase") != phase:
        return False
    if metadata.get("variant") != variant:
        return False

    expected_layers = [asdict(spec) for spec in effective_specs]
    if metadata.get("layers") != expected_layers:
        return False

    expected_block_dirs = [
        normalize_rel_path(str(Path("cache") / "blocks" / spec.block_id()))
        for spec in effective_specs
    ]
    actual_block_dirs = [
        normalize_rel_path(str(item))
        for item in metadata.get("block_dirs", [])
    ]
    if actual_block_dirs != expected_block_dirs:
        return False

    if should_copy_single_block_model(len(effective_specs)) != bool(metadata.get("single_block_copy", False)):
        return False

    for spec in effective_specs:
        block_dir = root_dir / "cache" / "blocks" / spec.block_id()
        if not block_cache_is_current(block_dir, spec):
            return False

    return True


def build_block(spec: BlockSpec, context: BuildContext, force: bool = False) -> Path:
    if spec.grain_code not in GRAIN_TEMPLATES:
        raise KeyError(f"Unsupported grain template: {spec.grain_code}")
    if spec.twin_code is not None and spec.twin_code not in TWIN_N111:
        raise KeyError(f"Unsupported twin code: {spec.twin_code}")

    root_dir = context.root_dir.resolve()
    block_dir = root_dir / "cache" / "blocks" / spec.block_id()
    final_lmp = block_dir / "wrapped.lmp"
    if block_cache_is_current(block_dir, spec) and not force:
        return block_dir

    if block_dir.exists() and force:
        shutil.rmtree(block_dir)
    elif block_dir.exists():
        shutil.rmtree(block_dir)
    block_dir.mkdir(parents=True, exist_ok=True)

    lx = nm_to_angstrom(spec.box_xy_nm)
    ly = nm_to_angstrom(spec.box_xy_nm)
    lz = nm_to_angstrom(spec.thickness_nm)

    seed_path = block_dir / "seed.xsf"
    if should_use_direct_block_build(spec) and spec.twin_code is not None:
        seed_info = write_orthogonal_single_twin_seed(
            output_path=seed_path,
            a0=context.a0,
            element=context.element,
            total_thickness_angstrom=lz,
            twin_code=spec.twin_code,
            box_xy_nm=spec.box_xy_nm,
            start_label=spec.start_label,
            start_direction=spec.start_direction,
        )
        seed_info["mode"] = "orthogonal_nanotwin"
    elif spec.twin_code is None:
        write_fcc_seed(seed_path, context.atomsk, context.a0, context.element)
        seed_info = {
            "mode": "fcc",
            "actual_thickness_angstrom": lz,
        }
    else:
        seed_info = write_nanotwin_seed(
            output_path=seed_path,
            a0=context.a0,
            element=context.element,
            total_thickness_angstrom=lz,
            twin_code=spec.twin_code,
            nx_inplane=context.twin_nx,
            ny_inplane=context.twin_ny,
            start_label=spec.start_label,
            start_direction=spec.start_direction,
        )
        seed_info["mode"] = "nanotwin"

    if should_use_direct_block_build(spec):
        # Pure twin / pure single-crystal blocks should be written directly from the seed,
        # avoiding the 1-grain Voronoi route that can introduce fake boundaries.
        run([context.atomsk, "seed.xsf", "-wrap", "wrapped.cfg"], cwd=block_dir)
    else:
        nodes = generate_honeycomb_nodes(
            nx_cols=GRAIN_TEMPLATES[spec.grain_code]["nx"],
            ny_rows=GRAIN_TEMPLATES[spec.grain_code]["ny"],
            lx=lx,
            ly=ly,
            lz=lz,
            x_shift_factor=spec.shiftx,
            y_shift_factor=spec.shifty,
            packing_scale=spec.packing_scale,
            pattern_rotation_deg=spec.pattern_rotation_deg,
            grain_angle_offset_deg=spec.grain_angle_offset_deg,
        )
        write_poly_file(block_dir / "poly_hex_regular.txt", lx, ly, lz, nodes)

        run(
            [
                context.atomsk,
                "--polycrystal",
                "seed.xsf",
                "poly_hex_regular.txt",
                "poly.cfg",
            ],
            cwd=block_dir,
        )
        run([context.atomsk, "poly.cfg", "-wrap", "wrapped.cfg"], cwd=block_dir)

    if should_unskew_for_lammps_export(spec):
        run([context.atomsk, "wrapped.cfg", "-unskew", "lammps"], cwd=block_dir)
    else:
        run([context.atomsk, "wrapped.cfg", "lammps"], cwd=block_dir)

    metadata = {
        "block_spec": asdict(spec),
        "grain_template": GRAIN_TEMPLATES[spec.grain_code],
        "seed": seed_info,
        "box_angstrom": {"lx": lx, "ly": ly, "lz": lz},
        "outputs": {
            "seed": "seed.xsf",
            "poly": "poly_hex_regular.txt",
            "cfg": "wrapped.cfg",
            "lmp": "wrapped.lmp",
        },
        "header": parse_lammps_header(final_lmp),
    }
    write_json(block_dir / "metadata.json", metadata)
    return block_dir


def attach_block_edge_phases(block_specs: list[BlockSpec]) -> list[BlockSpec]:
    adjusted_specs: list[BlockSpec] = []
    for spec in block_specs:
        adjusted_specs.append(
            BlockSpec(
            grain_code=spec.grain_code,
            twin_code=spec.twin_code,
            thickness_nm=spec.thickness_nm,
            box_xy_nm=spec.box_xy_nm,
            packing_scale=spec.packing_scale,
            shiftx=spec.shiftx,
            shifty=spec.shifty,
            pattern_rotation_deg=spec.pattern_rotation_deg,
            grain_angle_offset_deg=spec.grain_angle_offset_deg,
            start_label=spec.start_label % 3,
            start_direction=normalize_direction(spec.start_direction),
        )
        )

    return adjusted_specs


def attach_family_specific_phase_rules(
    block_specs: list[BlockSpec],
    family: str,
    a0: float,
) -> list[BlockSpec]:
    adjusted = attach_block_edge_phases(block_specs)
    if family != "TTG":
        return adjusted

    if len(adjusted) <= 1:
        return adjusted

    selected_specs: list[BlockSpec] = [
        BlockSpec(
            grain_code=adjusted[0].grain_code,
            twin_code=adjusted[0].twin_code,
            thickness_nm=adjusted[0].thickness_nm,
            box_xy_nm=adjusted[0].box_xy_nm,
            packing_scale=adjusted[0].packing_scale,
            shiftx=adjusted[0].shiftx,
            shifty=adjusted[0].shifty,
            pattern_rotation_deg=adjusted[0].pattern_rotation_deg,
            grain_angle_offset_deg=adjusted[0].grain_angle_offset_deg,
            start_label=0,
            start_direction=1,
        )
    ]

    for spec in adjusted[1:]:
        previous = selected_specs[-1]
        selected_label, selected_direction = select_ttg_ab_interface_phase(previous, spec, a0)
        selected_specs.append(
            BlockSpec(
                grain_code=spec.grain_code,
                twin_code=spec.twin_code,
                thickness_nm=spec.thickness_nm,
                box_xy_nm=spec.box_xy_nm,
                packing_scale=spec.packing_scale,
                shiftx=spec.shiftx,
                shifty=spec.shifty,
                pattern_rotation_deg=spec.pattern_rotation_deg,
                grain_angle_offset_deg=spec.grain_angle_offset_deg,
                start_label=selected_label,
                start_direction=selected_direction,
            )
        )

    return selected_specs


def prepare_merge_paths(block_dirs: list[Path], model_dir: Path) -> tuple[list[str], str, str]:
    input_cfgs = [
        os.path.relpath(str((block_dir / "wrapped.cfg").resolve()), start=str(model_dir.resolve()))
        for block_dir in block_dirs
    ]
    return input_cfgs, "merged.cfg", "wrapped.cfg"


def build_model_from_blocks(
    model_name: str,
    family: str,
    phase: str,
    variant: str,
    block_specs: list[BlockSpec],
    context: BuildContext,
    force: bool = False,
) -> Path:
    root_dir = context.root_dir.resolve()
    model_dir = root_dir / "models" / family / model_name / variant
    final_lmp = model_dir / "wrapped.lmp"
    effective_specs = attach_family_specific_phase_rules(block_specs, family, context.a0)

    if (
        final_lmp.exists()
        and not force
        and model_output_is_current(
            model_dir=model_dir,
            effective_specs=effective_specs,
            family=family,
            phase=phase,
            variant=variant,
            root_dir=root_dir,
        )
    ):
        return model_dir

    if model_dir.exists():
        shutil.rmtree(model_dir)
    model_dir.mkdir(parents=True, exist_ok=True)

    block_dirs = [build_block(spec, context, force=force) for spec in effective_specs]

    if should_copy_single_block_model(len(block_dirs)):
        source_block = block_dirs[0]
        shutil.copy2(source_block / "wrapped.cfg", model_dir / "wrapped.cfg")
        shutil.copy2(source_block / "wrapped.lmp", model_dir / "wrapped.lmp")

        metadata = {
            "model_name": model_name,
            "family": family,
            "phase": phase,
            "variant": variant,
            "layers": [asdict(spec) for spec in effective_specs],
            "block_dirs": [relative_to_root(path, root_dir) for path in block_dirs],
            "outputs": {
                "cfg": "wrapped.cfg",
                "lmp": "wrapped.lmp",
            },
            "header": parse_lammps_header(model_dir / "wrapped.lmp"),
            "single_block_copy": True,
        }
        write_json(model_dir / "metadata.json", metadata)
        return model_dir

    block_cfgs, output_cfg, output_wrapped = prepare_merge_paths(block_dirs, model_dir)

    merge_cmd = [
        context.atomsk,
        "--merge",
        "Z",
        str(len(block_cfgs)),
        *block_cfgs,
        output_cfg,
    ]
    run(merge_cmd, cwd=model_dir)
    run([context.atomsk, output_cfg, "-wrap", output_wrapped], cwd=model_dir)
    run([context.atomsk, output_wrapped, "lammps"], cwd=model_dir)

    metadata = {
        "model_name": model_name,
        "family": family,
        "phase": phase,
        "variant": variant,
        "layers": [asdict(spec) for spec in effective_specs],
        "block_dirs": [relative_to_root(path, root_dir) for path in block_dirs],
        "outputs": {
            "cfg": "wrapped.cfg",
            "lmp": "wrapped.lmp",
        },
        "header": parse_lammps_header(final_lmp),
    }
    write_json(model_dir / "metadata.json", metadata)
    return model_dir


def gradient_variant_specs(model: dict[str, str], variant: str, context: BuildContext) -> list[BlockSpec]:
    family = model.get("family", "")

    if variant == "AB":
        layout = [
            (model["grain_a"], model["twin_a"], LAYER_THICKNESS_NM["AB_A"], "A"),
            (model["grain_b"], model["twin_b"], LAYER_THICKNESS_NM["AB_B"], "B"),
        ]
    elif variant == "ABAB":
        layout = [
            (model["grain_a"], model["twin_a"], LAYER_THICKNESS_NM["ABAB_A"], "A"),
            (model["grain_b"], model["twin_b"], LAYER_THICKNESS_NM["ABAB_B"], "B"),
            (model["grain_a"], model["twin_a"], LAYER_THICKNESS_NM["ABAB_A"], "A"),
            (model["grain_b"], model["twin_b"], LAYER_THICKNESS_NM["ABAB_B"], "B"),
        ]
    else:
        raise ValueError(f"Unsupported variant: {variant}")

    normalized_layout = []
    for grain_code, twin_code, thickness_nm, layer_role in layout:
        if family == "TTG":
            grain_code = "single"
        elif family == "GSG":
            twin_code = None

        if family == "TTG":
            pattern_rotation_deg = 0.0
            grain_angle_offset_deg = 0.0
        else:
            pattern_rotation_deg = 90.0 if layer_role == "B" else 0.0
            grain_angle_offset_deg = 90.0 if layer_role == "B" else 0.0

        normalized_layout.append(
            (
                grain_code,
                twin_code,
                thickness_nm,
                pattern_rotation_deg,
                grain_angle_offset_deg,
            )
        )

    return [
        BlockSpec(
            grain_code=grain_code,
            twin_code=twin_code,
            thickness_nm=thickness_nm,
            box_xy_nm=context.box_xy_nm,
            packing_scale=context.packing_scale,
            shiftx=context.shiftx,
            shifty=context.shifty,
            pattern_rotation_deg=pattern_rotation_deg,
            grain_angle_offset_deg=grain_angle_offset_deg,
        )
        for grain_code, twin_code, thickness_nm, pattern_rotation_deg, grain_angle_offset_deg in normalized_layout
    ]


def collect_pure_components(models: Iterable[dict[str, str]], context: BuildContext) -> list[tuple[str, str | None]]:
    seen: set[tuple[str, str | None]] = set()
    ordered: list[tuple[str, str | None]] = []
    for model in models:
        specs = gradient_variant_specs(model, "AB", context)
        for pair in ((specs[0].grain_code, specs[0].twin_code), (specs[1].grain_code, specs[1].twin_code)):
            if pair not in seen:
                seen.add(pair)
                ordered.append(pair)
    return ordered


def classify_component(grain_code: str, twin_code: str | None) -> tuple[str, str, str]:
    if grain_code == "single" and twin_code is not None:
        return "PURE_TWIN", f"single_{twin_code}", f"PURE_TWIN__single__{twin_code}"
    if twin_code is None:
        return "PURE_GRAIN", f"{grain_code}_FCC", f"PURE_GRAIN__{grain_code}__FCC"
    twin_tag = twin_code or "FCC"
    return "COMPONENT", f"{grain_code}_{twin_tag}", f"COMPONENT__{grain_code}__{twin_tag}"


def build_plan_models(
    context: BuildContext,
    phase: str,
    include_pure: bool,
    families: set[str] | None = None,
    model_names: set[str] | None = None,
    variants: set[str] | None = None,
    force: bool = False,
) -> dict[str, object]:
    selected_models = [
        model
        for model in PLAN_MODELS
        if phase == "all" or model["phase"] == phase
    ]
    if families:
        selected_models = [
            model for model in selected_models if model["family"] in families
        ]
    if model_names:
        selected_models = [
            model for model in selected_models if model["name"] in model_names
        ]

    built_gradient: list[dict[str, str]] = []
    for model in selected_models:
        for variant in ("AB", "ABAB"):
            if variants and variant not in variants:
                continue
            build_model_from_blocks(
                model_name=model["name"],
                family=model["family"],
                phase=model["phase"],
                variant=variant,
                block_specs=gradient_variant_specs(model, variant, context),
                context=context,
                force=force,
            )
            built_gradient.append(
                {
                    "name": model["name"],
                    "family": model["family"],
                    "phase": model["phase"],
                    "variant": variant,
                }
            )

    built_pure: list[dict[str, str]] = []
    if include_pure:
        pure_components = collect_pure_components(selected_models if phase != "P1" else PLAN_MODELS, context)
        for grain_code, twin_code in pure_components:
            pure_family, model_name, pure_job_name = classify_component(grain_code, twin_code)
            block_specs = [
                BlockSpec(
                    grain_code=grain_code,
                    twin_code=twin_code,
                    thickness_nm=LAYER_THICKNESS_NM["PURE"],
                    box_xy_nm=context.box_xy_nm,
                    packing_scale=context.packing_scale,
                    shiftx=context.shiftx,
                    shifty=context.shifty,
                )
            ]
            build_model_from_blocks(
                model_name=model_name,
                family=pure_family,
                phase="P1",
                variant="single",
                block_specs=block_specs,
                context=context,
                force=force,
            )
            built_pure.append({"name": pure_job_name, "variant": "single", "phase": "P1", "family": pure_family})

    manifest = {
        "phase": phase,
        "include_pure": include_pure,
        "families": sorted(families) if families else [],
        "model_names": sorted(model_names) if model_names else [],
        "variants": sorted(variants) if variants else [],
        "gradient_models": built_gradient,
        "pure_models": built_pure,
    }
    manifests_dir = context.root_dir / "manifests"
    manifests_dir.mkdir(parents=True, exist_ok=True)
    write_json(manifests_dir / f"build_manifest_{phase}.json", manifest)
    return manifest


def parse_layers(text: str, context: BuildContext) -> list[BlockSpec]:
    specs: list[BlockSpec] = []
    for chunk in [item.strip() for item in text.split(",") if item.strip()]:
        parts = [item.strip() for item in chunk.split(":")]
        if len(parts) != 3:
            raise argparse.ArgumentTypeError(
                "Each layer must use grain:twin:thickness_nm or grain:none:thickness_nm"
            )
        grain_code, twin_code, thickness_text = parts
        twin_value = None if twin_code.lower() == "none" else twin_code
        try:
            thickness_nm = float(thickness_text)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(
                f"Invalid thickness value: {thickness_text!r}"
            ) from exc
        specs.append(
            BlockSpec(
                grain_code=grain_code,
                twin_code=twin_value,
                thickness_nm=thickness_nm,
                box_xy_nm=context.box_xy_nm,
                packing_scale=context.packing_scale,
                shiftx=context.shiftx,
                shifty=context.shifty,
            )
        )
    if not specs:
        raise argparse.ArgumentTypeError("At least one layer must be provided")
    return specs


DEFAULT_CONFIG = Path(__file__).with_name("config.toml")


def _strip_toml_comment(line: str) -> str:
    """删除字符串外的 # 注释，满足本项目配置文件的 TOML 写法。"""
    quote = None
    for index, char in enumerate(line):
        if char in {"'", '\"'}:
            if quote is None:
                quote = char
            elif quote == char:
                quote = None
        elif char == "#" and quote is None:
            return line[:index]
    return line


def _load_simple_toml(text: str) -> dict:
    """无需第三方包的简化 TOML 读取器，覆盖本项目 config.toml 使用的语法。"""
    result: dict = {}
    current: dict | None = result
    for raw_line in text.splitlines():
        line = _strip_toml_comment(raw_line).strip()
        if not line:
            continue
        if line.startswith("[[") and line.endswith("]]" ):
            section = line[2:-2].strip()
            result.setdefault(section, []).append({})
            current = result[section][-1]
            continue
        if line.startswith("[") and line.endswith("]"):
            current = result
            for part in line[1:-1].strip().split("."):
                current = current.setdefault(part, {})
            continue
        if "=" not in line or current is None:
            raise ValueError(f"无法读取配置行：{raw_line}")
        key, raw_value = (item.strip() for item in line.split("=", 1))
        normalized = raw_value.replace("true", "True").replace("false", "False")
        try:
            current[key] = ast.literal_eval(normalized)
        except (SyntaxError, ValueError) as exc:
            raise ValueError(f"配置值格式错误：{raw_line}") from exc
    return result


def load_configuration(config_path: Path) -> dict:
    """读取 TOML 配置，并在缺少关键段落时给出清晰错误。"""
    if not config_path.is_file():
        raise FileNotFoundError(f"配置文件不存在：{config_path}")
    if tomllib is not None:
        with config_path.open("rb") as handle:
            config = tomllib.load(handle)
    else:
        config = _load_simple_toml(config_path.read_text(encoding="utf-8"))
    required = {"paths", "build", "layer_thickness_nm", "twin_n111", "grain_templates", "plan_models", "run"}
    missing = required.difference(config)
    if missing:
        raise ValueError(f"配置文件缺少必要段落：{', '.join(sorted(missing))}")
    return config


def apply_configuration(config: dict, config_path: Path) -> None:
    """把配置写入建模器运行时参数；模型函数随后统一使用这些值。"""
    global DEFAULT_ATOMSK, DEFAULT_OUTPUT_ROOT, LAYER_THICKNESS_NM, TWIN_N111, GRAIN_TEMPLATES, PLAN_MODELS
    paths = config["paths"]
    build = config["build"]
    DEFAULT_ATOMSK = str(paths["atomsk"])
    output_root = Path(str(paths["output_root"]))
    if not output_root.is_absolute():
        output_root = (config_path.parent / output_root).resolve()
    DEFAULT_OUTPUT_ROOT = str(output_root)
    LAYER_THICKNESS_NM = {key: float(value) for key, value in config["layer_thickness_nm"].items()}
    TWIN_N111 = {key: int(value) for key, value in config["twin_n111"].items()}
    GRAIN_TEMPLATES = {
        key: {"nx": int(value["nx"]), "ny": int(value["ny"]), "label": str(value["label"])}
        for key, value in config["grain_templates"].items()
    }
    PLAN_MODELS = [dict(item) for item in config["plan_models"]]
    # build 段由 build_parser 读取；在这里保留以便验证所有公共参数均存在。
    required_build = {"box_xy_nm", "a0_angstrom", "element", "packing_scale", "shiftx", "shifty", "twin_nx", "twin_ny"}
    missing = required_build.difference(build)
    if missing:
        raise ValueError(f"[build] 缺少参数：{', '.join(sorted(missing))}")


def run_arguments_from_config(config: dict) -> list[str]:
    """在不写命令行参数时，使用 [run] 段决定本次构建任务。"""
    run = config["run"]
    command = str(run.get("command", "plan"))
    common = [command]
    if command == "block":
        common.extend(["--grain", str(run["block_grain"]), "--twin", str(run["block_twin"]), "--thickness", str(run["block_thickness_nm"])])
    elif command == "merge":
        layers = run["merge_layers"]
        if isinstance(layers, list):
            layers = ",".join(str(item) for item in layers)
        common.extend(["--name", str(run["merge_name"]), "--family", str(run["merge_family"]), "--phase", str(run["merge_phase"]), "--layers", str(layers)])
    elif command == "plan":
        common.extend(["--phase", str(run.get("plan_phase", "all"))])
        for field, flag in (("families", "--families"), ("models", "--models"), ("variants", "--variants")):
            values = run.get(field, [])
            if values:
                common.extend([flag, ",".join(str(item) for item in values)])
        if not bool(run.get("include_pure", True)):
            common.append("--skip-pure")
    else:
        raise ValueError("[run].command 只能是 block、merge 或 plan")
    if bool(run.get("force", False)):
        common.insert(0, "--force")
    return common


def configuration_summary(config: dict, config_path: Path) -> str:
    build = config["build"]
    run = config["run"]
    return (
        f"已读取配置：{config_path}\n"
        f"a0 = {build['a0_angstrom']} angstrom；盒子面内尺寸 = {build['box_xy_nm']} nm；"
        f"本次任务 = {run.get('command', 'plan')}"
    )

CONFIG_BUILD: dict = {
    "box_xy_nm": 18.0,
    "a0_angstrom": 3.633,
    "element": "Cu",
    "packing_scale": 1.0,
    "shiftx": 0.25,
    "shifty": 0.35,
    "twin_nx": 8,
    "twin_ny": 8,
}

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            "Unified builder for grain-only, nanotwin polycrystal, custom AB/ABAB stacking, "
            "and the v5 double-gradient model matrix."
        )
    )
    parser.add_argument("--atomsk", default=DEFAULT_ATOMSK, help="Atomsk 可执行文件路径；填 auto 时从 PATH 查找")
    parser.add_argument("--root", default=DEFAULT_OUTPUT_ROOT, help="输出根目录")
    parser.add_argument("--box-xy", type=float, default=CONFIG_BUILD["box_xy_nm"], help="面内盒子尺寸（nm）")
    parser.add_argument("--a0", type=float, default=CONFIG_BUILD["a0_angstrom"], help="晶格常数（Å）")
    parser.add_argument("--element", default=CONFIG_BUILD["element"], help="元素符号")
    parser.add_argument("--packing-scale", type=float, default=CONFIG_BUILD["packing_scale"], help="种子打包系数")
    parser.add_argument("--shiftx", type=float, default=CONFIG_BUILD["shiftx"], help="种子 x 方向偏移系数")
    parser.add_argument("--shifty", type=float, default=CONFIG_BUILD["shifty"], help="种子 y 方向偏移系数")
    parser.add_argument("--twin-nx", type=int, default=CONFIG_BUILD["twin_nx"], help="孪晶种子 x 向重复数")
    parser.add_argument("--twin-ny", type=int, default=CONFIG_BUILD["twin_ny"], help="孪晶种子 y 向重复数")
    parser.add_argument("--force", action="store_true", help="Rebuild even if outputs already exist")

    subparsers = parser.add_subparsers(dest="command", required=True)

    block_parser = subparsers.add_parser(
        "block",
        allow_abbrev=False,
        help="Build one single-layer block. Use twin=none for grain-only output.",
    )
    block_parser.add_argument("--grain", required=True, choices=sorted(GRAIN_TEMPLATES.keys()))
    block_parser.add_argument("--twin", required=True, help="Twin code (T4/T6/.../T30) or none")
    block_parser.add_argument("--thickness", required=True, type=float, help="Layer thickness in nm")

    merge_parser = subparsers.add_parser(
        "merge",
        allow_abbrev=False,
        help="Build a custom merged model from a comma-separated layer list.",
    )
    merge_parser.add_argument("--name", required=True, help="Model name")
    merge_parser.add_argument("--family", default="CUSTOM", help="Model family")
    merge_parser.add_argument("--phase", default="custom", help="Phase label")
    merge_parser.add_argument(
        "--layers",
        required=True,
        help="Comma-separated layers: grain:twin:thickness_nm",
    )

    plan_parser = subparsers.add_parser(
        "plan",
        allow_abbrev=False,
        help="Build models from 双层梯度纳米孪晶铜_建模计划_v5.md presets.",
    )
    plan_parser.add_argument(
        "--phase",
        choices=["all", "P0", "P1", "P2"],
        default="all",
        help="Which build phase to generate",
    )
    plan_parser.add_argument(
        "--skip-pure",
        action="store_true",
        help="Skip Pure baseline models",
    )
    plan_parser.add_argument(
        "--families",
        default="",
        help="Optional comma-separated family filter, e.g. TTG,GSG,DG,REF,PURE",
    )
    plan_parser.add_argument(
        "--models",
        default="",
        help="Optional comma-separated model names, e.g. TTG-2,DG-A,REF",
    )
    plan_parser.add_argument(
        "--variants",
        default="",
        help="Optional comma-separated variant filter, e.g. AB or ABAB",
    )

    return parser


def build_context_from_args(args: argparse.Namespace) -> BuildContext:
    return BuildContext(
        atomsk=ensure_atomsk(args.atomsk),
        root_dir=Path(args.root).resolve(),
        element=args.element,
        a0=args.a0,
        box_xy_nm=args.box_xy,
        packing_scale=args.packing_scale,
        shiftx=args.shiftx,
        shifty=args.shifty,
        twin_nx=args.twin_nx,
        twin_ny=args.twin_ny,
    )


def parse_csv_filter(text: str) -> set[str] | None:
    values = {item.strip() for item in text.split(",") if item.strip()}
    return values or None


def write_plan_reference(context: BuildContext) -> None:
    manifests_dir = context.root_dir / "manifests"
    manifests_dir.mkdir(parents=True, exist_ok=True)
    plan_reference = {
        "box_xy_nm": context.box_xy_nm,
        "layer_thickness_nm": LAYER_THICKNESS_NM,
        "twin_codes": {
            code: {
                "n111": n111,
                "spacing_nm": round(n111 * context.a0 / math.sqrt(3.0) / 10.0, 4),
            }
            for code, n111 in TWIN_N111.items()
        },
        "grain_templates": GRAIN_TEMPLATES,
        "plan_models": PLAN_MODELS,
    }
    write_json(manifests_dir / "plan_reference.json", plan_reference)


def main() -> int:
    bootstrap = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    bootstrap.add_argument("--config", default=str(DEFAULT_CONFIG), help="TOML 配置文件路径")
    bootstrap.add_argument("--show-config", action="store_true", help="仅显示配置摘要，不建模")
    bootstrap_args, remaining = bootstrap.parse_known_args()
    config_path = Path(bootstrap_args.config).expanduser().resolve()
    config = load_configuration(config_path)
    apply_configuration(config, config_path)
    global CONFIG_BUILD
    CONFIG_BUILD = config["build"]
    if bootstrap_args.show_config:
        print(configuration_summary(config, config_path))
        return 0
    parser = build_parser()
    # 不输入子命令时，完全按 [run] 段执行；命令行参数仍可覆盖本次任务选择。
    args = parser.parse_args(remaining or run_arguments_from_config(config))
    print(configuration_summary(config, config_path))
    context = build_context_from_args(args)
    context.root_dir.mkdir(parents=True, exist_ok=True)
    write_plan_reference(context)

    if args.command == "block":
        twin_code = None if args.twin.lower() == "none" else args.twin
        spec = BlockSpec(
            grain_code=args.grain,
            twin_code=twin_code,
            thickness_nm=args.thickness,
            box_xy_nm=context.box_xy_nm,
            packing_scale=context.packing_scale,
            shiftx=context.shiftx,
            shifty=context.shifty,
        )
        block_dir = build_block(spec, context, force=args.force)
        print(f"Built block: {block_dir}")
        return 0

    if args.command == "merge":
        block_specs = parse_layers(args.layers, context)
        model_dir = build_model_from_blocks(
            model_name=args.name,
            family=args.family,
            phase=args.phase,
            variant="custom",
            block_specs=block_specs,
            context=context,
            force=args.force,
        )
        print(f"Built merged model: {model_dir}")
        return 0

    if args.command == "plan":
        manifest = build_plan_models(
            context=context,
            phase=args.phase,
            include_pure=not args.skip_pure,
            families=parse_csv_filter(args.families),
            model_names=parse_csv_filter(args.models),
            variants=parse_csv_filter(args.variants),
            force=args.force,
        )
        print(json.dumps(manifest, indent=2, ensure_ascii=False))
        return 0

    parser.error("Unsupported command")
    return 2


if __name__ == "__main__":
    sys.exit(main())

