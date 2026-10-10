from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# The user asked to compare the 5 strained VASP structures.
# strain_00 is the undeformed reference and is not plotted here.
VASP_POINTS = np.array(
    [
        [0.05, 5.139227],
        [0.10, 11.690390],
        [0.15, 15.658513],
        [0.20, 11.730121],
        [0.25, 9.785866],
    ],
    dtype=float,
)


def read_curve(path: Path) -> np.ndarray:
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if len(parts) < 4:
            continue
        rows.append([float(value) for value in parts[:4]])
    return np.asarray(rows, dtype=float)


def nearest_stress(curve: np.ndarray, strain: float) -> tuple[float, float]:
    index = int(np.argmin(np.abs(curve[:, 0] - strain)))
    return float(curve[index, 0]), float(curve[index, 2])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("base_dir", type=Path, help="Directory containing EAM/SUS2 data folders")
    args = parser.parse_args()

    base_dir = args.base_dir.resolve()
    output_dir = base_dir / "对比"
    output_dir.mkdir(exist_ok=True)

    curve_files = {
        "EAM": base_dir / "EAM" / "results_Y" / "1_stress_strain.dat",
        "SUS2": base_dir / "SUS2" / "results_Y" / "1_stress_strain.dat",
    }

    curves = {name: read_curve(path) for name, path in curve_files.items()}

    plt.figure(figsize=(10, 6))
    colors = {"EAM": "#1f77b4", "SUS2": "#ff7f0e"}
    for name, curve in curves.items():
        plt.plot(curve[:, 0], curve[:, 2], label=name, linewidth=1.8, alpha=0.9, color=colors[name])

    plt.scatter(
        VASP_POINTS[:, 0],
        VASP_POINTS[:, 1],
        label="VASP (5 points)",
        color="black",
        marker="D",
        s=58,
        zorder=5,
    )

    for strain, stress in VASP_POINTS:
        plt.annotate(
            f"{strain:.2f}",
            (strain, stress),
            textcoords="offset points",
            xytext=(5, 6),
            fontsize=8,
            color="black",
        )

    plt.title("Stress-Strain Comparison: Direction Y, Structure 1")
    plt.xlabel("Strain")
    plt.ylabel("Stress (Direction Y) [GPa]")
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.legend()
    plt.tight_layout()

    image_path = output_dir / "comparison_Y_1_with_vasp_points.png"
    plt.savefig(image_path, dpi=300)
    plt.close()

    summary_rows: list[dict[str, float | str]] = []
    for strain, vasp_stress in VASP_POINTS:
        row: dict[str, float | str] = {
            "strain": float(strain),
            "vasp_stress_gpa": float(vasp_stress),
        }
        for name, curve in curves.items():
            nearest_strain, stress = nearest_stress(curve, float(strain))
            row[f"{name.lower()}_nearest_strain"] = nearest_strain
            row[f"{name.lower()}_stress_gpa"] = stress
            row[f"{name.lower()}_error_gpa"] = stress - float(vasp_stress)
            row[f"{name.lower()}_abs_error_gpa"] = abs(stress - float(vasp_stress))
        summary_rows.append(row)

    csv_path = output_dir / "comparison_Y_1_with_vasp_points.csv"
    with csv_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    eam_errors = np.array([float(row["eam_error_gpa"]) for row in summary_rows], dtype=float)
    sus2_errors = np.array([float(row["sus2_error_gpa"]) for row in summary_rows], dtype=float)
    summary_text = "\n".join(
        [
            "VASP points were taken from work_Y/strain_05,10,15,20,25.",
            "VASP stress uses -sigma_yy from OUTCAR 'in kB', converted to GPa by dividing by 10.",
            "",
            f"EAM MAE  = {np.mean(np.abs(eam_errors)):.6f} GPa",
            f"EAM RMSE = {np.sqrt(np.mean(eam_errors ** 2)):.6f} GPa",
            f"SUS2 MAE  = {np.mean(np.abs(sus2_errors)):.6f} GPa",
            f"SUS2 RMSE = {np.sqrt(np.mean(sus2_errors ** 2)):.6f} GPa",
        ]
    )
    accuracy_path = output_dir / "comparison_Y_1_accuracy.txt"
    accuracy_path.write_text(summary_text, encoding="utf-8")

    print(f"Saved plot: {image_path}")
    print(f"Saved table: {csv_path}")
    print(f"Saved summary: {accuracy_path}")
    print(summary_text)


if __name__ == "__main__":
    main()
