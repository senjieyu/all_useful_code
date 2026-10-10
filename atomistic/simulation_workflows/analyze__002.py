import sys

import matplotlib.pyplot as plt
import numpy as np


def analyze(data_file, plot_file, result_file, mode):
    data = np.loadtxt(data_file, comments="#")

    if data.ndim == 1:
        data = data.reshape(1, -1)

    strain = data[:, 0]
    stress = data[:, 1]

    plt.figure()
    plt.plot(strain, stress, linewidth=1.5, label=f"Shear {mode.upper()}")
    plt.xlabel("Shear Strain")
    plt.ylabel("Shear Stress (GPa)")
    plt.title(f"Shear Stress-Strain Curve ({mode.upper()})")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_file, dpi=200)
    plt.close()

    max_idx = int(np.argmax(np.abs(stress)))
    max_stress = float(stress[max_idx])
    yield_strain = float(strain[max_idx])

    with open(result_file, "a", encoding="utf-8") as handle:
        handle.write(f"Mode: {mode}\n")
        handle.write(f"Peak Shear Stress: {max_stress:.6f} GPa\n")
        handle.write(f"Shear Strain At Peak: {yield_strain:.6f}\n")
        handle.write("-" * 20 + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python analyze.py <data_file> <plot_file> <result_file> <mode>")
        sys.exit(1)

    analyze(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
