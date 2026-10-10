#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess


def run_cmd(cmd):
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def build_periodic_nanotwin(
    atomsk_cmd="atomsk",
    a=3.615,
    nx=3,
    ny_twin=2,
    nz=3,
    n_domains=4,
    output_prefix="Cu_periodic_nanotwin",
    output_format="cfg",
    cleanup=False
):
    """
    Build a periodic uniformly spaced nanotwin structure in Cu using Atomsk.

    Orientation:
      X = [11-2]
      Y = [111]   # twin-plane normal
      Z = [-110]

    Structure:
      A / B / A / B / ...
    where B is the mirrored variant of A.

    Parameters
    ----------
    nx, ny_twin, nz : int
        Duplications of one twin lamella block.
    n_domains : int
        Number of lamellae stacked along Y.
        Use an even number for symmetric A/B alternation.
    """

    if n_domains < 2:
        raise ValueError("n_domains must be >= 2")
    if ny_twin < 1:
        raise ValueError("ny_twin must be >= 1")

    base_file = f"{output_prefix}_A.xsf"
    twin_file = f"{output_prefix}_B.xsf"
    final_file = f"{output_prefix}.{output_format}"

    # 1) Build one Cu lamella A
    # Atomsk create + orient + duplicate are standard operations
    run_cmd([
        atomsk_cmd,
        "--create", "fcc", str(a), "Cu",
        "orient", "[11-2]", "[111]", "[-110]",
        "-duplicate", str(nx), str(ny_twin), str(nz),
        base_file
    ])

    # 2) Mirror A to obtain twin-related variant B
    # Official twin tutorial uses -mirror 0 Y -wrap
    run_cmd([
        atomsk_cmd,
        base_file,
        "-mirror", "0", "Y",
        "-wrap",
        twin_file
    ])

    # 3) Stack A/B/A/B/... periodically along Y
    merge_cmd = [atomsk_cmd, "--merge", "Y", str(n_domains)]

    for i in range(n_domains):
        merge_cmd.append(base_file if i % 2 == 0 else twin_file)

    merge_cmd.append(final_file)
    run_cmd(merge_cmd)

    print("\nDone.")
    print(f"Output: {final_file}")
    print(f"Lamella thickness control (duplicate along Y): {ny_twin}")
    print(f"Number of lamellae: {n_domains}")
    print("Structure sequence:", " / ".join("A" if i % 2 == 0 else "B" for i in range(n_domains)))

    if cleanup:
        for f in [base_file, twin_file]:
            if os.path.exists(f):
                os.remove(f)
        print("Temporary files removed.")


if __name__ == "__main__":
    # =====================================
    # User parameters: edit here directly
    # =====================================
    atomsk_cmd = "atomsk"     # or full path to atomsk

    a = 3.615                 # Cu lattice constant, Angstrom

    # Small test system (try to keep under ~1000 atoms first)
    nx = 3                    # in-plane size along X
    ny_twin = 2               # thickness of one twin lamella block
    nz = 3                    # in-plane size along Z

    n_domains = 4             # A/B/A/B -> 4 lamellae, periodic nanotwin
    output_prefix = "Cu_periodic_nanotwin_small"
    output_format = "cfg"     # cfg / xsf / lmp etc. if Atomsk supports it
    cleanup = False

    build_periodic_nanotwin(
        atomsk_cmd=atomsk_cmd,
        a=a,
        nx=nx,
        ny_twin=ny_twin,
        nz=nz,
        n_domains=n_domains,
        output_prefix=output_prefix,
        output_format=output_format,
        cleanup=cleanup
    )
