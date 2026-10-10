#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess


def run_cmd(cmd):
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def build_uniform_twin_region(atomsk_cmd, a, nx, nz, thickness, n_pairs, prefix):
    """
    Build one uniform nanotwin region:
    A(thickness) / B(thickness) repeated n_pairs times
    """
    a_file = f"{prefix}_A.xsf"
    b_file = f"{prefix}_B.xsf"
    out_file = f"{prefix}.xsf"

    # build A block
    run_cmd([
        atomsk_cmd,
        "--create", "fcc", str(a), "Cu",
        "orient", "[11-2]", "[111]", "[-110]",
        "-duplicate", str(nx), str(thickness), str(nz),
        a_file
    ])

    # build B block by mirroring A
    run_cmd([
        atomsk_cmd,
        a_file,
        "-mirror", "0", "Y",
        "-wrap",
        b_file
    ])

    merge_cmd = [atomsk_cmd, "--merge", "Y", str(2 * n_pairs)]
    for i in range(n_pairs):
        merge_cmd.append(a_file)
        merge_cmd.append(b_file)
    merge_cmd.append(out_file)

    run_cmd(merge_cmd)
    return a_file, b_file, out_file


def build_two_region_gradient(
    atomsk_cmd="atomsk",
    a=3.615,
    nx=2,
    nz=2,
    thick_A=4,
    n_pairs_A=2,
    thick_B=6,
    n_pairs_B=3,
    output_file="Cu_AB_gradient.cfg",
    cleanup=False
):
    """
    Region A: uniform twin thickness = thick_A
    Region B: uniform twin thickness = thick_B
    Final structure = region A + region B
    """

    # build region A
    files_A = build_uniform_twin_region(
        atomsk_cmd=atomsk_cmd,
        a=a,
        nx=nx,
        nz=nz,
        thickness=thick_A,
        n_pairs=n_pairs_A,
        prefix="region_A"
    )

    # build region B
    files_B = build_uniform_twin_region(
        atomsk_cmd=atomsk_cmd,
        a=a,
        nx=nx,
        nz=nz,
        thickness=thick_B,
        n_pairs=n_pairs_B,
        prefix="region_B"
    )

    regionA_file = files_A[2]
    regionB_file = files_B[2]

    # merge A region + B region
    run_cmd([
        atomsk_cmd,
        "--merge", "Y", "2",
        regionA_file, regionB_file,
        output_file
    ])

    print("\nDone.")
    print(f"Output: {output_file}")
    print(f"Region A: thickness = {thick_A}, pairs = {n_pairs_A}")
    print(f"Region B: thickness = {thick_B}, pairs = {n_pairs_B}")

    if cleanup:
        for f in list(files_A) + list(files_B):
            if os.path.exists(f):
                os.remove(f)
        print("Temporary files removed.")


if __name__ == "__main__":
    # =====================================
    # only edit these parameters
    # =====================================
    atomsk_cmd = "atomsk"
    a = 3.615

    nx = 2
    nz = 2

    # region A: fine twins
    thick_A = 4
    n_pairs_A = 2

    # region B: coarse twins
    thick_B = 6
    n_pairs_B = 3

    output_file = "Cu_AB_gradient.cfg"
    cleanup = False

    build_two_region_gradient(
        atomsk_cmd=atomsk_cmd,
        a=a,
        nx=nx,
        nz=nz,
        thick_A=thick_A,
        n_pairs_A=n_pairs_A,
        thick_B=thick_B,
        n_pairs_B=n_pairs_B,
        output_file=output_file,
        cleanup=cleanup
    )
