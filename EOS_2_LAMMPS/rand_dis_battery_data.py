#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
from ase.io import write, iread
from ase.build import make_supercell


def main():
    if len(sys.argv) < 3:
        print("Usage: python rand_dis_battery_10_30_to_data.py <input.xyz> <out_dir>")
        print("Example: python rand_dis_battery_10_30_to_data.py CU331.xyz eos_data")
        sys.exit(1)

    in_files = sys.argv[1]
    out_dir = sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)

    # Identity supercell (same as your original)
    sc = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

    print(f"in_files: {in_files}  out_dir: {out_dir}")

    # EOS sampling by volume ratio V/V0
    vol_scale = [0.94,0.96,0.98,1.00,1.02,1.04,1.06,1.08]


    n_written = 0

    # If input has multiple frames, we’ll write files for each frame too.
    for frame_i, f in enumerate(iread(in_files)):
        aa = f
        a0 = make_supercell(aa, sc)

        spos = a0.get_scaled_positions()
        cell0 = a0.get_cell()

        for vs in vol_scale:
            ls = vs ** (1.0 / 3.0)  # linear scale so that volume scales by vs
            a = a0.copy()
            a.set_cell(cell0 * ls, scale_atoms=False)
            a.set_scaled_positions(spos)

            # File name
            # Example: data_V1.0200_frame0
            fname = os.path.join(out_dir, f"data_V{vs:.4f}_frame{frame_i}.lmp")

            # IMPORTANT:
            # - atom_style 'atomic' is what we want for EAM / box/relax / minimize
            # - specorder keeps element order stable (especially important if multi-species)
            # Here you are Cu, but leave it in for safety.
            write(
                fname,
                a,
                format="lammps-data",
                atom_style="atomic",
                specorder=["Cu"],
            )

            n_written += 1

    print("Wrote LAMMPS data files:", n_written)
    print("Example output:", os.path.join(out_dir, "data_V1.0000_frame0.lmp"))


if __name__ == "__main__":
    main()

