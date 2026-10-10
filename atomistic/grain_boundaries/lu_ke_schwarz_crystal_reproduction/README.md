# Lu Ke Schwarz-crystal reproduction

Reproducible Atomsk/LAMMPS inputs for a 16-grain copper polycrystal model used in the Lu Ke Schwarz-crystal reproduction folder.

## Contents

- `atomsk_modeling_commands.txt` creates eight oriented FCC-Cu components, generates the selected grains, merges them, and writes a LAMMPS structure.
- `polycrystal_nodes.txt` defines the 130.1 Å cubic box and the 16 Voronoi grain nodes.
- `melt.in` is the LAMMPS heating and constant-pressure input (1 K to 1401 K).
- `submit_lammps.slurm` is a portable SLURM submission template for the LAMMPS input.

## Run outline

1. Install [Atomsk](https://atomsk.univ-lille.fr/) and run the commands in `atomsk_modeling_commands.txt` from this directory.
2. Copy the resulting `final.lmp` into the LAMMPS run directory.
3. Supply a compatible copper EAM/FS potential as `MCu3.eam.fs` (not versioned here), then submit `submit_lammps.slurm`.

Generated XSF/CFG/LAMMPS structures, compressed trajectories, scheduler logs, potential files, and paper/document attachments are intentionally excluded from the source catalog.
