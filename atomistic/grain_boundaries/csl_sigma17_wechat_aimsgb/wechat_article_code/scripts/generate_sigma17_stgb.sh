#!/usr/bin/env bash
# Build a Sigma17{223}[110] FCC Cu symmetric tilt grain boundary with Atomsk.
# Grain 2 is written as Al only to retain a distinct LAMMPS atom type; the
# LAMMPS input maps both types to the same Cu EAM potential.

set -euo pipefail

ATOMSK_BIN="${ATOMSK_BIN:-atomsk}"
LATTICE_CONSTANT="${LATTICE_CONSTANT:-3.597}"
DUPLICATE_X="${DUPLICATE_X:-20}"
DUPLICATE_Y="${DUPLICATE_Y:-20}"
DUPLICATE_Z="${DUPLICATE_Z:-1}"
OUTPUT="${OUTPUT:-sigma17_stgb_initial.lmp}"

"${ATOMSK_BIN}" --create fcc "${LATTICE_CONSTANT}" Cu orient [-33-4] [2-2-3] [110] \
  -duplicate "${DUPLICATE_X}" "${DUPLICATE_Y}" "${DUPLICATE_Z}" grain_1.xsf
"${ATOMSK_BIN}" --create fcc "${LATTICE_CONSTANT}" Al orient [-334] [-22-3] [110] \
  -duplicate "${DUPLICATE_X}" "${DUPLICATE_Y}" "${DUPLICATE_Z}" grain_2.xsf
"${ATOMSK_BIN}" --merge Y 2 grain_1.xsf grain_2.xsf "${OUTPUT}"

echo "Created ${OUTPUT}"
