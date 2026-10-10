#!/usr/bin/env bash
# Extract lines written as "@ i j energy mJ/m^2" by in.optimize_stgb.lmp.

set -euo pipefail

LOG_FILE="${1:-log.lammps}"
grep '^@' "${LOG_FILE}" | tr '@' ' ' | sort -n -k 3
