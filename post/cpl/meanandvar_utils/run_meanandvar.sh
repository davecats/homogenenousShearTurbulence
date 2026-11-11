#!/usr/bin/env bash


MEANANDVAR_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="/home/ws/th5982/research/projects/HST-S1S2/fortran-code/box_ar-effect/steady_state/data/Re-1000/B-4-3/A-0"

# Compute mean and uncertainty of Runtimedata
cpl run "$MEANANDVAR_DIR/meanandvar" < "$DATA_DIR/statistics/Runtimedata"
