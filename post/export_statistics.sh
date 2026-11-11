#!/usr/bin/env bash

# ============================================================================
# User configuration
# ============================================================================
# This wrapper can be run from any working directory with:
#
#     ./post/export_statistics.sh
#
# Change the variables below for the case that should be exported.

# Directory containing this repository.  This script is stored in ``post/``,
# so the repository root is its parent directory.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd -- "$SCRIPT_DIR/.." && pwd)"

# Directory containing hst.in and the statistics/ directory.
RUN_DIR="/home/ws/th5982/research/projects/HST-S1S2/fortran-code/parameter_studies/steady_state/productionrun/data/Re-1000/B-3-2/A-0_75/T-14"

# Output file and format.  Supported formats are "netcdf4" and "json".
OUTPUT_FILE="$RUN_DIR/statistics.nc"
OUTPUT_FORMAT="netcdf4"

# Initial transient time to exclude from the statistics.
# Set this to an empty string to preserve the default behavior (skip no rows).
T_INIT_STAT="40.0"

# Set to 1 if the flow relaminarized: all turbulent statistics are then stored as exactly zero.
RELAMINARIZED=0

# Virtual-environment settings.
# Set INITIALIZE_VIRTUAL_ENV=1 for the first run, or whenever the package set
# needs to be installed/updated.  Set it to 0 for normal runs afterwards.
USE_VIRTUAL_ENV=1
INITIALIZE_VIRTUAL_ENV=0
VIRTUAL_ENV_DIR="$PROJECT_DIR/.venv-post"
VIRTUAL_ENV_PYTHON="$VIRTUAL_ENV_DIR/bin/python"

# Packages required by the Python statistics exporter.  meanandvar itself is
# the separate CPL executable configured below and is not a pip package.
PYTHON_PACKAGES=(numpy xarray netCDF4)

# System Python used to create the virtual environment and, when
# USE_VIRTUAL_ENV=0, to run the exporter directly.
PYTHON_EXECUTABLE="python3"
MEANANDVAR_EXECUTABLE="$PROJECT_DIR/post/cpl/meanandvar_utils/meanandvar"

# ============================================================================
# Execution
# ============================================================================
set -euo pipefail

die() {
    echo "export_statistics.sh: $*" >&2
    exit 1
}

[[ -d "$PROJECT_DIR" ]] || die "project directory does not exist: $PROJECT_DIR"
[[ -d "$RUN_DIR" ]] || die "run directory does not exist: $RUN_DIR"
[[ -f "$RUN_DIR/hst.in" ]] || die "missing input deck: $RUN_DIR/hst.in"
[[ -f "$MEANANDVAR_EXECUTABLE" ]] || die "missing meanandvar executable: $MEANANDVAR_EXECUTABLE"

case "$USE_VIRTUAL_ENV" in
    0|1) ;;
    *) die "USE_VIRTUAL_ENV must be 0 or 1 (got: $USE_VIRTUAL_ENV)" ;;
esac
case "$INITIALIZE_VIRTUAL_ENV" in
    0|1) ;;
    *) die "INITIALIZE_VIRTUAL_ENV must be 0 or 1 (got: $INITIALIZE_VIRTUAL_ENV)" ;;
esac

if (( USE_VIRTUAL_ENV )); then
    if (( INITIALIZE_VIRTUAL_ENV )); then
        echo "Initializing virtual environment: $VIRTUAL_ENV_DIR"
        "$PYTHON_EXECUTABLE" -m venv "$VIRTUAL_ENV_DIR"
        "$VIRTUAL_ENV_PYTHON" -m pip install --upgrade pip
        "$VIRTUAL_ENV_PYTHON" -m pip install --upgrade "${PYTHON_PACKAGES[@]}"
    fi
    [[ -x "$VIRTUAL_ENV_PYTHON" ]] || die \
        "virtual environment not found: $VIRTUAL_ENV_DIR (set INITIALIZE_VIRTUAL_ENV=1 first)"
    PYTHON_EXECUTABLE="$VIRTUAL_ENV_PYTHON"
fi

case "$OUTPUT_FORMAT" in
    netcdf4|json) ;;
    *) die "OUTPUT_FORMAT must be netcdf4 or json (got: $OUTPUT_FORMAT)" ;;
esac

if [[ "$OUTPUT_FILE" != /* ]]; then
    OUTPUT_FILE="$PROJECT_DIR/$OUTPUT_FILE"
fi

ARGS=(
    -m post.python
    export_statistics
    "$RUN_DIR"
    --format "$OUTPUT_FORMAT"
    --output "$OUTPUT_FILE"
    --meanandvar "$MEANANDVAR_EXECUTABLE"
)

if [[ -n "$T_INIT_STAT" ]]; then
    ARGS+=(--t-init-stat "$T_INIT_STAT")
fi
if [[ "$RELAMINARIZED" == "1" ]]; then
    ARGS+=(--relaminarized)
fi

echo "Exporting statistics"
echo "  run directory: $RUN_DIR"
echo "  output:       $OUTPUT_FILE"
echo "  format:       $OUTPUT_FORMAT"
if [[ -n "$T_INIT_STAT" ]]; then
    echo "  t_init_stat:  $T_INIT_STAT"
else
    echo "  t_init_stat:  disabled (no rows skipped)"
fi

cd "$PROJECT_DIR"
PYTHONPATH="$PROJECT_DIR${PYTHONPATH:+:$PYTHONPATH}" \
    "$PYTHON_EXECUTABLE" "${ARGS[@]}"
