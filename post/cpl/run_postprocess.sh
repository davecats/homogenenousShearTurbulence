#!/usr/bin/env bash

# Post-process one HST case.  The only case-specific values normally needing
# changes are in this block (or DATA_DIR can be supplied as the first argument).
DATA_DIR="/home/ws/th5982/research/codes/homogenenousShearTurbulence/run"
# DATA_DIR=/home/ws/th5982/research/codes/homogenenousShearTurbulence/run-taylor_green
NF_MIN=1
NF_MAX=10
DN=1
NPROCS=8

# Set BUILD=0 when postprocess is already compiled and only the case changes.
BUILD=1

set -euo pipefail

# Directory of this script with symlinks resolved, so it can be run from any
# folder (e.g. ./path/to/run_postprocess.sh [DATA_DIR], or through a symlink).
POSTPROCESS_DIR="$(cd -- "$(dirname -- "$(realpath -- "${BASH_SOURCE[0]}")")" && pwd)"
DATA_DIR="$(realpath -- "${1:-$DATA_DIR}")"
POSTPROCESS_EXE="$POSTPROCESS_DIR/postprocess"

die() {
    echo "run_postprocess.sh: $*" >&2
    exit 1
}

[[ -f "$DATA_DIR/hst.in" ]] || die "missing $DATA_DIR/hst.in"
for required_dir in fields p_fields statistics; do
    [[ -d "$DATA_DIR/$required_dir" ]] || die "missing $DATA_DIR/$required_dir/"
done

if (( BUILD )); then
    (cd "$POSTPROCESS_DIR" && mpicpl make postprocess.cpl)
fi
[[ -x "$POSTPROCESS_EXE" ]] || die "post-processing executable not found: $POSTPROCESS_EXE"

# Keep generated inputs out of both the case and post/cpl/.  The CPL program
# accepts the post-processing input as argv[1] and the simulation input as
# argv[2] (readinput.cpl), so there is no need to copy either file manually.
WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/hst-postprocess.XXXXXX")"
trap 'rm -rf -- "$WORK_DIR"' EXIT

python3 "$POSTPROCESS_DIR/src/convert_hstin.py" \
    "$DATA_DIR/hst.in" "$WORK_DIR/readinput.in"

cat > "$WORK_DIR/postpro.in" <<EOF
nfmin=$NF_MIN
nfmax=$NF_MAX
dn=$DN

path_name=$DATA_DIR/
EOF

echo "Post-processing case: $DATA_DIR"
echo "Fields: $NF_MIN..$NF_MAX (dn=$DN), MPI ranks: $NPROCS"

cd "$POSTPROCESS_DIR"
mpirun -np "$NPROCS" "$POSTPROCESS_EXE" \
    "$WORK_DIR/postpro.in" "$WORK_DIR/readinput.in"
