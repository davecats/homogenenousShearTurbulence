#!/bin/bash
# Regression against the reference fields in tests/reference/.
#   tests/regression.sh [build-cpu|build-gpu] [nranks] [--update]
# Runs each deck of DECKS for its 50 steps and compares the final field with
# the stored reference at 1e-10 relative (CPU and GPU builds differ by
# ~1e-14).  --update rewrites the references from this build: do that only
# on purpose, after a change that is meant to alter the numbers.
set -u
here=$(cd "$(dirname "$0")/.." && pwd)
build=${1:-build-cpu}; np=${2:-1}; update=${3:-}
DECKS="small small_s2x small_stokes"
status=0
for d in $DECKS; do
  work=$(mktemp -d); cp "$here/tests/decks/$d.in" "$work/hst.in"
  ( cd "$work" && mpirun -np "$np" "$here/$build/hst" hst.in > run.log 2>&1 ) || { echo "$d: run failed"; status=1; continue; }
  if [ "$update" = "--update" ]; then
    cp "$work/Dati.cart.out" "$here/tests/reference/$d.fld"; echo "$d: reference updated"
  else
    echo "--- $d"
    python3 "$here/tests/compare_fields.py" "$work/Dati.cart.out" "$here/tests/reference/$d.fld" 1e-10 | tail -2 || status=1
  fi
  rm -rf "$work"
done
[ $status -eq 0 ] && echo "REGRESSION OK" || echo "REGRESSION FAILED"
exit $status
