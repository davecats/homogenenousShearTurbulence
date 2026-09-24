#!/bin/bash
# Run the test programs of a build.   usage:  tests/run_tests.sh [build-cpu|build-gpu] [nranks]
# Each test prints PASSED or FAILED; the script exits non-zero if any failed.
set -u
here=$(cd "$(dirname "$0")/.." && pwd)
build=${1:-build-cpu}
np=${2:-1}
work=$(mktemp -d)
cd "$work"
status=0
run() {   # run <ranks> <exe> <deck>
  echo "--- $2 on $1 rank(s), deck $3"
  mpirun -np "$1" "$here/$build/$2" "$here/tests/decks/$3" 2>&1 | grep -v -E "Authorization|hcoll|^ *$"
  [ "${PIPESTATUS[0]}" -eq 0 ] || status=1
}
run "$np" test_roundtrip small.in
run "$np" test_linsolve small.in
run 1     test_kelvin kelvin.in
rm -rf "$work"
[ $status -eq 0 ] && echo "ALL PASSED" || echo "SOME FAILED"
exit $status
