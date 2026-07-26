#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
exe="$root/hevea_reduced_sphere"
test -x "$exe"
tmpdir=$(mktemp -d /tmp/hevea-sphere-baseline.XXXXXX)
trap 'rm -rf -- "$tmpdir"' EXIT
set +e
(cd "$tmpdir" && "$exe" 160 800 9 30 90 0.5 0.9 0.08 > run.log 2>&1)
status=$?
set -e
log="$tmpdir/run.log"
grep -q 'ERROR: uncertified initial profile parameters' "$log"
printf 'expected_failure=pass exit=%s reason=uncertified_initial_profile\n' "$status"
test "$status" -eq 2
