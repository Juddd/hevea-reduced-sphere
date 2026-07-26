#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
out="$(mktemp -d "${TMPDIR:-/tmp}/reduced-sphere-fast-XXXXXX")"
log="$out/native.log"
(cd "$out" && HEVEA_PREVIEW_NX=300 HEVEA_PREVIEW_NY=600 OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}" "$root/hevea_reduced_sphere" 4000 4000 21 142 997 .52 .5 .237) | tee "$log"
test "$(grep -c '^METRIC Stage=' "$log")" -eq 4
! grep -Eq 'MinPrimitive=-|MinFlowJacobian=-|nan|inf' "$log"
test "$(find "$out" -name '*.vtk' | wc -l)" -eq 8
printf 'fast_e2e=pass output=%s vtk=8 stages=4\n' "$out"
