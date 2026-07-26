#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
paper="${HEVEA_PAPER_OUTPUT_DIRECTORY:-$root/artifacts/paper-grid}"
visual="$paper/figure9-mma-v5/figure9-visual-report.json"
audit_tmp="$(mktemp -d "${TMPDIR:-/tmp}/reduced-sphere-final-audit-XXXXXX")"
trap 'rm -rf -- "$audit_tmp"' EXIT

OMP_NUM_THREADS="${OMP_NUM_THREADS:-16}" "$root/tests/run_paper_grid.sh"
"$root/tools/verify_paper_metrics" --manifest "$paper/manifest.json"
"$root/tools/audit_mesh" --manifest "$paper/manifest.json"
(cd "$paper" && sha256sum -c SHA256SUMS)

test -f "$visual"
jq -e '.ClosedSphere == true and .ConnectedComponents == 1 and
  .InternalHoles == 0 and .AspectRelativeError < 0.10 and
  .RidgesExpected == [21,142,997] and .RidgesCounted == [21,142,997] and
  (.EnvelopeComparison | contains("envelope_compare=pass")) and
  (.ImageSHA256 | length) == 8' "$visual" >/dev/null
jq -r '.ImageSHA256 | to_entries[] | "\(.value)  \(.key)"' "$visual" > "$audit_tmp/IMAGE_SHA256SUMS"
sha256sum -c "$audit_tmp/IMAGE_SHA256SUMS"

printf '{broken' > "$audit_tmp/corrupt-manifest.json"
if "$root/tools/export_paper_meshes" "$audit_tmp/corrupt-manifest.json" "$audit_tmp/export" \
    >"$audit_tmp/corrupt.out" 2>"$audit_tmp/corrupt.err"; then
  printf 'corrupt_manifest=unexpected-success\n' >&2
  exit 1
fi
printf 'corrupt_manifest=structured-failure\n'

if ps -eo args= | awk -v executable="$root/hevea_reduced_sphere" '$1==executable{found=1} END{exit !found}'; then
  printf 'residual_native_process=fail\n' >&2
  exit 1
fi

test -f "$root/README.md"
rg -q 'f_\{1,3\}' "$root/README.md"
rg -q 'f_∞' "$root/README.md"
rg -q '"OutputMode" -> "Preview"' "$root/wolfram/CPHeveaReducedSphere.wl"
if rg -n 'TO''DO|FIX''ME' "$root/hevea_reduced_sphere.cpp" \
    "$root/hevea_reduced_sphere_wolfram.cpp" "$root/src" "$root/tools" \
    "$root/tests" "$root/docs"; then
  printf 'cleanliness=fail\n' >&2
  exit 1
fi

printf 'residual_native_process=pass visual_hashes=8 docs=pass cleanliness=pass final_audit=pass\n'
