#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
out="${HEVEA_PAPER_OUTPUT_DIRECTORY:-$root/artifacts/paper-grid}"
mkdir -p "$out"
exe_hash="$(sha256sum "$root/hevea_reduced_sphere" | cut -d' ' -f1)"
ribbon_fractions="0.41148336662952684,0.83293586428607069,0.89850213833742598,0.96908999242130689"
lambda_powers="0.60979123180328654,0,0"
transition_splits="0.80193065649632189,0.5,0.5"
profile_hash="profile-20260725-author-envelope-paper-flow-certified-v5"
config="4000x20000-ridges=21,142,997-radius=.52-eta=.5-fraction=.237-ribbon=$ribbon_fractions-lambda=$lambda_powers-splits=$transition_splits-profile=$profile_hash-stage3-anchor=equator-exe=$exe_hash"
config_hash="$(printf '%s' "$config" | sha256sum | cut -d' ' -f1)"
if [[ -f "$out/manifest.json" ]] && grep -q "$config_hash" "$out/manifest.json" &&
   [[ -f "$out/reduced_sphere_stage=3_dir=2_ridges=997.bin" ]] &&
   [[ -f "$out/SHA256SUMS" ]] && (cd "$out" && sha256sum -c SHA256SUMS >/dev/null); then
  printf 'paper_grid=cache-hit config_hash=%s\n' "$config_hash"
  exit 0
fi
stamp="$(date +%Y%m%d-%H%M%S)"
[[ ! -f "$out/manifest.json" ]] || mv "$out/manifest.json" "$out/manifest.stale-$stamp.json"
[[ ! -f "$out/SHA256SUMS" ]] || mv "$out/SHA256SUMS" "$out/SHA256SUMS.stale-$stamp"
manifest_tmp="$out/manifest.json.tmp.$$"
trap 'rm -f "$manifest_tmp"' EXIT
cd "$out"
HEVEA_DIAGNOSTICS_ONLY=1 HEVEA_BINARY_OUTPUT=1 OMP_NUM_THREADS="${OMP_NUM_THREADS:-16}" /usr/bin/time -v -o resource.txt stdbuf -oL -eL "$root/hevea_reduced_sphere" 4000 20000 21 142 997 .52 .5 .237 | tee native.log
printf '{\n  "output_directory": "%s",\n  "grid": [4000, 20000],\n  "ridges": [21, 142, 997],\n  "ball_radius": 0.52,\n  "eta": 0.5,\n  "target_fraction": 0.237,\n  "ribbon_fractions": [%s],\n  "lambda_powers": [%s],\n  "transition_splits": [%s],\n  "profile_hash": "%s",\n  "stage3_anchor": "equator",\n  "executable_sha256": "%s",\n  "config_hash": "%s",\n  "native_log": "%s/native.log"\n}\n' "$out" "$ribbon_fractions" "$lambda_powers" "$transition_splits" "$profile_hash" "$exe_hash" "$config_hash" "$out" > "$manifest_tmp"
mv "$manifest_tmp" manifest.json
sha256sum ./*.bin native.log manifest.json > SHA256SUMS
trap - EXIT
printf 'paper_grid=pass manifest=%s/manifest.json config_hash=%s\n' "$out" "$config_hash"
