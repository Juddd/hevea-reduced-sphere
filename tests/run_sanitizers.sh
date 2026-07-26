#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
build="$(mktemp -d "${TMPDIR:-/tmp}/reduced-sphere-sanitize-build-XXXXXX")"
cmake -S "$root" -B "$build" -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS='-fsanitize=address,undefined -fno-omit-frame-pointer' -DCMAKE_EXE_LINKER_FLAGS='-fsanitize=address,undefined' >/dev/null
cmake --build "$build" -j 4 --target corrugation_tests numeric_kernel_tests >/dev/null
ASAN_OPTIONS=detect_leaks=1 "$build/corrugation_tests"
ASAN_OPTIONS=detect_leaks=1 "$build/numeric_kernel_tests"
printf 'asan_ubsan=pass\n'
