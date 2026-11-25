#!/bin/bash
set -euo pipefail

rm -rf build/CMakeCache.txt build/CMakeFiles build/external_build build/generated
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX"
cmake --build build -j "$CPU_COUNT"
cmake --install build
