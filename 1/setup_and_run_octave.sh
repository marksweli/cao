#!/usr/bin/env bash
# Prepare redbKIT for GNU Octave (Linux) and run the Re = 100 cylinder case.
# In MATLAB this script is NOT needed: download redbKIT, run `make` in its
# root folder, set REDBKIT_PATH and run main_cylinder_Re100.m.
#
#   usage:  bash setup_and_run_octave.sh [redbKIT_dir]
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
RK="${1:-$HERE/../../redbKIT}"

command -v octave >/dev/null || sudo apt-get install -y octave octave-dev gnuplot-nox xvfb

if [ ! -d "$RK" ]; then
  git clone https://github.com/redbKIT/redbKIT.git "$RK"
  git -C "$RK" checkout 59b8ba99469e7cdf1030aeeaf85da82f16764526
fi
RK="$(cd "$RK" && pwd)"

# redbKIT's C assemblers include MATLAB's "blas.h" but do not use BLAS:
# an empty header is enough to build them with Octave's mkoctfile.
INC="$(mktemp -d)"; echo "/* dummy */" > "$INC/blas.h"
build() { (cd "$RK/$1" && mkoctfile --mex -I"$INC" -O3 -fopenmp -lgomp "$2" "${@:3}"); }
build FEM_library/Models/ADR ADR_assembler_C_omp.c
build FEM_library/Models/ADR Mass_assembler_C_omp.c
build FEM_library/Models/ADR ADR_SUPGassembler_C_omp.c ../../Core/Tools.c
build FEM_library/Models/CFD CFD_assembler_C_omp.c ../../Core/Tools.c
build FEM_library/Models/CFD CFD_assembler_ExtForces.c

cd "$HERE"
# solver (about 1.5-2 h on 2 cores; log in results/run_log.txt); figures are then regenerated with the
# Qt toolkit through a virtual X server for better quality
REDBKIT_PATH="$RK" octave --no-gui --eval "main_cylinder_Re100"
xvfb-run -a octave --no-gui --eval "graphics_toolkit qt; postprocess_cylinder_Re100"
