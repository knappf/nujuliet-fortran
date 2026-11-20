#!/bin/bash
# Stage-4 unit check of the HF-ScRPA building blocks: Julia dumps (scrpa_blocks.jl) vs the Fortran test program
# (scrpa_blocks.f90, linked against the objects of BUILD). Usage: run_unit.sh [BUILD=build_s4] [NT=8]
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../../.." && pwd)
BUILD=${1:-build_s4}
NT=${2:-8}
JULIA=${JULIA:-$HOME/.julia/juliaup/julia-1.10.11+0.x64.linux.gnu/bin/julia}
WORK=$HERE/../work/unit_s4
REFIO=$HERE/../work/julia_t8/bundled/io/He4/HF
mkdir -p "$WORK/io/He4"
rm -rf "$WORK/io/He4/HF"; cp -r "$REFIO" "$WORK/io/He4/HF"
cd "$WORK"
if [ ! -f dump/energy.txt ] || [ "${REDO_JULIA:-0}" = 1 ]; then
    OPENBLAS_NUM_THREADS=$NT "$JULIA" --project=@v1.10 -t "$NT" "$HERE/scrpa_blocks.jl" > julia.log 2>&1 || { tail -30 julia.log; exit 1; }
fi
set +eu; source /etc/profile.d/30_meta_modules.sh >/dev/null 2>&1; module load intel-oneapi-hpc-compilers >/dev/null 2>&1; set -eu
OBJS=$(ls "$ROOT/$BUILD"/core/*.o "$ROOT/$BUILD"/hf/*.o "$ROOT/$BUILD"/app/nj_input.o "$ROOT/$BUILD"/rpa/*.o | grep -v nj_dispatch_rpa)
MKLLIB=$(dirname "$(which ifx)")/../../../mkl/latest/lib
FCLIB=$(dirname "$(which ifx)")/../lib
ifx -O2 -march=core-avx2 -qopenmp -I "$ROOT/$BUILD" "$HERE/scrpa_blocks.f90" $OBJS -o scrpa_blocks \
    -L"$MKLLIB" -Wl,-rpath,"$MKLLIB" -Wl,-rpath,"$FCLIB" -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm -ldl
ulimit -s unlimited; OMP_STACKSIZE=512M OMP_NUM_THREADS=$NT ./scrpa_blocks
