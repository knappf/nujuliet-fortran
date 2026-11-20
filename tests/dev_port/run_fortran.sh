#!/bin/bash
# Run one dev-port test case with the Fortran port (namelists in tests/dev_port/fortran/<case>/, run in file order).
#
#   tests/dev_port/run_fortran.sh <case> [NTHREADS=8]
#
# Work directory: tests/dev_port/work/fortran_t<NTHREADS>/<case>/ (wiped first), outputs in its io/, log in run.log.
# Environment: NUJULIET (binary, default ./nujuliet of the repository).
set -euo pipefail
CASE=${1:?usage: run_fortran.sh <case> [NTHREADS]}
NT=${2:-8}
HERE=$(cd "$(dirname "$0")" && pwd)
NUJULIET=${NUJULIET:-$HERE/../../nujuliet}
[ -d "$HERE/fortran/$CASE" ] || { echo "unknown case $CASE"; exit 2; }
(( NT <= 32 )) || { echo "at most 32 threads"; exit 2; }

WORK=$HERE/work/fortran_t$NT/$CASE
rm -rf "$WORK"
mkdir -p "$WORK/io/parse"
cd "$WORK"
cp "$HERE/fortran/$CASE"/*.nml .
export OMP_NUM_THREADS=$NT MKL_NUM_THREADS=$NT OMP_STACKSIZE=512M
ulimit -s unlimited
echo "fortran case $CASE, $NT threads, work dir $WORK, binary $NUJULIET"
/usr/bin/time -v "$NUJULIET" $(ls *.nml | sort) > run.log 2> time.log || {
    echo "FAILED: see $WORK/run.log and $WORK/time.log"; tail -n 30 run.log; exit 1; }
grep -E "Elapsed|Maximum resident" time.log
