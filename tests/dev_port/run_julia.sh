#!/bin/bash
# Run one dev-port test case with the patched Julia reference (branch port-ref).
#
#   tests/dev_port/run_julia.sh <case> [NTHREADS=8]
#
# Work directory: tests/dev_port/work/julia_t<NTHREADS>/<case>/ (wiped first), outputs in its io/, log in run.log,
# step times in timings_julia.txt. Environment: JULIA (binary), NUJULIET_REF (patched reference worktree).
set -euo pipefail
CASE=${1:?usage: run_julia.sh <case> [NTHREADS]}
NT=${2:-8}
HERE=$(cd "$(dirname "$0")" && pwd)
JULIA=${JULIA:-$HOME/.julia/juliaup/julia-1.10.11+0.x64.linux.gnu/bin/julia}
export NUJULIET_REF=${NUJULIET_REF:-/scratch/local.ssd/knapp/job_25503048.pbs-m1/NuJuliet_dev_ref}
[ -f "$HERE/julia/$CASE.jl" ] || { echo "unknown case $CASE"; exit 2; }
(( NT <= 32 )) || { echo "at most 32 threads"; exit 2; }

WORK=$HERE/work/julia_t$NT/$CASE
rm -rf "$WORK"
mkdir -p "$WORK/io"
cd "$WORK"
export OPENBLAS_NUM_THREADS=$NT
echo "julia case $CASE, $NT threads, work dir $WORK, ref $NUJULIET_REF ($(git -C "$NUJULIET_REF" rev-parse --short HEAD))"
/usr/bin/time -v "$JULIA" --project=@v1.10 -t "$NT,1" "$HERE/julia/$CASE.jl" > run.log 2> time.log || {
    echo "FAILED: see $WORK/run.log and $WORK/time.log"; tail -n 30 time.log; exit 1; }
grep -E "Elapsed|Maximum resident" time.log
grep -E "TOTAL_TIME|CASE_DONE" run.log
