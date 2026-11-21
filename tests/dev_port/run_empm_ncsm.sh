#!/bin/bash
# EMPM / NCSM consumer checks (stage 5). EMPM and NCSM have no Julia counterpart: the reference is a LEGACY BASELINE made
# with the frozen stage-4 binary (legacy HF + EMPM/NCSM reading io/<Path>/Bin/*), which the dev-form chain (dev HF ->
# EMPM/NCSM reading io/<Path>/<Reference>/bin/*) must reproduce.
#
#   tests/dev_port/run_empm_ncsm.sh baseline [NT=16]   legacy inputs fortran/empm_ncsm_legacy/*.nml with NUJULIET_LEGACY
#                                                       -> work/empm_baseline/{io,log_NN.txt,times.txt}
#   tests/dev_port/run_empm_ncsm.sh dev [NT=16]        dev inputs fortran/empm_ncsm/*.nml with NUJULIET, then the error
#                                                       cases fortran/empm_ncsm/errors/*.nml (each must fail with the text
#                                                       of its '! EXPECT:' line; '! MKDIR: io/...' lines create stale
#                                                       directories first) -> work/empm_dev_t<NT>/
#   tests/dev_port/run_empm_ncsm.sh compare [NT=16]    python3 -I empm_compare.py work/empm_baseline/io work/empm_dev_t<NT>/io
#
# TAG=<suffix> puts the dev run into work/empm_dev<suffix>_t<NT> (e.g. for a second binary). Every namelist file runs in its own nujuliet invocation (wall time per file in times.txt). Environment: NUJULIET
# (default ./nujuliet of the repository), NUJULIET_LEGACY (default the frozen binary of commit d753f29, kept as
# runs/Pb208_13sh_hw12_scrpa/nujuliet_d753f29; to rebuild: git worktree add ../snap d753f29 && make there).
set -euo pipefail
MODE=${1:?usage: run_empm_ncsm.sh baseline|dev|compare [NTHREADS]}
NT=${2:-16}
HERE=$(cd "$(dirname "$0")" && pwd)
(( NT <= 32 )) || { echo "at most 32 threads"; exit 2; }
export OMP_NUM_THREADS=$NT MKL_NUM_THREADS=$NT OMP_STACKSIZE=512M
ulimit -s unlimited

run_dir() {   # binary, input dir, work dir
    local BIN=$1 IN=$2 WORK=$3
    rm -rf "$WORK"
    mkdir -p "$WORK/io"
    cd "$WORK"
    cp "$IN"/*.nml .
    : > times.txt
    echo "inputs $IN, $NT threads, work dir $WORK, binary $BIN"
    for f in $(ls *.nml | sort); do
        local tag=${f%.nml}
        if ! /usr/bin/time -f "$tag %e s %M KB" -o time_$tag.txt "$BIN" "$f" > log_$tag.txt 2>&1; then
            cat time_$tag.txt >> times.txt
            echo "FAILED: $f (see $WORK/log_$tag.txt)"; tail -n 20 log_$tag.txt; exit 1
        fi
        cat time_$tag.txt >> times.txt
        rm -f time_$tag.txt
        tail -n 1 times.txt
    done
}

case $MODE in
baseline)
    BIN=${NUJULIET_LEGACY:-$(cd "$(dirname "$0")/../.." && pwd)/runs/Pb208_13sh_hw12_scrpa/nujuliet_d753f29}
    run_dir "$BIN" "$HERE/fortran/empm_ncsm_legacy" "$HERE/work/empm_baseline"
    ;;
dev)
    BIN=${NUJULIET:-$HERE/../../nujuliet}
    WORK=$HERE/work/empm_dev${TAG:-}_t$NT
    run_dir "$BIN" "$HERE/fortran/empm_ncsm" "$WORK"
    # error paths: each file in a fresh io/ of its own, must exit non-zero (or print the error and skip) with the EXPECT text
    mkdir -p "$WORK/errors"
    for f in $(ls "$HERE/fortran/empm_ncsm/errors"/*.nml | sort); do
        tag=$(basename "$f" .nml)
        exp=$(sed -n 's/^! EXPECT: //p' "$f" | head -n 1)
        rm -rf "$WORK/errors/$tag"; mkdir -p "$WORK/errors/$tag/io"; cd "$WORK/errors/$tag"
        cp "$f" .
        # '! MKDIR: <dir below io/>' lines: directories left over from an earlier run (e.g. a stale io/<Path>/EMPM/bin)
        for d in $(sed -n 's/^! MKDIR: //p' "$f"); do
            [[ $d == io/* && $d != *..* ]] || { echo "bad MKDIR $d in $f"; exit 2; }
            mkdir -p "$d"
        done
        set +e; "$BIN" "$(basename "$f")" > log.txt 2>&1; st=$?; set -e
        if grep -qF -- "$exp" log.txt; then
            echo "error case $tag: ok (exit $st) ... $exp"
        else
            echo "error case $tag: FAILED (exit $st), expected: $exp"; tail -n 5 log.txt; exit 1
        fi
    done
    ;;
compare)
    python3 -I "$HERE/empm_compare.py" "$HERE/work/empm_baseline/io" "$HERE/work/empm_dev${TAG:-}_t$NT/io"
    ;;
*)
    echo "unknown mode $MODE"; exit 2 ;;
esac
