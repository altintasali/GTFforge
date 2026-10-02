#!/usr/bin/env bash
# Guard 06: unstranded groups use gffcompare -M, stranded groups keep stranded single exons
#
# Run on its own:   .tests/guards/06_single_exon_policy_follows_strandedness.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

fixture_rundir
sed -i 's|,unstranded$|,reverse|' "$R/input/samples.csv"
sed -i 's|^\(trt_[0-9]\),\(.*\),reverse$|\1,\2,unstranded|' "$R/input/samples.csv"
(cd "$R" && snakemake -n -p --cores 1 > run.log 2>&1)
ctrl="$(grep -E 'gffcompare .*support/ctrl/cmp ' "$R/run.log")"
trt="$(grep -E 'gffcompare .*support/trt/cmp ' "$R/run.log")"
case "$ctrl" in *" -M "*) echo "ERROR: stranded group ctrl still drops single exons"; FAIL=1 ;; esac
case "$trt" in *" -M "*) ;; *) echo "ERROR: unstranded group trt does not use -M"; FAIL=1 ;; esac
[ -n "$ctrl" ] && [ -n "$trt" ] || { echo "ERROR: gffcompare commands not found"; FAIL=1; }

exit $FAIL
