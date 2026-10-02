#!/usr/bin/env bash
# Guard 02: a sample row with both gtf and bam is rejected
#
# Run on its own:   .tests/guards/02_gtf_and_bam_in_one_row_is_rejected.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

fixture_rundir
sed -i 's|^ctrl_1,ctrl,\([^,]*\),,|ctrl_1,ctrl,\1,/x.bam,|' "$R/input/samples.csv"
if (cd "$R" && snakemake -n --cores 1 > run.log 2>&1); then
  echo "ERROR: dry-run accepted a row with both gtf and bam"; FAIL=1
elif ! grep -q "exactly one of 'gtf'" "$R/run.log"; then
  echo "ERROR: no clear message for gtf+bam"; tail -5 "$R/run.log"; FAIL=1
fi

exit $FAIL
