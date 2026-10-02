#!/usr/bin/env bash
# Guard 03: a group with a single sample is rejected
#
# Run on its own:   .tests/guards/03_group_with_one_sample_is_rejected.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

fixture_rundir
sed -i 's|^ctrl_1,ctrl,|ctrl_1,lonely,|' "$R/input/samples.csv"
if (cd "$R" && snakemake -n --cores 1 > run.log 2>&1); then
  echo "ERROR: dry-run accepted a one-sample group"; FAIL=1
elif ! grep -q "lonely (1)" "$R/run.log"; then
  echo "ERROR: the message does not name the group"; tail -5 "$R/run.log"; FAIL=1
fi

exit $FAIL
