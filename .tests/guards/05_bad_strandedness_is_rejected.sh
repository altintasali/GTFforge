#!/usr/bin/env bash
# Guard 05: an unknown strandedness value fails schema validation
#
# Run on its own:   .tests/guards/05_bad_strandedness_is_rejected.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

fixture_rundir
sed -i 's|^ctrl_1,\(.*\),unstranded$|ctrl_1,\1,sideways|' "$R/input/samples.csv"
if (cd "$R" && snakemake -n --cores 1 > run.log 2>&1); then
  echo "ERROR: dry-run accepted strandedness 'sideways'"; FAIL=1
fi

exit $FAIL
