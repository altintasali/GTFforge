#!/usr/bin/env bash
# Guard 04: a min_samples larger than a group warns at parse time
#
# Run on its own:   .tests/guards/04_unreachable_support_threshold_warns.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

fixture_rundir "$(printf 'support:\n  min_samples: 5\n  min_fraction: 0.33\n  single_exon: auto')"
if ! (cd "$R" && snakemake -n --cores 1 > run.log 2>&1); then
  echo "ERROR: dry-run failed"; tail -10 "$R/run.log"; FAIL=1
elif ! grep -q "none of their transcripts can pass" "$R/run.log"; then
  echo "ERROR: no warning for an unreachable threshold"; FAIL=1
fi

exit $FAIL
