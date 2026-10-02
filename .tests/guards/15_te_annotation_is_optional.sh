#!/usr/bin/env bash
# Guard 15: without ref.te_gtf the TE stage is skipped and the table has no TE columns
#
# Run on its own:   .tests/guards/15_te_annotation_is_optional.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

fixture_rundir
sed -i 's|  te_gtf: .*|  te_gtf: ""|' "$R/input/config.yaml"
if ! run_workflow; then
  echo "ERROR: run without te_gtf failed"; tail -20 "$R/run.log"; FAIL=1
elif head -1 "$R/results/novel_transcripts.tsv" | grep -q te_at_tss; then
  echo "ERROR: TE columns present without a TE annotation"; FAIL=1
fi

exit $FAIL
