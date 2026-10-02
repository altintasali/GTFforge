#!/usr/bin/env bash
# Guard 07: filter.keep_contigs_regex removes scaffold transcripts, reference and novel
#
# Run on its own:   .tests/guards/07_contig_regex_drops_unplaced_scaffolds.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

fixture_rundir "$(printf 'filter:\n  keep_classes: "jkxuioy"\n  keep_contigs_regex: "^chr"')"
if ! run_workflow; then
  echo "ERROR: run failed"; tail -20 "$R/run.log"; FAIL=1
elif grep -q '^unplacedT' "$R/results/gtfforge.gtf"; then
  echo "ERROR: unplacedT features survived the contig filter"; FAIL=1
elif [ "$(novel_fates | grep -c .)" -ne 4 ]; then
  echo "ERROR: expected 4 novel transcripts without the scaffold one"; novel_fates; FAIL=1
fi

exit $FAIL
