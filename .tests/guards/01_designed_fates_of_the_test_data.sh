#!/usr/bin/env bash
# Guard 01: every novel structure of the test data gets its designed fate
#
# Run on its own:   .tests/guards/01_designed_fates_of_the_test_data.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

# The synthetic data (.tests/generate_test_data.py) encodes one case per
# filter. Any change to support, merge, class or contig handling shows up
# here as a different set of kept transcripts.
fixture_rundir
if ! run_workflow; then
  echo "ERROR: test run failed"; tail -30 "$R/run.log"; FAIL=1
else
  expected="$(printf 'j\tbam,ctrl,trt\nu\tbam\nu\tctrl\nu\ttrt\nx\tctrl\n')"
  got="$(novel_fates)"
  if [ "$got" != "$expected" ]; then
    echo "ERROR: novel transcripts differ from the design"
    diff <(echo "$expected") <(echo "$got"); FAIL=1
  fi
  # N_j is a new isoform of G1 and must keep G1's gene_id
  if ! awk -F'\t' '$5=="j" && $2=="ENSGT0000001.1"' "$R/results/novel_transcripts.tsv" | grep -q .; then
    echo "ERROR: the j isoform did not join its reference gene"; FAIL=1
  fi
  # N_rare (1 replicate) and the unstranded single exon never reach the GTF
  if grep -q '24001' "$R/results/gtfforge.gtf" || grep -q '27001' "$R/results/gtfforge.gtf"; then
    echo "ERROR: an unsupported or single-exon transfrag reached the GTF"; FAIL=1
  fi
fi

exit $FAIL
