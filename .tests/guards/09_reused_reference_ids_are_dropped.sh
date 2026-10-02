#!/usr/bin/env bash
# Guard 09: a transcript_id reused on two loci (UCSC refGene style) is dropped and the GTF still validates
#
# Run on its own:   .tests/guards/09_reused_reference_ids_are_dropped.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

fixture_rundir
{ cat .tests/resources/reference.gtf
  printf 'chrT\tHAVANA\texon\t33001\t33200\t.\t+\t.\tgene_id "ENSGT0000004.1"; transcript_id "ENSTT0000004.1"; exon_number "1";\n'
} > "$T/dup.gtf"
sed -i "s|  gtf: .*|  gtf: $T/dup.gtf|" "$R/input/config.yaml"
if ! run_workflow; then
  echo "ERROR: run failed (validation should have passed after dropping)"; tail -20 "$R/run.log"; FAIL=1
else
  if grep -q 'ENSTT0000004.1' "$R/results/gtfforge.gtf"; then
    echo "ERROR: the reused transcript_id survived"; FAIL=1
  fi
  if ! grep -q '"ref_transcript_multi_locus": 1' "$R/results/finalize/finalize_stats.json"; then
    echo "ERROR: the drop is not reported"; FAIL=1
  fi
fi

exit $FAIL
