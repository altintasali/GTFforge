#!/usr/bin/env bash
# Guard 10: validate_gtf.py fails on a broken GTF and names the problem
#
# Run on its own:   .tests/guards/10_validation_rejects_broken_gtfs.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

printf 'chr1\tx\texon\t10\t20\t.\t+\t.\tgene_id "G"; transcript_id "T";\nchr2\tx\texon\t10\t20\t.\t+\t.\tgene_id "G"; transcript_id "T";\n' > "$T/bad.gtf"
if python3 workflow/scripts/validate_gtf.py --gtf "$T/bad.gtf" --out "$T/v.json" 2> "$T/err"; then
  echo "ERROR: a transcript on two chromosomes passed validation"; FAIL=1
elif ! grep -q transcript_on_two_loci "$T/err"; then
  echo "ERROR: the failure does not name the check"; cat "$T/err"; FAIL=1
fi

exit $FAIL
