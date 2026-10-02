#!/usr/bin/env bash
# Guard 11: gtfforge samples rebuilds the GTF rows of .tests/samples.csv
#
# Run on its own:   .tests/guards/11_cli_samples_reproduces_the_fixture.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

python3 workflow/scripts/gtfforge samples --gtf-dir .tests/assemblies --dry-run \
  > "$T/sheet.csv" 2> "$T/err" || { echo "ERROR: samples failed"; cat "$T/err"; FAIL=1; }
norm() { grep -v '^#' "$1" | grep -v '^sample,' | grep -v ',bam,' \
  | awk -F, '{n=split($3,p,"/"); print $1","$2","p[n]","$5}' | sort; }
if ! diff <(norm "$T/sheet.csv") <(norm .tests/samples.csv) > "$T/diff"; then
  echo "ERROR: generated sheet differs from the fixture"; cat "$T/diff"; FAIL=1
fi

exit $FAIL
