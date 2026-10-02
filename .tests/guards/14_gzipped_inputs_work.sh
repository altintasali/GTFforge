#!/usr/bin/env bash
# Guard 14: gzipped reference, TE GTF and sample assemblies are read transparently
#
# Run on its own:   .tests/guards/14_gzipped_inputs_work.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

fixture_rundir
for f in reference.gtf te.gtf genome.fa; do gzip -c ".tests/resources/$f" > "$T/$f.gz"; done
mkdir -p "$T/asm"; for f in .tests/assemblies/*.gtf; do gzip -c "$f" > "$T/asm/$(basename "$f").gz"; done
sed -i -e "s|  gtf: .*|  gtf: $T/reference.gtf.gz|" -e "s|  te_gtf: .*|  te_gtf: $T/te.gtf.gz|" \
       -e "s|  fasta: .*|  fasta: $T/genome.fa.gz|" "$R/input/config.yaml"
sed -i -E "s|,[^,]*/assemblies/([^,]*\.gtf),|,$T/asm/\1.gz,|" "$R/input/samples.csv"
if ! run_workflow; then
  echo "ERROR: run with gzipped inputs failed"; tail -20 "$R/run.log"; FAIL=1
elif [ "$(novel_fates | grep -c .)" -ne 5 ]; then
  echo "ERROR: gzipped inputs changed the result"; novel_fates; FAIL=1
fi

exit $FAIL
