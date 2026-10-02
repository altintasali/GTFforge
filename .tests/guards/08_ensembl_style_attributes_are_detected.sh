#!/usr/bin/env bash
# Guard 08: an Ensembl-style reference keeps gene_biotype and needs no extra nf-core flag
#
# Run on its own:   .tests/guards/08_ensembl_style_attributes_are_detected.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

fixture_rundir
sed -e 's/gene_type/gene_biotype/; s/transcript_type/transcript_biotype/' \
  .tests/resources/reference.gtf > "$T/ensembl.gtf"
sed -i "s|  gtf: .*|  gtf: $T/ensembl.gtf|" "$R/input/config.yaml"
if ! run_workflow; then
  echo "ERROR: run failed"; tail -20 "$R/run.log"; FAIL=1
else
  if ! grep 'novel "1"' "$R/results/gtfforge.gtf" | grep -q 'gene_biotype "novel"'; then
    echo "ERROR: novel genes do not use gene_biotype"; FAIL=1
  fi
  if grep -q gene_type "$R/results/gtfforge.gtf"; then
    echo "ERROR: gene_type leaked into an Ensembl-style GTF"; FAIL=1
  fi
  if grep -q featurecounts_group_type "$R/results/README_nfcore.txt"; then
    echo "ERROR: README asks for --featurecounts_group_type with gene_biotype"; FAIL=1
  fi
fi

exit $FAIL
