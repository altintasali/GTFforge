#!/usr/bin/env bash
# Guard 13: the MultiQC report carries every GTFforge section
#
# Run on its own:   .tests/guards/13_report_has_every_section.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

fixture_rundir
if ! run_workflow; then
  echo "ERROR: run failed"; tail -20 "$R/run.log"; FAIL=1
else
  for s in overview support_histogram support_kept novel_classes novel_by_group dropped \
           accuracy te_families config_used; do
    if ! compgen -G "$R/results/qc/multiqc_report_data/multiqc_gtfforge_${s}_*.txt" > /dev/null; then
      echo "ERROR: report section missing: $s"; FAIL=1
    fi
  done
  if ! grep -q "GTFforge" "$R/results/qc/multiqc_report_data/multiqc_software_versions.txt"; then
    echo "ERROR: software versions do not list GTFforge"; FAIL=1
  fi
fi

exit $FAIL
