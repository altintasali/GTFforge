#!/usr/bin/env bash
# Guard 12: gtfforge init + run --dry-run work in a directory that is not a checkout
#
# Run on its own:   .tests/guards/12_cli_run_in_a_standalone_directory.sh
# Run all guards:   .tests/guards/run.sh
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
guard_init

mkdir -p "$T/analysis"
python3 workflow/scripts/gtfforge init --input-dir "$T/analysis/input" 2> "$T/init.err"
cp .tests/samples.csv "$T/analysis/input/samples.csv"
sed -i "s|\.tests/|$GUARD_REPO_ROOT/.tests/|" "$T/analysis/input/samples.csv"
sed -e "s|^samples: .*|samples: input/samples.csv|" -e "s|\.tests/|$GUARD_REPO_ROOT/.tests/|" \
  config/test.yaml > "$T/analysis/input/config.yaml"
if ! python3 workflow/scripts/gtfforge run --directory "$T/analysis" --dry-run > "$T/run.log" 2>&1; then
  echo "ERROR: run --dry-run failed"; tail -20 "$T/run.log"; FAIL=1
elif [ ! -L "$T/analysis/workflow" ]; then
  echo "ERROR: workflow/ was not linked into the analysis directory"; FAIL=1
fi

exit $FAIL
