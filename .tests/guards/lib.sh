#!/usr/bin/env bash
# Shared setup and fixtures for the guard scripts in this directory.
#
# Each guard is standalone: its OWN temp directory ($T), its own failures in
# $FAIL, non-zero exit if any check failed. Guards never depend on each
# other's leftovers; shared fixtures are functions here.

GUARD_REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$GUARD_REPO_ROOT" || exit 1

guard_init() {
  T="$(mktemp -d)"
  export T
  # shellcheck disable=SC2034  # read and set by the sourcing guard
  FAIL=0
  if [ -z "${GUARD_KEEP_TMP:-}" ]; then
    trap 'rm -rf "$T"' EXIT
  else
    trap 'echo "[guard] workdir kept: $T" >&2' EXIT
  fi
}

# A standalone analysis directory at $T/run running the bundled test data:
# workflow/ linked to this checkout, input/config.yaml + input/samples.csv
# with absolute paths. Extra lines in $1 (optional) are appended to the
# config, so a guard can override any key.
fixture_rundir() {
  local extra="${1:-}"
  R="$T/run"
  export R
  mkdir -p "$R/input"
  ln -s "$GUARD_REPO_ROOT/workflow" "$R/workflow"
  sed -e "s|\.tests/|$GUARD_REPO_ROOT/.tests/|" .tests/samples.csv > "$R/input/samples.csv"
  sed -e "s|^samples: .*|samples: input/samples.csv|" \
      -e "s|\.tests/|$GUARD_REPO_ROOT/.tests/|" \
      config/test.yaml > "$R/input/config.yaml"
  if [ -n "$extra" ]; then
    printf '%s\n' "$extra" >> "$R/input/config.yaml"
  fi
}

# Run the workflow in $R; returns snakemake's exit status, log in $R/run.log.
run_workflow() {
  (cd "$R" && snakemake --cores 2 "$@" > run.log 2>&1)
}

# The novel transcripts of a finished run as "class<TAB>support_groups",
# one per line, sorted -- MSTRG ids are not stable, structures are.
novel_fates() {
  awk -F'\t' 'NR > 1 {print $5 "\t" $7}' "$R/results/novel_transcripts.tsv" | sort
}
