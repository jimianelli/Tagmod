#!/bin/sh
set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
project_dir=$(CDPATH= cd -- "$script_dir/.." && pwd)

"$script_dir/build.sh"
run_output=$("$script_dir/run-model.sh" both)
result_dir=$(printf '%s\n' "$run_output" | sed -n 's/^Results: //p')

if [ -z "$result_dir" ] || [ ! -f "$result_dir/tm.rep" ]; then
  echo "error: model did not produce tm.rep" >&2
  exit 1
fi

biomass=$(awk '$0 == "Biomass" {getline; print; exit}' "$result_dir/tm.rep")
if ! awk -v value="$biomass" 'BEGIN {exit !(value + 0 > 0)}'; then
  echo "error: expected positive biomass, got: $biomass" >&2
  exit 1
fi

echo "Smoke test passed; Biomass=$biomass"
echo "Results: $result_dir"
