#!/bin/sh
set -eu

usage() {
  echo "usage: $0 DATASET [--mcmc]" >&2
  echo "       DATASET is a .dat filename stem in the data directory" >&2
}

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
  usage
  exit 2
fi

dataset=${1%.dat}
mode=mpd
if [ "$#" -eq 2 ]; then
  if [ "$2" != "--mcmc" ]; then
    usage
    exit 2
  fi
  mode=mcmc
fi

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
project_dir=$(CDPATH= cd -- "$script_dir/.." && pwd)
data_file="$project_dir/data/$dataset.dat"
executable="$project_dir/src/tm"

if [ ! -f "$data_file" ]; then
  echo "error: dataset not found: $data_file" >&2
  exit 2
fi

if [ ! -x "$executable" ]; then
  "$script_dir/build.sh"
fi

timestamp=$(date -u '+%Y%m%dT%H%M%SZ')
result_dir="$project_dir/runs/generated/${dataset}_${mode}_$timestamp"
work_dir="$result_dir/work"
mkdir -p "$work_dir"

cp "$data_file" "$work_dir/input.dat"
printf '%s\n' "input.dat" > "$work_dir/tm.dat"

echo "Running $dataset ($mode)"
if [ "$mode" = "mcmc" ]; then
  iterations=${MCMC_ITERATIONS:-1000000}
  save_interval=${MCMC_SAVE:-200}
  (cd "$work_dir" && "$executable" -nox -iprint 50 \
    -mcmc "$iterations" -mcsave "$save_interval")
  (cd "$work_dir" && "$executable" -nox -iprint 50 -mceval)
else
  (cd "$work_dir" && "$executable" -nox -iprint 50)
fi

for output in tm.par tm.std tm.cor tm.rep tm.psv mcout.rep checkfile.rep; do
  if [ -f "$work_dir/$output" ]; then
    cp "$work_dir/$output" "$result_dir/$output"
  fi
done

cat > "$result_dir/run-metadata.txt" <<EOF
dataset=$dataset
mode=$mode
created_utc=$timestamp
admb=$(command -v admb)
mcmc_iterations=${MCMC_ITERATIONS:-1000000}
mcmc_save=${MCMC_SAVE:-200}
EOF

echo "Results: $result_dir"
