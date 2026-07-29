#!/bin/sh
set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
project_dir=$(CDPATH= cd -- "$script_dir/.." && pwd)

datasets="
Seamounts
Seamounts_no_RF
Seguam
Seguam1
Seguam1b
Seguam34
Seguam_M
Seguam_no_RF
Seguam_no_RF_M
both
both_no_RF
nearshore
nearshore_no_RF
"

for dataset in $datasets; do
  "$script_dir/run-model.sh" "$dataset"
done

Rscript "$project_dir/R/validate-archive.R" "$project_dir"
