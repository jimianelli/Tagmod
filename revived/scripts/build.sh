#!/bin/sh
set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
project_dir=$(CDPATH= cd -- "$script_dir/.." && pwd)
src_dir="$project_dir/src"

if ! command -v admb >/dev/null 2>&1; then
  echo "error: ADMB is not available on PATH" >&2
  exit 127
fi

echo "Building tm from $src_dir/tm.tpl"
(cd "$src_dir" && admb -f tm.tpl)
file "$src_dir/tm"
