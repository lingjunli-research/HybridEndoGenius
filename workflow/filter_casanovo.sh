#!/bin/bash
set -euo pipefail

echo "PWD:"
pwd

echo "Running filter_casanovo.py"

TARGET="$(basename "$1")"
DECOY="$(basename "$2")"
FDR="$3"
DATABASE="$(basename "$4")"

mkdir -p "$PWD/filter"

python /app/filter_casanovo.py \
  --target "${TARGET}" \
  --decoy "${DECOY}" \
  --score_col "search_engine_score[1]" \
  --fdr "${FDR}" \
  --outdir "$PWD/filter" \
  --db "${DATABASE}"

echo "Finished filter_casanovo.py"
