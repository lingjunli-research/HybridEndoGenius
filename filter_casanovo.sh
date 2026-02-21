#!/bin/bash
set -euo pipefail

echo "Running filter_casanovo.py"

TARGET_MZTAB="$1"
DECOY_MZTAB="$2"
FDR_VALUE="$3"
DATABASE="$4"

python filter_casanovo.py \
  --target "${TARGET_MZTAB}" \
  --decoy "${DECOY_MZTAB}" \
  --score_col "search_engine_score[1]" \
  --fdr "${FDR_VALUE}" \
  --outdir results \
  --db "${DATABASE}"

echo "Finished filter_casanovo.py"
