#!/bin/bash
set -euo pipefail

echo "Running filter_casanovo.py"

TARGET_MZTAB="$1"
DECOY_MZTAB="$2"

python filter_casanovo.py \
  --target "${TARGET_MZTAB}" \
  --decoy "${DECOY_MZTAB}" \
  --score_col "search_engine_score[1]" \
  --fdr 0.05 \
  --outdir results \
  --db1 motif_DB.fasta \
  --db2 new_NPs_since_2010.fasta \
  --db3 np_database.fasta

echo "End"
