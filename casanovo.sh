#!/bin/bash
set -e

MGF_FILE="$(basename "$1")"
YAML_FILE="$(basename "$2")"
CKPT_FILE="$(basename "$3")"
OUT_DIR="$4"
SAMPLE="$5"

# Define the expected output name
MZTAB_OUT="${SAMPLE}.mztab"

echo "Running Casanovo"
echo "MGF: $MGF_FILE"
echo "YAML: $YAML_FILE"
echo "Checkpoint: $CKPT_FILE"
# Run Casanovo inside the container
singularity exec docker://mriffle/casanovo:4.2.1 \
  casanovo sequence \
    -o "$MZTAB_OUT" \
    -m "$CKPT_FILE" \
    -c "$YAML_FILE" \
    "$MGF_FILE"
