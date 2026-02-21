#!/bin/bash
set -e

MGF_FILE="$1"
YAML_FILE="$2"
CKPT_FILE="$3"

echo "Running Casanovo"
echo "MGF: $MGF_FILE"
echo "YAML: $YAML_FILE"
echo "Checkpoint: $CKPT_FILE"
# Run Casanovo inside the container
singularity exec docker://mriffle/casanovo:4.2.1 \
  casanovo sequence \
    -o ${MGF_FILE%.mgf}.mztab \
    -m "$CKPT_FILE" \
    -c "$YAML_FILE" \
    "$MGF_FILE"
