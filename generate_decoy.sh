#!/bin/bash
set -e

INPUT="$(basename "$1")"
SAMPLE="$2"
OUTPUT="${SAMPLE}_decoy.mgf"

python3 generate_decoy.py $INPUT $OUTPUT
