#!/bin/bash

INPUT="$1"
BASENAME=$(basename "$INPUT" .mgf)
OUTPUT="${BASENAME}_decoy.mgf"

python3 generate_decoy.py $INPUT $OUTPUT
