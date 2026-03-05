#!/bin/bash
set -e

JSON="$(basename "$1")"
# Run EndoGenius with JSON config
echo "Running EndoGenius inside Docker container"
echo "JSON: $JSON"
python /app/EndoGenius/endogenius_cli.py --config "$JSON"
