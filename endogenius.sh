#!/bin/bash
set -e

# Run EndoGenius with JSON config
singularity exec endogenius.sif python /app/EndoGenius/endogenius_cli.py --config $1
