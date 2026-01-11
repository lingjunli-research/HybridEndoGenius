#!/bin/bash

# Run Casanovo inside the container
singularity exec docker://mriffle/casanovo:4.2.1 \
  casanovo sequence \
    -o ${1%.mgf}.mztab \
    -m casanovo_v4_2_0.ckpt \
    -c casanovo.yaml \
    $1
