#!/bin/bash
# Convert Dockerfile to a Singularity recipe.
# Requires Singularity Python (spython)

TMPRECIPE=$( mktemp )

spython recipe Dockerfile > $TMPRECIPE 2> /dev/null

cat $TMPRECIPE \
    | grep -v ^CUDA \
    | grep -v ^SHELL \
    | sed -e 's;\${CUDA/\./-};11-0;' \
    | sed -e 's;\${CUDA/\./};110;'\
    | sed -e 's;\${CUDA_VERSION};11.0.3;' \
    > Singularity.def

/bin/rm -f $TMPRECIPE
