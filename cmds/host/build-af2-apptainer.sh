#!/bin/zsh

# Aim: build AF2 image
# Input:
#  - Dockerfile.runtime
#  - .devcontainer/
#  - alphafold/
# Output:
# Usage:
# Example:
# Dependencies:

set -e

# allow change name from cli
IMAGE_NAME=${1:-"alphafold2.3.sif"}

# Set variables
BASE=$(realpath $(dirname $0))  # ~/alphafold/cmds/host
WD=$(dirname $(dirname $BASE))  # /host/path/to/alphafold

pushd $WD
apptainer build $IMAGE_NAME .devcontainer/apptainer.runtime.def
