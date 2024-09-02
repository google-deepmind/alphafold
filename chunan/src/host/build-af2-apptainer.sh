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
BASE=$(realpath $(dirname $0))             # ~/alphafold/chunan/src/host
WD=$(dirname $(dirname $(dirname $BASE)))  # ~/alphafold
defFile=$(dirname $(dirname $BASE))/container/apptainer.runtime.def

pushd $WD
apptainer build $IMAGE_NAME $defFile
