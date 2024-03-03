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

# Set variables
WD=$(dirname $(dirname $(realpath $0)))  # /host/path/to/alphafold/
DOCKER_IMAGE_NAME="$USER/alphafold2.3:base"

pushd $WD
docker build --tag $DOCKER_IMAGE_NAME -f .devcontainer/Dockerfile.runtime .
