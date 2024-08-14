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
DOCKER_IMAGE_NAME=${1:-"$USER/alphafold2.3:base"}

# Set variables
WD=$(dirname $(dirname $(dirname $(realpath $0))))  # /host/path/to/alphafold/

pushd $WD
docker build --tag $DOCKER_IMAGE_NAME -f .devcontainer/Dockerfile.runtime .
