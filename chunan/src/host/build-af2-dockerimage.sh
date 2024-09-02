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
BASE=$(realpath $(dirname $0))             # ~/alphafold/chunan/src/host
WD=$(dirname $(dirname $(dirname $BASE)))  # ~/alphafold
dockerFile=$(dirname $(dirname $BASE))/container/Dockerfile.runtime

pushd $WD
docker build --tag $DOCKER_IMAGE_NAME -f $dockerFile .
