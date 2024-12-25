#!/bin/bash
#
# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Downloads and unzips the Small BFD database for AlphaFold.
#
# Usage: bash download_small_bfd.sh /path/to/download/directory
set -e

if [[ $# -eq 0 ]]; then
    echo "Error: download directory must be provided as an input argument."
    exit 1
fi

if ! command -v aria2c &> /dev/null ; then
    echo "Error: aria2c could not be found. Please install aria2c (sudo apt install aria2)."
    exit 1
fi

DOWNLOAD_DIR="$1"
ROOT_DIR="${DOWNLOAD_DIR}/small_bfd"
SOURCE_URL="https://storage.googleapis.com/alphafold-databases/reduced_dbs/bfd-first_non_consensus_sequences.fasta.gz"
BASENAME=$(basename "${SOURCE_URL}")

mkdir --parents "${ROOT_DIR}"

# Check if the file already exists
if [[ -f "${ROOT_DIR}/${BASENAME}" ]]; then
    echo "File ${ROOT_DIR}/${BASENAME} already exists. Skipping download."
else
    aria2c "${SOURCE_URL}" --dir="${ROOT_DIR}" --max-connection-per-server=16 --split=16 --min-split-size=1M
fi

# Check if the unzipped file exists
if [[ -f "${ROOT_DIR}/${BASENAME%.gz}" ]]; then
    echo "File ${ROOT_DIR}/${BASENAME%.gz} already exists. Skipping decompression."
else
    pushd "${ROOT_DIR}"
    gunzip "${ROOT_DIR}/${BASENAME}"
    popd
fi
