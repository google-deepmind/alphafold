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
# Downloads and unzips all required data for AlphaFold.
#
# Usage: bash download_all_data.sh /path/to/download/directory
set -e

if [[ $# -eq 0 ]]; then
    echo "Error: download directory must be provided as an input argument."
    exit 1
fi

if ! command -v rsync &> /dev/null ; then
    echo "Error: rsync could not be found. Please install rsync."
    exit
fi

if ! command -v aria2c &> /dev/null ; then
    echo "Error: aria2c could not be found. Please install aria2c (sudo apt install aria2)."
    exit
fi

DOWNLOAD_DIR="$1"
SCRIPT_DIR="$(dirname "$(realpath "$0")")"

echo "Downloading AlphaFold parameters..."
bash "${SCRIPT_DIR}/download_alphafold_params.sh" "${DOWNLOAD_DIR}"

echo "Downloading BFD..."
bash "${SCRIPT_DIR}/download_bfd.sh" "${DOWNLOAD_DIR}"

echo "Downloading MGnify..."
bash "${SCRIPT_DIR}/download_mgnify.sh" "${DOWNLOAD_DIR}"

echo "Downloading PDB70..."
bash "${SCRIPT_DIR}/download_pdb70.sh" "${DOWNLOAD_DIR}"

echo "Downloading PDB mmCIF files..."
bash "${SCRIPT_DIR}/download_pdb_mmcif.sh" "${DOWNLOAD_DIR}"

echo "Downloading Uniclust30..."
bash "${SCRIPT_DIR}/download_uniclust30.sh" "${DOWNLOAD_DIR}"

echo "Downloading Uniref90..."
bash "${SCRIPT_DIR}/download_uniref90.sh" "${DOWNLOAD_DIR}"

echo "All data downloaded."
