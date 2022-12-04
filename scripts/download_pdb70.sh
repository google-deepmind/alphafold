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
# Downloads and unzips the PDB70 database for AlphaFold.
#
# Usage: bash download_pdb70.sh /path/to/download/directory
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
ROOT_DIR="${DOWNLOAD_DIR}/pdb70"
SOURCE_URL="http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/old-releases/pdb70_from_mmcif_200401.tar.gz"
BASENAME=$(basename "${SOURCE_URL}")

if [ -d "${ROOT_DIR}" ]; then
    echo "WARNING: Destination directory '${ROOT_DIR}' does already exist."
    read -p "Proceed by deleting existing download directory? [Y/n]" -n1 -r
    echo
    if [[ ! $REPLY =~ ^[Nn]$ ]]; then
        echo "INFO: Deleting previous download directory: '${ROOT_DIR}'"
        rm -rf "${ROOT_DIR}"
    else
        echo "Aborting download."
        exit 0
    fi
fi

mkdir --parents "${ROOT_DIR}"
aria2c "${SOURCE_URL}" --dir="${ROOT_DIR}"

if ! command -v pigz &> /dev/null
then
    tar --extract --verbose --file="${ROOT_DIR}/${BASENAME}" \
        --directory="${ROOT_DIR}"
else
    tar -I pigz --extract --verbose --file="${ROOT_DIR}/${BASENAME}" \
        --directory="${ROOT_DIR}"
fi

# the files are only user-readable. For multi-user systems we need to correct:
chmod 444 "${ROOT_DIR}"/*

rm "${ROOT_DIR}/${BASENAME}"
