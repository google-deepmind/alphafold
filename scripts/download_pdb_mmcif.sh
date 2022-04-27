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
# Downloads, unzips and flattens the PDB database for AlphaFold.
#
# Usage: bash download_pdb_mmcif.sh /path/to/download/directory
set -e

if [[ $# -eq 0 ]]; then
    echo "Error: download directory must be provided as an input argument."
    exit 1
fi

if ! command -v aria2c &> /dev/null ; then
    echo "Error: aria2c could not be found. Please install aria2c (sudo apt install aria2)."
    exit 1
fi

if ! command -v rsync &> /dev/null ; then
    echo "Error: rsync could not be found. Please install rsync."
    exit 1
fi

DOWNLOAD_DIR="$1"
ROOT_DIR="${DOWNLOAD_DIR}/pdb_mmcif"
RAW_DIR="${ROOT_DIR}/raw"
MMCIF_DIR="${ROOT_DIR}/mmcif_files"

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

# as we want to catch timeout errors, we turn abort-upon error off:
set +e

echo "Running rsync to fetch all mmCIF files (note that the rsync progress estimate might be inaccurate)..."
echo "If the download speed is too slow, try changing the mirror to:"
echo "  * rsync.ebi.ac.uk::pub/databases/pdb/data/structures/divided/mmCIF/ (Europe)"
echo "  * ftp.pdbj.org::ftp_data/structures/divided/mmCIF/ (Asia)"
echo "or see https://www.wwpdb.org/ftp/pdb-ftp-sites for more download options."
mkdir --parents "${RAW_DIR}"
rsync --recursive --links --perms --times --compress --info=progress2 --delete --port=33444 \
  rsync.rcsb.org::ftp_data/structures/divided/mmCIF/ \
  "${RAW_DIR}"

# we save this rsync's return code for multiple use
rsync_ret_code=$?

# now we test whether we have run into a timeout
# if this happens, give a reasonable hint to the user
if [ $rsync_ret_code -eq 10 ]; then
    echo
    echo "ERROR: rsync ran into a timeout. Possible reasons include"
    echo "       not setting the RSYNC_PROXY variable in an environment"
    echo "       where a webproxy is used."
    exit $rsync_ret_code
elif [ $rsync_ret_code -ne 0 ]; then
    echo 
    echo "ERROR: caught unkwon rsync-download error"
    exit $rsync_ret_code
fi

# finally, we turn abort-upon-error on again
set -e

# if we have pigz in PATH, we can attempt to decompress in parallel
if ! command -v unpigz &> /dev/null
then
    uncompress_cmd=gunzip
else
    uncompress_cmd=unpigz
fi

echo "Unzipping all mmCIF files..."
# unsure how many processors may be used, yet 2 is faster than 1 in any case
find "${RAW_DIR}/" -type f -iname "*.gz" -print0 | xargs -0 -P2 "${uncompress_cmd}"

echo "Flattening all mmCIF files..."
mkdir --parents "${MMCIF_DIR}"
find "${RAW_DIR}" -type d -empty -delete  # Delete empty directories.
for subdir in "${RAW_DIR}"/*; do
  mv "${subdir}/"*.cif "${MMCIF_DIR}"
done

# Delete empty download directory structure.
find "${RAW_DIR}" -type d -empty -delete

aria2c "ftp://ftp.wwpdb.org/pub/pdb/data/status/obsolete.dat" --dir="${ROOT_DIR}"
