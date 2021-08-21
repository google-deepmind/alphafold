#!/bin/bash
# Description: AlphaFold non-docker version
# Author: Sanjay Kumar Srikakulam
# https://github.com/kalininalab/alphafold_non_docker/blob/main/run_alphafold.sh
ldconfig

# Author: WTTAT
usage() {
    echo " "
    echo "something wrong"
    echo " "
    exit 1
}

while getopts ":f:m:d:p:b" i; do
        case "${i}" in
        f)
                fasta_paths=$OPTARG
        ;;
        m)
                model_names=$OPTARG
        ;;
        d)
                max_template_date=$OPTARG
        ;;
        p)
                preset=$OPTARG
        ;;
        b)
                benchmark=$OPTARG
        ;;
        esac
done

echo "BATCH_BUCKET : $BATCH_BUCKET"
echo "REGION : $REGION"

# Parse input and set defaults
if [[ "$model_names" == "" || "$fasta_paths" == "" || "$max_template_date" == "" ]] ; then
    usage
fi

if [[ "$benchmark" == "" ]] ; then
    benchmark=false
fi

if [[ "$preset" == "" ]] ; then
    preset="full_dbs"
fi

if [[ "$preset" != "full_dbs" && "$preset" != "casp14" && "$preset" != "reduced_dbs" ]] ; then
    echo "Unknown preset! Using default ('full_dbs')"
    preset="full_dbs"
fi


# This bash script looks for the run_alphafold.py script in its current working directory, if it does not exist then exits
# current_working_dir=$(pwd)
# alphafold_script="$current_working_dir/run_alphafold.py"

# Export ENVIRONMENT variables and set CUDA devices for use
# CUDA GPU control

# no need for batch
# export CUDA_VISIBLE_DEVICES=-1
# if [[ "$use_gpu" == true ]] ; then
#     export CUDA_VISIBLE_DEVICES=0

#     if [[ "$gpu_devices" ]] ; then
#         export CUDA_VISIBLE_DEVICES=$gpu_devices
#     fi
# fi

# # OpenMM threads control
# if [[ "$openmm_threads" ]] ; then
#     export OPENMM_CPU_THREADS=$openmm_threads
# fi

# data_dir="/mnt/dataset/"
data_dir="/fsx/dataset/"
# Path and user config (change me if required)

bfd_database_path="$data_dir/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
small_bfd_database_path="$data_dir/small_bfd/bfd-first_non_consensus_sequences.fasta"
# mgnify_database_path="$data_dir/mgnify/mgy_clusters.fa"
# template_mmcif_dir="$data_dir/pdb_mmcif/mmcif_files"
# obsolete_pdbs_path="$data_dir/pdb_mmcif/obsolete.dat"
# pdb70_database_path="$data_dir/pdb70/pdb70"
uniclust30_database_path="$data_dir/uniclust30/uniclust30_2018_08/uniclust30_2018_08"
# uniref90_database_path="$data_dir/uniref90/uniref90.fasta"

# Binary path (change me if required)
# hhblits_binary_path=$(which hhblits)
# hhsearch_binary_path=$(which hhsearch)
# jackhmmer_binary_path=$(which jackhmmer)
# kalign_binary_path=$(which kalign)



######
#  if S3 URL，download fasta，change to file name
#  only support one file
#  by WTTAT
echo "start downloading"
aws s3 cp s3://input/$fasta_paths ./ --region $REGION
# fasta_paths="${fasta_paths##*/}"
# echo "fasta_paths changed to $fasta_paths"

# ######
echo "start running af2"
# Run AlphaFold with required parameters
# 'reduced_dbs' preset does not use bfd and uniclust30 databases
if [[ "$preset" == "reduced_dbs" ]]; then
    $(python /app/alphafold/run_alphafold.py --BATCH_BUCKET="$BATCH_BUCKET" --small_bfd_database_path="$small_bfd_database_path" --fasta_paths="$fasta_paths" --model_names="$model_names" --max_template_date="$max_template_date" --preset="$preset" --benchmark="$benchmark" --logtostderr)
else
    $(python /app/alphafold/run_alphafold.py  --BATCH_BUCKET="$BATCH_BUCKET"  --bfd_database_path="$bfd_database_path" --uniclust30_database_path="$uniclust30_database_path" --fasta_paths="$fasta_paths" --model_names="$model_names" --max_template_date="$max_template_date" --preset="$preset" --benchmark="$benchmark" --logtostderr)
fi