#!/bin/bash
# Description: AlphaFold non-docker version
# Author: Sanjay Kumar Srikakulam
# https://github.com/kalininalab/alphafold_non_docker/blob/main/run_alphafold.sh
ldconfig

# Author: WTTAT
usage() {
    echo " "
    echo "六个环境变量"
    echo "model names : $model_names"
    echo "fasta_path : $fasta_path"
    echo "max_template_date : $max_template_date"
    echo "preset : $preset"
    echo "benchmark : $benchmark"
    echo "BATCH_BUCKET : $BATCH_BUCKET"
    echo " "
    exit 1
}
# Parse input and set defaults
if [[ "$model_names" == "" || "$fasta_path" == "" || "$max_template_date" == "" ]] ; then
    usage
fi

if [[ "$benchmark" == "" ]] ; then
    benchmark=false
fi

if [[ "$use_gpu" == "" ]] ; then
    use_gpu=true
fi

if [[ "$gpu_devices" == "" ]] ; then
    gpu_devices=0
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

# TensorFlow control
export TF_FORCE_UNIFIED_MEMORY='1'

# JAX control
export XLA_PYTHON_CLIENT_MEM_FRACTION='4.0'

data_dir="/mnt/dataset/"
# Path and user config (change me if required)

bfd_database_path="$data_dir/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
small_bfd_database_path="$data_dir/small_bfd/bfd-first_non_consensus_sequences.fasta"
mgnify_database_path="$data_dir/mgnify/mgy_clusters.fa"
template_mmcif_dir="$data_dir/pdb_mmcif/mmcif_files"
obsolete_pdbs_path="$data_dir/pdb_mmcif/obsolete.dat"
pdb70_database_path="$data_dir/pdb70/pdb70"
uniclust30_database_path="$data_dir/uniclust30/uniclust30_2018_08/uniclust30_2018_08"
uniref90_database_path="$data_dir/uniref90/uniref90.fasta"

# Binary path (change me if required)
# hhblits_binary_path=$(which hhblits)
# hhsearch_binary_path=$(which hhsearch)
# jackhmmer_binary_path=$(which jackhmmer)
# kalign_binary_path=$(which kalign)

# Run AlphaFold with required parameters
# 'reduced_dbs' preset does not use bfd and uniclust30 databases
if [[ "$preset" == "reduced_dbs" ]]; then
    $(python /app/alphafold/run_alphafold.py --hhblits_binary_path=$hhblits_binary_path --hhsearch_binary_path=$hhsearch_binary_path --jackhmmer_binary_path=$jackhmmer_binary_path --kalign_binary_path=$kalign_binary_path --small_bfd_database_path=$small_bfd_database_path --mgnify_database_path=$mgnify_database_path --template_mmcif_dir=$template_mmcif_dir --obsolete_pdbs_path=$obsolete_pdbs_path --pdb70_database_path=$pdb70_database_path --uniref90_database_path=$uniref90_database_path --data_dir=$data_dir --output_dir=$output_dir --fasta_paths=$fasta_path --model_names=$model_names --max_template_date=$max_template_date --preset=$preset --benchmark=$benchmark --logtostderr)
else
    $(python /app/alphafold/run_alphafold.py --hhblits_binary_path=$hhblits_binary_path --hhsearch_binary_path=$hhsearch_binary_path --jackhmmer_binary_path=$jackhmmer_binary_path --kalign_binary_path=$kalign_binary_path --bfd_database_path=$bfd_database_path --mgnify_database_path=$mgnify_database_path --template_mmcif_dir=$template_mmcif_dir --obsolete_pdbs_path=$obsolete_pdbs_path --pdb70_database_path=$pdb70_database_path --uniclust30_database_path=$uniclust30_database_path --uniref90_database_path=$uniref90_database_path --data_dir=$data_dir --output_dir=$output_dir --fasta_paths=$fasta_path --model_names=$model_names --max_template_date=$max_template_date --preset=$preset --benchmark=$benchmark --logtostderr)
fi