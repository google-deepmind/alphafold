#!/bin/bash
ldconfig

# Author: WTTAT
# AWS Batch start script.
usage() {
    echo " "
    echo "something wrong"
    echo " "
    exit 1
}

# while getopts ":f:m:d:p:b:" i; do
while getopts ":f:m:d:p:h" i; do

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
        h)
                echo "usage: -f fasta -m model -d template -p preset corresponding to af2 parameters"
                exit 1

        # b)
        #         benchmark=$OPTARG
        # ;;
        esac
done

echo "BATCH_BUCKET : $BATCH_BUCKET"
echo "REGION : $REGION"
echo "fasta_paths : $fasta_paths"
echo "model_names : $model_names"
echo "max_template_date : $max_template_date"
echo "preset : $preset"
# echo "benchmark : $benchmark"

pwd

# Parse input and set defaults
if [[ "$model_names" == "" || "$fasta_paths" == "" || "$max_template_date" == "" ]] ; then
    usage
fi

if [[ "$benchmark" == "" ]] ; then
    benchmark=false
fi

if [[ "$preset" == "full" ]] ; then
    preset="full_dbs"
fi
if [[ "$preset" == "reduced" ]] ; then
    preset="reduced_dbs"
fi

if [[ "$preset" != "full_dbs" && "$preset" != "casp14" && "$preset" != "reduced_dbs" ]] ; then
    echo "Unknown preset! Using default ('full_dbs')"
    preset="full_dbs"
fi

echo "preset reset : $preset"
echo "benchmark reset : $benchmark"


# data_dir="/mnt/dataset/"
data_dir="/fsx/dataset"
# Path and user config (change me if required)

bfd_database_path="$data_dir/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
small_bfd_database_path="$data_dir/small_bfd/bfd-first_non_consensus_sequences.fasta"
uniclust30_database_path="$data_dir/uniclust30/uniclust30_2018_08/uniclust30_2018_08"

# mgnify_database_path="$data_dir/mgnify/mgy_clusters.fa"
# template_mmcif_dir="$data_dir/pdb_mmcif/mmcif_files"
# obsolete_pdbs_path="$data_dir/pdb_mmcif/obsolete.dat"
# pdb70_database_path="$data_dir/pdb70/pdb70"
# uniref90_database_path="$data_dir/uniref90/uniref90.fasta"

# Binary path (change me if required)
# hhblits_binary_path=$(which hhblits)
# hhsearch_binary_path=$(which hhsearch)
# jackhmmer_binary_path=$(which jackhmmer)
# kalign_binary_path=$(which kalign)


######
#  only support one file
#  by WTTAT
echo "start downloading"
aws s3 cp s3://$BATCH_BUCKET/$BATCH_DIR_PREFIX/$fasta_paths ./ --region $REGION

####
# get vCPU setting
vcpu=$[$(curl -s $ECS_CONTAINER_METADATA_URI | jq '.Limits.CPU')/1024]
echo "get vCPU : $vcpu"


# ######
echo "start running af2"
# Run AlphaFold with required parameters
# 'reduced_dbs' preset does not use bfd and uniclust30 databases
if [[ "$preset" == "reduced_dbs" ]]; then
    $(python /app/alphafold/run_alphafold.py --vcpu=$vcpu --BATCH_BUCKET="$BATCH_BUCKET" --small_bfd_database_path="$small_bfd_database_path" --fasta_paths="$fasta_paths" --model_names="$model_names" --max_template_date="$max_template_date" --preset="$preset" --benchmark="$benchmark" --logtostderr)
    # echo "running command : python /app/alphafold/run_alphafold.py --vcpu=$vcpu --BATCH_BUCKET="$BATCH_BUCKET" --small_bfd_database_path="$small_bfd_database_path" --fasta_paths="$fasta_paths" --model_names="$model_names" --max_template_date="$max_template_date" --preset="$preset" --benchmark="$benchmark" --logtostderr"
else
    $(python /app/alphafold/run_alphafold.py  --vcpu=$vcpu --BATCH_BUCKET="$BATCH_BUCKET"  --bfd_database_path="$bfd_database_path" --uniclust30_database_path="$uniclust30_database_path" --fasta_paths="$fasta_paths" --model_names="$model_names" --max_template_date="$max_template_date" --preset="$preset" --benchmark="$benchmark" --logtostderr)
    # echo  "running command : python /app/alphafold/run_alphafold.py  --vcpu=$vcpu --BATCH_BUCKET="$BATCH_BUCKET"  --bfd_database_path="$bfd_database_path" --uniclust30_database_path="$uniclust30_database_path" --fasta_paths="$fasta_paths" --model_names="$model_names" --max_template_date="$max_template_date" --preset="$preset" --benchmark="$benchmark" --logtostderr"
fi

echo "start ziping"

fasta_name=${fasta_paths%.*}

cd /app/output/
tar -zcvf $fasta_name.tar.gz $fasta_name/
# mv $fasta_name.tar.gz /app/output/$fasta_name/

echo "start uploading"
aws s3 sync /app/output/$fasta_name s3://$BATCH_BUCKET/output/$fasta_name  --region $REGION
# aws s3 cp /app/output/$fasta_name.tar.gz s3://$BATCH_BUCKET/output/  --region $REGION

# add metadata
aws s3 cp /app/output/$fasta_name.tar.gz s3://$BATCH_BUCKET/output/  --metadata {'"id"':'"'$id'"'} --region $REGION

echo "all done"