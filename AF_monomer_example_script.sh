#!/bin/bash

DOWNLOAD_DIR=/nvme1/common/af_data
UNIPROT_TO_NCBI_PATH=/nvme1/common/uniprot_to_ncbi.pkl
AFM_DIR="path_to_folder_with_results"
FASTAS_DIR="path_to_folder_with_queries"

# CHECK THESE CAREFULLY!!!
output_dir=$AFM_DIR
target_ids="filename_of_your_query(no_extension)"
models_to_relax=none
use_precomputed_msas=True
gpu_devices="0"  # "all"
# To use a custom MSA you also need use_precomputed_msas to be set to True and
# you need to have the 3 precomputed msas files (even if they are not used).
# The file custom_msa_out_path should be an a3m file with query in the first row
# and it should be in the same folder as the other 3 precomputed msas.
custom_msa_out_path="filename_custom_msa.a3m"

for target_id in $target_ids
do
  echo $target_id
  fasta_path=$FASTAS_DIR/$target_id.fasta
  results_dir=$output_dir/$target_id
  
  mkdir -p $results_dir

  python docker/run_docker.py \
    --fasta_paths=$fasta_path \
    --output_dir=$output_dir \
    --max_template_date=1000-01-01 \#to not use templates use this date
    --model_preset=monomer \
    --use_precomputed_msas=$use_precomputed_msas \
    --uniprot_to_ncbi_path=$UNIPROT_TO_NCBI_PATH \
    --models_to_relax=$models_to_relax \
    --data_dir=$DOWNLOAD_DIR \
    --gpu_devices=$gpu_devices \
    --custom_msa_out_path=$custom_msa_out_path \
    --num_recycle=3 \
  &> $results_dir/af_OX_plus_mnemonics.log
done
