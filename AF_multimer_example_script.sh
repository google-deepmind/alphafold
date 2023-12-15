#!/bin/bash

DOWNLOAD_DIR=/nvme1/common/af_data
UNIPROT_TO_NCBI_PATH=/nvme1/common/uniprot_to_ncbi.pkl
AFM_DIR="path_to_folder_with_results"
FASTAS_DIR="path_to_folder_with_queries"

# CHECK THESE CAREFULLY!!!
output_dir=$AFM_DIR
target_ids="filename_of_your_query(no_extension)"
models_to_relax=all
use_precomputed_msas=False
gpu_devices="0"  # "all"
# Basename of pickled dictionary containing externally matched species.
# If you don't want to use externally matched species use empty string.
externally_matched_species_dict_basename="externally_matched_species_dict.pkl"
# Basename of pickled collection containing species to which 
# many-to-some pairing should be restricted
many_to_some_species_to_pair_basename="many_to_some_species_to_pair.pkl"
# Basename of pickled dictionary containing the confidence values of 
# the pairings given by externally_matched_species_dict_basename.
confidences_externally_matched_species_basename="confidences_externally_matched_species.pkl"
# Lowest acceptable confidence value for pairing sequences
min_confidence=0
# Only match orthologs to query chains in multimer mode.
match_only_orthologs=False
# Whether to stop after input features are created, but before models are run.
stop_at_etl=False
# Number of recycles
num_recycle=3
# Tolerance for early stopping during recycling.
recycle_early_stop_tolerance=0.5


for target_id in $target_ids
do
  fasta_path=$FASTAS_DIR/$target_id.fasta
  results_dir=$output_dir/$target_id
  
  mkdir -p $results_dir

    python docker/run_docker.py \
      --fasta_paths=$fasta_path \
      --output_dir=$output_dir \
      --max_template_date=1000-01-01 \#to not use templates use this date
      --model_preset=multimer \
      --use_precomputed_msas=$use_precomputed_msas \
      --models_to_relax=$models_to_relax \
      --data_dir=$DOWNLOAD_DIR \
      --uniprot_to_ncbi_path=$UNIPROT_TO_NCBI_PATH \
      --externally_matched_species_dict_basename="$externally_matched_species_dict_basename" \
      --many_to_some_species_to_pair_basename="$many_to_some_species_to_pair_basename" \
      --confidences_externally_matched_species_basename="$confidences_externally_matched_species_basename" \
      --min_confidence=$min_confidence \
      --match_only_orthologs=$match_only_orthologs \
      --stop_at_etl=$stop_at_etl \
      --num_recycle=$num_recycle \
      --recycle_early_stop_tolerance=$recycle_early_stop_tolerance \
      --gpu_devices=$gpu_devices \
  &> $results_dir/af_OX_plus_mnemonics.log
done
