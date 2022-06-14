#!/bin/bash

# Set the path to download dir
DOWNLOAD_DIR=/datashare/alphafold
INPUT_DIR=~/scratch/fred862/code/alphafold/alphafold/input
REPO_DIR=~/scratch/fred862/code/alphafold

# Load your modules as before
# ON NARVAL USE cuda/11.4 !!!
#module load gcc/9.3.0 openmpi/4.0.3 cuda/11.2.2 cudnn/8.2.0
#module load kalign/2.03 hmmer/3.2.1 openmm-alphafold/7.5.1 hh-suite/3.3.0

#cd $SCRATCH # Set the appropriate folder where the repo is contained

# Generate your virtual environment in $SLURM_TMPDIR
#virtualenv --no-download ${SLURM_TMPDIR}/my_env
source ~/env/alphafold_env/bin/activate
pip show numpy

# Install alphafold and dependencies
#pip install --no-index pdbfixer==1.7 alphafold==2.1.1

# Run your commands
#python ${REPO_DIR}/run_alphafold.py \
python ${REPO_DIR}/run_alphafold.py \
   --data_dir=${DOWNLOAD_DIR} \
   --fasta_paths=${INPUT_DIR}/5VAY.fasta \
   --bfd_database_path=${DOWNLOAD_DIR}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
   --pdb70_database_path=${DOWNLOAD_DIR}/pdb70/pdb70 \
   --template_mmcif_dir=${DOWNLOAD_DIR}/pdb_mmcif/mmcif_files \
   --uniclust30_database_path=${DOWNLOAD_DIR}/uniclust30/uniclust30_2018_08/uniclust30_2018_08  \
   --uniref90_database_path=${DOWNLOAD_DIR}/uniref90/uniref90.fasta  \
   --hhblits_binary_path=${EBROOTHHMINSUITE}/bin/hhblits \
   --hhsearch_binary_path=${EBROOTHHMINSUITE}/bin/hhsearch \
   --jackhmmer_binary_path=${EBROOTHMMER}/bin/jackhmmer \
   --kalign_binary_path=${EBROOTKALIGN}/bin/kalign \
   --mgnify_database_path=${DOWNLOAD_DIR}/mgnify/mgy_clusters_2018_12.fa \
   --output_dir=${SCRATCH}/alphafold_output \
   --obsolete_pdbs_path=${DOWNLOAD_DIR}/pdb_mmcif/obsolete.dat \
   --max_template_date=2020-05-14 \
   --model_preset=monomer_casp14 \
   --use_gpu_relax=True
