#!/bin/zsh

# Aim: Create MSAs
# Input:
#   - FASTA: path to the fasta file
#   - OUTDIR: path to the output directory
#   - WAITPID: PID to wait for
# Output:
# Usage:
# Example:
# Dependencies:

# Function to wait for a process to complete
wait_for_pid() {
    local pid=$1
    while kill -0 "$pid" 2> /dev/null; do
        sleep 1
    done
}

# Check if a PID is provided as the first command-line argument
if [ -n "$1" ]; then
    # Wait for the specified PID to complete
    echo waiting for PID $1 to finish ...
    wait_for_pid $1
    echo PID $1 completed
    echo
fi


BASE=/workspaces/alphafold/
WD=$(dirname $(realpath $0))
OUTDIR=$WD
DATA=/mnt/bob/alphafold
TODAY=$(date +%Y-%m-%d)  # 2021-07-01
FASTA=$WD/abag.fasta
pushd $WD

python $BASE/run_alphafold_msa.py \
    --fasta_paths=$FASTA \
    --model_preset=multimer \
    --max_template_date=$TODAY \
    --output_dir=$OUTDIR \
    --data_dir=$DATA \
    --bfd_database_path=$DATA/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
    --uniref30_database_path=$DATA/uniref30/UniRef30_2021_03 \
    --uniref90_database_path=$DATA/uniref90/uniref90.fasta \
    --mgnify_database_path=$DATA/mgnify/mgy_clusters_2022_05.fa \
    --template_mmcif_dir=$DATA/pdb_mmcif/mmcif_files \
    --obsolete_pdbs_path=$DATA/pdb_mmcif/obsolete.dat \
    --pdb_seqres_database_path=$DATA/pdb_seqres/pdb_seqres.txt \
    --uniprot_database_path=$DATA/uniprot/uniprot.fasta \
    --use_gpu_relax=true  > $OUTDIR/run_msa.log 2>&1