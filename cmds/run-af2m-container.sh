#!/bin/zsh

# Aim: Initialize the alphafold2-multimer container

AF2MDATA=/mnt/bob/shared/alphafold

# change OUTDIR accordingly
HOSTOUTDIR=$(dirname $(dirname $(realpath $0)))/out
containerName="af2mrun"

# To use the second GPU (device=0 is the first GPU)
# replace all with '"device=1"' (including the quotes)
docker run --name $containerName \
    --gpus all \
    --network none \
    -itd \
    -v $AF2MDATA:/mnt/data/alphafold \
    -v $HOSTOUTDIR:/home/vscode/out \
    chunan/alphafold2.3:runtime

# # example exec for computing MSAs
# docker exec $containerName \
#     zsh /home/vscode/alphafold/cmds/run-af2m-msa.sh \
#     --fasta /home/vscode/out/agName/abName/abag.fasta \
#     --outdir /home/vscode/out/agName/abName/ \
#     --data /mnt/data/alphafold

# # example exec for predicting structures
# docker exec $containerName \
#     zsh /home/vscode/alphafold/cmds/run-af2m-struct.sh \
#     --fasta /home/vscode/out/agName/abName/abag.fasta \
#     --outdir /home/vscode/out/agName/abName/ \
#     --data /mnt/data/alphafold

# NOTE:
# The arguments for both steps are identical.
# The only difference is the script name.
# This is because predicting structure script makes use of the precomputed MSAs.
# And the location of precomputed MSAs are determined by the input FASTA file name and output directory.
