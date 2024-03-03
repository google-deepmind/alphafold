#!/bin/zsh

# Aim: Initialize the alphafold2-multimer container
# Input:
# Output:
# Usage:
# Example:
# Dependencies:

AF2MDATA=/mnt/bob/shared/alphafold

# replace all with '"device=1"' (including the quotes) to use the second GPU (device=0 is the first GPU)
# docker run --gpus all --rm \
#     -v $AF2MDATA:/mnt/data/alphafold \
#     -v $OUTDIR:/home/vscode/output \
#     $USER/chunan/alphafold2.3:runtime zsh

# change OUTDIR accordingly
HOSTOUTDIR=$(dirname $(dirname $(realpath $0)))/out
containerName="af2mrun"

docker run --name $containerName \
    --gpus all --network none \
    -itd \
    -v $AF2MDATA:/mnt/data/alphafold \
    -v $HOSTOUTDIR:/home/vscode/out \
    chunan/alphafold2.3:runtime

# example exec
docker exec $containerName \
    zsh /home/vscode/alphafold/cmds/run-af2m-msa.sh \
    --fasta /home/vscode/out/ADAM28/ILX510_57175103_57175103/abag.fasta \
    --outdir /home/vscode/out/ADAM28/ILX510_57175103_57175103/ \
    --data /mnt/data/alphafold
