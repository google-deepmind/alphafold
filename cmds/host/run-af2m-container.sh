#!/bin/zsh

# Aim: Initialize the alphafold2-multimer container

# ------------------------------------------------------------------------------
# FUNCTION
# ------------------------------------------------------------------------------
# Define the usage function
function usage() {
    echo "Usage: $(basename $0) [OPTIONS]"
    echo "  --outdir, -o <path>    : path to the HOST output directory (default: $(dirname $(dirname $(realpath $0)))/out)"
    echo "  --data, -d <path>      : path to the HOST AF2M data directory (default: /mnt/bob/shared/alphafold)"
    echo "  --container, -c <name> : name of the container (default: af2mrun)"
    echo "  --gpus <gpuIndex>      : index of the GPU to use (default: no GPU), accepts numbers or 'all'"
    echo "  --help, -h             : display this help"
    exit 1
}

function setGPUIndex() {
    if [ -z "$gpuIndex" ]; then
        # echo "No GPU specified, running without GPU access"
        echo ""
    elif [ "$gpuIndex" = "all" ]; then
        echo "device=all"
    else
        # echo "Using GPU $gpuIndex"
        # gpuIndex="device=$gpuIndex"
        echo "device=$gpuIndex"
    fi
}

# ------------------------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------------------------
# default values
hostOUTDIR=$(dirname $(dirname $(realpath $0)))/out
containerName="af2mrun"
AF2MDATA=/mnt/bob/shared/alphafold
gpuIndex=""  # Empty by default, meaning no GPU access

# --------------------
# get command line arguments and flags
# --------------------
while [[ $# -gt 1 ]]
do
    key="$1"
    case $key in
        --outdir|-o)
            hostOUTDIR="$2"
            shift 2;;
        --data|-d)
            AF2MDATA="$2"
            shift 2;;
        --container|-c)
            containerName="$2"
            shift 2;;
        --gpus)
            gpuIndex="$2"
            shift 2;;
        --help|-h)
            usage
            shift # past argument
            exit 1
            ;;
        *)
            echo "Illegal option: $key"
            usage
            exit 1
        ;;
    esac
done

echo "----------------------------------------"
echo "Inputs:"
echo "  hostOUTDIR   : $hostOUTDIR"
echo "  AF2MDATA     : $AF2MDATA"
echo "  containerName: $containerName"
echo "  gpuIndex     : $gpuIndex"
echo "----------------------------------------"
echo

gpuArg="$(setGPUIndex)"

# Add --gpus option only if gpuIndex is not empty
if [ -z "$gpuIndex" ]; then
    docker run --name $containerName \
        --network none \
        -itd \
        -v $AF2MDATA:/mnt/data/alphafold \
        -v $hostOUTDIR:/home/vscode/out \
        chunan/alphafold2.3:base
else
    docker run --name $containerName \
        --gpus $gpuArg \
        --network none \
        -itd \
        -v $AF2MDATA:/mnt/data/alphafold \
        -v $hostOUTDIR:/home/vscode/out \
        chunan/alphafold2.3:base
fi

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
