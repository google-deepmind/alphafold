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
        echo "No GPU specified, running without GPU access"
    elif [ "$gpuIndex" = "all" ]; then
        echo "Using all GPUs"
    else
        echo "Using GPU $gpuIndex"
        gpuIndex="device=$gpuIndex"
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

setGPUIndex

# Add --gpus option only if gpuIndex is not empty
docker run --name $containerName \
    ${gpuIndex:+--gpus $gpuIndex} \
    --network none \
    -itd \
    -v $AF2MDATA:/mnt/data/alphafold \
    -v $hostOUTDIR:/home/vscode/out \
    chunan/alphafold2.3:base