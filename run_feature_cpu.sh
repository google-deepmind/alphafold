#!/bin/bash
# Description: AlphaFold non-docker version
# Author: Sanjay Kumar Srikakulam
# Modified by Bozitao Zhong
# Edited by Yinying Yao for AF-Multimer version

usage() {
        echo ""
        echo "Usage: $0 <OPTIONS>"
        echo "Required Parameters:"
        echo "-d <data_dir>                                 Path to directory of supporting data"
        echo "-P <pretrained_data_dir>                      Choose db_preset model configuration - no ensembling (full_dbs) or 8 model ensemblings (casp14) (default: 'full_dbs')"
        echo "-o <output_dir>                               Path to a directory that will store the results."
        echo "-e <num_ensemble>                             Override default num_ensemble."
        # edited by Yinying
        echo "-m <model_preset>                             Choose preset model configuration - the monomer model, the monomer model with extra ensembling, monomer model with pTM head, or multimer model"
        echo "-n <num_multimer_predictions_per_model>       How many predictions (each with a different random seed) will be generated per model"
        echo "-j <nproc>                                    How many processors (each with a different random seed) will be used in the feature construction"
        echo "-f <fasta_path>                               Path to a FASTA file containing one sequence"
        echo "-t <max_template_date>                        Maximum template release date to consider (ISO-8601 format - i.e. YYYY-MM-DD). Important if folding historical test sets"
        echo "Optional Parameters:"
        echo "-b <benchmark>                                Run multiple JAX model evaluations to obtain a timing that excludes the compilation time, which should be more indicative of the time required for inferencing many
    proteins (default: 'False')"
        echo "-g <use_gpu>                                  Enable NVIDIA runtime to run with GPUs (default: 'True')"
        echo "-a <gpu_devices>                              Comma separated list of devices to pass to 'CUDA_VISIBLE_DEVICES' (default: 'all')"
        echo "-p <db_preset>                                Choose db_preset model configuration - no ensembling (full_dbs) or 8 model ensemblings (casp14) (default: 'full_dbs')"
        echo ""
        exit 1
}

while getopts ":d:P:o:m:e:f:it:a:n:j:p:bg" x; do
        case "${x}" in
        d)
                data_dir=$OPTARG
        ;;
        o)
                output_dir=$OPTARG
        ;;
        m)
                model_preset=$OPTARG
        ;;
        e)
                num_ensemble=$OPTARG
        ;;
        f)
                fasta_path=$OPTARG
        ;;
        n)
                num_multimer_predictions_per_model=$OPTARG
        ;;
        j)
                nproc=$OPTARG
        ;;
        t)
                max_template_date=$OPTARG
        ;;
        b)
                benchmark=true
        ;;
        g)
                use_gpu=true
        ;;
        a)
                gpu_devices=$OPTARG
        ;;
        p)
                db_preset=$OPTARG
        ;;
        P)
                pretrained_data_dir=$OPTARG
        ;;
        *)
                echo Unknown argument! exiting ...
                usage
        ;;
        esac
done

# Parse input and set defaults
if [[ "$data_dir" == "" || "$output_dir" == "" || "$model_preset" == "" || "$fasta_path" == "" || "$pretrained_data_dir" == "" ]] ; then
    usage
fi

if [[ "$max_template_date" == "" ]] ; then
    max_template_date=None
fi

if [[ "$benchmark" == "" ]] ; then
    benchmark=false
fi

if [[ "$use_gpu" == "" ]] ; then
    use_gpu=false
fi

if [[ "$nproc" == "" ]] ; then
    nproc=8
fi

if [[ "$gpu_devices" == "" ]] ; then
    gpu_devices="all"
fi

if [[ "$num_multimer_predictions_per_model" == "" ]] ; then
    num_multimer_predictions_per_model=5
fi

if [[ "$db_preset" == "" ]] ; then
    db_preset="full_dbs"
fi

if [[ "$db_preset" != "full_dbs" && "$db_preset" != "casp14" ]] ; then
    echo "Unknown db_preset! Using default ('full_dbs')"
    db_preset="full_dbs"
fi

# This bash script looks for the run_feature_cpu.py script in its current working directory, if it does not exist then exits
alphafold_script="$(readlink -f $(dirname $0))/run_feature_cpu.py"  

if [ ! -f "$alphafold_script" ]; then
    echo "Alphafold python script $alphafold_script does not exist."
    exit 1
fi

# Export ENVIRONMENT variables and set CUDA devices for use
if [[ "$use_gpu" == true ]] ; then
    export CUDA_VISIBLE_DEVICES=0

    if [[ "$gpu_devices" ]] ; then
        export CUDA_VISIBLE_DEVICES=0
    fi
fi

export TF_FORCE_UNIFIED_MEMORY='1'
export XLA_PYTHON_CLIENT_MEM_FRACTION='4.0'

# Path and user config (change me if required)
bfd_database_path="$data_dir/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
mgnify_database_path="$data_dir/mgnify/mgy_clusters.fa"
template_mmcif_dir="$data_dir/pdb_mmcif/mmcif_files"
obsolete_pdbs_path="$data_dir/pdb_mmcif/obsolete.dat"

uniref30_database_path="$data_dir/uniref30_uc30/UniRef30_2022_02/UniRef30_2022_02"
uniref90_database_path="$data_dir/uniref90/uniref90.fasta"


# Binary path (change me if required)
hhblits_binary_path=$(which hhblits)
hhsearch_binary_path=$(which hhsearch)
jackhmmer_binary_path=$(which jackhmmer)
kalign_binary_path=$(which kalign)

# added by Yinying for AF multimer

if [[ "$model_preset" == "" ]] ; then
    model_preset="monomer"
fi

if [[ "$model_preset" == "monomer" || "$model_preset" == "monomer_casp14" || "$model_preset" == "monomer_ptm" ]] ; then
    pdb70_database_path="$data_dir/pdb70/pdb70"
    cmd="python $alphafold_script \
        --use_gpu_relax=$use_gpu \
        --hhblits_binary_path=$hhblits_binary_path \
        --hhsearch_binary_path=$hhsearch_binary_path \
        --jackhmmer_binary_path=$jackhmmer_binary_path \
        --kalign_binary_path=$kalign_binary_path \
        --bfd_database_path=$bfd_database_path \
        --mgnify_database_path=$mgnify_database_path \
        --template_mmcif_dir=$template_mmcif_dir \
        --obsolete_pdbs_path=$obsolete_pdbs_path \
        --pdb70_database_path=$pdb70_database_path \
        --uniref30_database_path=$uniref30_database_path \
        --uniref90_database_path=$uniref90_database_path \
        --data_dir=$pretrained_data_dir \
        --output_dir=$output_dir \
        --random_seed=$(($RANDOM*$RANDOM*$RANDOM)) \
        --fasta_paths=$fasta_path \
        --model_preset=$model_preset \
        --max_template_date=$max_template_date \
        --db_preset=$db_preset \
        --benchmark=$benchmark \
        --num_threads=$nproc \
        --logtostderr"
    echo "$cmd"
    eval "$cmd"
elif [[  "$model_preset" =~ "multimer" ]] ; then

    uniprot_database_path="$data_dir/uniprot/uniprot.fasta"
    pdb_seqres_database_path="$data_dir/pdb_seqres/pdb_seqres.txt"
    cmd="python $alphafold_script \
        --use_gpu_relax=$use_gpu \
        --hhblits_binary_path=$hhblits_binary_path \
        --hhsearch_binary_path=$hhsearch_binary_path \
        --jackhmmer_binary_path=$jackhmmer_binary_path \
        --kalign_binary_path=$kalign_binary_path \
        --bfd_database_path=$bfd_database_path \
        --mgnify_database_path=$mgnify_database_path \
        --template_mmcif_dir=$template_mmcif_dir \
        --obsolete_pdbs_path=$obsolete_pdbs_path \
        --uniref30_database_path=$uniref30_database_path \
        --pdb_seqres_database_path=$pdb_seqres_database_path \
        --uniprot_database_path=$uniprot_database_path \
        --uniref90_database_path=$uniref90_database_path \
        --data_dir=$pretrained_data_dir \
        --output_dir=$output_dir \
        --fasta_paths=$fasta_path \
        --random_seed=$(($RANDOM*$RANDOM*$RANDOM)) \
        --model_preset=$model_preset \
        --num_multimer_predictions_per_model=$num_multimer_predictions_per_model \
        --max_template_date=$max_template_date \
        --db_preset=$db_preset \
        --benchmark=$benchmark \
        --num_threads=$nproc \
        --logtostderr"

    echo "$cmd"
    eval "$cmd"
elif [[ "$model_preset" != "monomer" && "$model_preset" != "monomer_casp14" && "$model_preset" != "monomer_ptm" && ! "$model_preset" =~ "multimer" ]] ; then
    echo "Unknown model_preset! "
    usage
fi



