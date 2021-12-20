#!/bin/bash
# Description: AlphaFold non-docker version
# Author: Sanjay Kumar Srikakulam
# Modified by Bozitao Zhong
# Edited by Yinying Yao for AF-Multimer version

usage() {
        echo ""
        echo "Usage: $0 <OPTIONS>"
        echo "Required Parameters:"
        echo "-d <data_dir>     Path to directory of supporting data"
        echo "-o <output_dir>   Path to a directory that will store the results."
        # edited by Yinying
        echo "-m <model_preset>  Choose preset model configuration - the monomer model, the monomer model with extra ensembling, monomer model with pTM head, or multimer model"
        echo "-f <fasta_path>   Path to a FASTA file containing one sequence"
        echo "-i is_prokaryote_list "
        echo "-t <max_template_date> Maximum template release date to consider (ISO-8601 format - i.e. YYYY-MM-DD). Important if folding historical test sets"
        echo "Optional Parameters:"
        echo "-b <benchmark>    Run multiple JAX model evaluations to obtain a timing that excludes the compilation time, which should be more indicative of the time required for inferencing many
    proteins (default: 'False')"
        echo "-g <use_gpu>      Enable NVIDIA runtime to run with GPUs (default: 'True')"
        echo "-a <gpu_devices>  Comma separated list of devices to pass to 'CUDA_VISIBLE_DEVICES' (default: 'all')"
        echo "-p <db_preset>       Choose db_preset model configuration - no ensembling (full_dbs) or 8 model ensemblings (casp14) (default: 'full_dbs')"
        echo ""
        exit 1
}

while getopts ":d:o:m:f:it:a:p:bg" x; do
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
        f)
                fasta_path=$OPTARG
        ;;
        i)
                is_prokaryote_list=true
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
        esac
done

# Parse input and set defaults
if [[ "$data_dir" == "" || "$output_dir" == "" || "$model_preset" == "" || "$fasta_path" == "" ]] ; then
    usage
fi

if [[ "$max_template_date" == "" ]] ; then
    max_template_date=None
fi

if [[ "$benchmark" == "" ]] ; then
    benchmark=false
fi

if [[ "$is_prokaryote_list" == "" ]] ; then
    is_prokaryote_list=false
fi

if [[ "$use_gpu" == "" ]] ; then
    use_gpu=true
fi

if [[ "$gpu_devices" == "" ]] ; then
    gpu_devices="all"
fi

if [[ "$db_preset" == "" ]] ; then
    db_preset="full_dbs"
fi

if [[ "$db_preset" != "full_dbs" && "$db_preset" != "casp14" ]] ; then
    echo "Unknown db_preset! Using default ('full_dbs')"
    db_preset="full_dbs"
fi


# This bash script looks for the run_alphafold_model.py script in its current working directory, if it does not exist then exits
alphafold_script="$(readlink -f $(dirname $0))/run_alphafold_model.py"

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

uniclust30_database_path="$data_dir/uniclust30/uniclust30_2018_08/uniclust30_2018_08"
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
    echo python $alphafold_script --hhblits_binary_path=$hhblits_binary_path --hhsearch_binary_path=$hhsearch_binary_path --jackhmmer_binary_path=$jackhmmer_binary_path --kalign_binary_path=$kalign_binary_path --bfd_database_path=$bfd_database_path --mgnify_database_path=$mgnify_database_path --template_mmcif_dir=$template_mmcif_dir --obsolete_pdbs_path=$obsolete_pdbs_path --pdb70_database_path=$pdb70_database_path --uniclust30_database_path=$uniclust30_database_path --uniref90_database_path=$uniref90_database_path --data_dir=$data_dir --output_dir=$output_dir --fasta_paths=$fasta_path --model_preset=$model_preset --max_template_date=$max_template_date --db_preset=$db_preset --benchmark=$benchmark --is_prokaryote_list=$is_prokaryote_list --logtostderr
    $(python $alphafold_script --hhblits_binary_path=$hhblits_binary_path --hhsearch_binary_path=$hhsearch_binary_path --jackhmmer_binary_path=$jackhmmer_binary_path --kalign_binary_path=$kalign_binary_path --bfd_database_path=$bfd_database_path --mgnify_database_path=$mgnify_database_path --template_mmcif_dir=$template_mmcif_dir --obsolete_pdbs_path=$obsolete_pdbs_path --pdb70_database_path=$pdb70_database_path --uniclust30_database_path=$uniclust30_database_path --uniref90_database_path=$uniref90_database_path --data_dir=$data_dir --output_dir=$output_dir --fasta_paths=$fasta_path --model_preset=$model_preset --max_template_date=$max_template_date --db_preset=$db_preset --benchmark=$benchmark --is_prokaryote_list=$is_prokaryote_list --logtostderr)
fi

if [[  "$model_preset" == "multimer" ]] ; then

    uniprot_database_path="$data_dir/uniprot/uniprot.fasta"
    pdb_seqres_database_path="$data_dir/pdb_seqres/pdb_seqres.txt"
    echo python $alphafold_script --hhblits_binary_path=$hhblits_binary_path --hhsearch_binary_path=$hhsearch_binary_path --jackhmmer_binary_path=$jackhmmer_binary_path --kalign_binary_path=$kalign_binary_path --bfd_database_path=$bfd_database_path --mgnify_database_path=$mgnify_database_path --template_mmcif_dir=$template_mmcif_dir --obsolete_pdbs_path=$obsolete_pdbs_path --uniclust30_database_path=$uniclust30_database_path --pdb_seqres_database_path=$pdb_seqres_database_path --uniprot_database_path=$uniprot_database_path --uniref90_database_path=$uniref90_database_path --data_dir=$data_dir --output_dir=$output_dir --fasta_paths=$fasta_path --model_preset=$model_preset --max_template_date=$max_template_date --db_preset=$db_preset --benchmark=$benchmark --is_prokaryote_list=$is_prokaryote_list --logtostderr
    $(python $alphafold_script --hhblits_binary_path=$hhblits_binary_path --hhsearch_binary_path=$hhsearch_binary_path --jackhmmer_binary_path=$jackhmmer_binary_path --kalign_binary_path=$kalign_binary_path --bfd_database_path=$bfd_database_path --mgnify_database_path=$mgnify_database_path --template_mmcif_dir=$template_mmcif_dir --obsolete_pdbs_path=$obsolete_pdbs_path --uniclust30_database_path=$uniclust30_database_path --pdb_seqres_database_path=$pdb_seqres_database_path --uniprot_database_path=$uniprot_database_path --uniref90_database_path=$uniref90_database_path --data_dir=$data_dir --output_dir=$output_dir --fasta_paths=$fasta_path --model_preset=$model_preset --max_template_date=$max_template_date --db_preset=$db_preset --benchmark=$benchmark --is_prokaryote_list=$is_prokaryote_list --logtostderr)

fi

if [[ "$model_preset" != "monomer" && "$model_preset" != "monomer_casp14" && "$model_preset" != "monomer_ptm" && "$model_preset" != "multimer" ]] ; then
    echo "Unknown model_preset! "
    usage
fi



