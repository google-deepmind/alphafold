
#!/bin/bash

# use traditional way for conda environment
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate alphafold


#setting environment for cuda-11.0 gpu2-5
export PATH=/usr/local/cuda-11.4/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda-11.4/lib64:$LD_LIBRARY_PATH

usage() {
        echo ""
        echo "Usage: $0 <OPTIONS>"
        echo "Required Parameters:"
        # edited by Yinying
        echo "-m <model_preset>  Choose preset model configuration - the monomer model, the monomer model with extra ensembling, monomer model with pTM head, or multimer model"
        echo "-n <num_multimer_predictions_per_model>       How many predictions (each with a different random seed) will be generated per model"
        echo "-t <template_date> Maximum template release date to consider (ISO-8601 format - i.e. YYYY-MM-DD). Important if folding historical test sets"
        echo ""
        exit 1
}

while getopts ":m:t:n:e:" i; do
        case "${i}" in

        m)
                model_preset=$OPTARG
        ;;

        t)
                max_template_date=$OPTARG
        ;;
        n)
                num_multimer_predictions_per_model=$OPTARG
        ;;
        e)
                num_ensemble=$OPTARG
        ;;
        *)
                echo Unknown argument!
                usage
        ;;

        esac
done

if [[ "$max_template_date" == "" ]] ; then
    template_date=2021-10-30
fi

if [[ "$max_template_date" == "no" ]] ; then
    template_date=1900-01-01
fi


# edited by Yinying
if [[ "$model_preset" == "" ]] ; then
    model_preset="monomer"
fi

# set default num_ensemble
if [[ "$model_preset" == "monomer" || "$model_preset" == "monomer_ptm" || "$model_preset" == "multimer" ]] ; then
    if [[ "$num_ensemble" == "" ]] ; then
        num_ensemble=1
    fi
fi

if [[ "$model_preset" == "monomer_casp14"  ]] ; then
    if [[ "$num_ensemble" == "" ]] ; then
        num_ensemble=8
    fi
fi

if [[ "$model_preset" == "multimer" ]] ; then
    if [[ "$num_multimer_predictions_per_model" == "" ]];then
        num_multimer_predictions_per_model=5
    fi
fi



if [[ "$model_preset" != "monomer" && "$model_preset" != "monomer_casp14" && "$model_preset" != "monomer_ptm" && "$model_preset" != "multimer" ]] ; then
    echo "Unknown model_preset! "
    usage
fi


echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"



# go to alphafold run_feature process
dir=`pwd`;
#run_docker_pth=/mnt/data/alphafold/docker;

af_official_repo=$(readlink -f $(dirname $0)) ;
#af_official_repo=/software/alphafold_multimer/alphafold ;
out_dir=$dir/output;
res_dir=$dir/res;
db_dir=/mnt/db;


mkdir $res_dir;
mkdir $out_dir;
mkdir $res_dir/models
mkdir $res_dir/lite;
mkdir $res_dir/full;
mkdir $dir/processed;


AF_process(){
	local dir=$1;
	local i=$2;
	local decoy_name=${i%.fasta};

	cd $af_official_repo; # fix error caused by some path configs.
    if [ ! -f $res_dir/lite/$decoy_name\_AF2_lite.tar.bz2 ]; then
		# not started yet??
		if [ ! -f $out_dir/$decoy_name/features.pkl ]; then
            echo File does not exist: $out_dir/$decoy_name/features.pkl;
        	echo Modeling is not started: $i;
	        cmd="bash $af_official_repo/run_feature_cpu.sh \
                -d $db_dir \
                -o $out_dir \
                -m $model_preset \
                -n $num_multimer_predictions_per_model \
                -f $dir/$i \
                -t $template_date";
	        echo "$cmd";eval "$cmd"
        else
            echo Find feature files in $out_dir/$decoy_name/features.pkl;
            echo Skip the run_feature process: $decoy_name
        fi
        # featuring is not started yet, Skip this sequence.
		if [ ! -f $out_dir/$decoy_name/features.pkl ]; then
            echo File does not exist: $out_dir/$decoy_name/features.pkl;
        	echo Modeling is not started: $i;
        	echo Skip Modeling: $i;
	    else
	        #
            echo Find feature files in $out_dir/$decoy_name/features.pkl;
            echo Runing modeling process : $decoy_name
            cmd="bash $af_official_repo/run_alphafold.sh \
                    -d $db_dir \
                    -o $out_dir \
                    -m $model_preset \
                    -f $dir/$i \
                    -n $num_multimer_predictions_per_model \
                    -t $template_date \
                    -e $num_ensemble";

		    echo "$cmd";eval "$cmd"

		    cd $out_dir && \
		    echo Selecting best_model ...
		    pushd $decoy_name
		        for f in ranked*.pdb;do cp ${f} ${res_dir}/models/${decoy_name}_${f};done
		    popd
		    #cp $decoy_name/ranked_0.pdb $res_dir/best_model/${decoy_name}_ranked_0.pdb &&
		    echo Collecting results files .... && \
		    tar jcf $decoy_name\_AF2_lite.tar.bz2  --exclude *.pkl --exclude $decoy_name/msas $decoy_name && \
		    mv $decoy_name\_AF2_lite.tar.bz2 $res_dir/lite && \
		    tar jcf $decoy_name\_AF2_full.tar.bz2 --remove-files $decoy_name &&  \
		    mv $decoy_name\_AF2_full.tar.bz2 $res_dir/full && \
		    mv $dir/$i $dir/processed &
        fi
    else
            echo Find modeling files in $res_dir/lite/$decoy_name\_AF2_lite.tar.bz2;
            echo Skip the run_alphafold process: $decoy_name
    fi
}

cd $dir;
echo We are now in `pwd`

total=`ls  |grep .fasta |wc -l`;
fin=0;
rest=$total;

for i in `ls  |grep .fasta`;
do
    echo $dir/$i;
    cd $dir;
     # run processing ...
    AF_process $dir $i ;
    # mv $dir/$i $dir/processed;
    let fin++;
    let rest--;

done

wait;

echo Sending final notify ....
python $af_official_repo/sms.py `whoami` `basename $dir`  $total $fin $rest;

echo "++++++++++++++++++++++++++++++++++++++++"

echo "process end at : "
date
cd $dir;

