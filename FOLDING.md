# Folding your protein 

Help message:
```shell
Usage: /repo/alphafold/run_my_alphafold_multimer.sh <OPTIONS>
Required Parameters:
-m <model_preset>  Choose preset model configuration - the monomer model, the monomer model with extra ensembling, monomer model with pTM head, or multimer model
-n <num_multimer_predictions_per_model>       How many predictions (each with a different random seed) will be generated per model
-j <nproc>  How many processors (each with a different random seed) will be used in the feature construction. Default: 8
-t <template_date>  Maximum template release date to consider (ISO-8601 format - i.e. YYYY-MM-DD). Important if folding historical test sets. Default: 2023-03-15
-e <num_ensemble>  Ensemble Number for pre-inference
-p <pretrained_data_date> Pretrained data release date to consider (ISO-8601 format - i.e. YYYY-MM-DD). Important if folding historical test sets. Default: 2022-12-06
-r <run_relax>  Run relax to {best, all, none} model(s). Default: best
-c <clean_run>  Make a clean run, full results (massive pkls) will be deleted. Default: false
```


```shell
cd <working-directory>
# 0) for monomer with: 
# 1) default template date 
# 2) the latest pretrain weight release 
# 3) only relax the best decoy 
# 4) use 8 processors for msa building 
# 5) fewest ensemble number 
# 6) will not clean the full results
bash /repo/alphafold/run_my_alphafold_multimer.sh

# for monomer with no template
bash /repo/alphafold/run_my_alphafold_multimer.sh -t no

# for monomer with custom template date
bash /repo/alphafold/run_my_alphafold_multimer.sh -t 1994-01-16

# run w/ a customized pretrain dataset (2022-03-02) for historical reproducibility
bash /repo/alphafold/run_my_alphafold_multimer.sh -p 2022-03-02 

# for  the latest multimer_v3
bash /repo/alphafold/run_my_alphafold_multimer.sh -m multimer

# for deprecated multimer_v2 (released in 2022-03-02)
bash /repo/alphafold/run_my_alphafold_multimer.sh -m multimer_v2 -p 2022-03-02

# for deprecated multimer_v1 (released in 2021-10-27 and 2022-01-19)
bash /repo/alphafold/run_my_alphafold_multimer.sh -m multimer_v1 -p 2022-01-19

# relax all decoys
bash /repo/alphafold/run_my_alphafold_multimer.sh -r all

# don't relax any decoys
bash /repo/alphafold/run_my_alphafold_multimer.sh -r none

# run w/ customized ensemble numbers
bash /repo/alphafold/run_my_alphafold_multimer.sh -e 6

# run a clean mode (no pkls will be reserved!)
bash /repo/alphafold/run_my_alphafold_multimer.sh -c true

# run to a specific fasta file
bash /repo/alphafold/run_my_alphafold_snakemake.sh -f ./S4_nosig.fasta 

# run MSA building w/ customized number of processor
bash /repo/alphafold/run_my_alphafold_multimer.sh -j 32

```