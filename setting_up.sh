# for A100 server jpas
# run the following as root
wget https://mirrors.bfsu.edu.cn/anaconda/archive/Anaconda3-2021.05-Linux-x86_64.sh

# install anaconda3, and set conda path to /opt/anaconda3/
# switch conda repo to mirrors.bfsu.edu.cn
conda create -y --name alphafold python==3.8

conda activate alphafold
conda install -y -c conda-forge openmm==7.5.1 cudnn==8.2.1.32 cudatoolkit==11.4.2 pdbfixer==1.7
conda install -y -c bioconda hmmer==3.3.2 hhsuite==3.3.0 kalign2==2.04

pip install  absl-py==0.13.0 biopython==1.79 chex==0.0.7 dm-haiku==0.0.4 dm-tree==0.1.6 immutabledict==2.0.0 jax==0.2.14 ml-collections==0.1.0 numpy==1.19.5 scipy==1.7.0 tensorflow-cpu==2.5.0

# fix CUBLAS_STATUS_EXECUTION_FAILED error for agis slurm with multiple GPU devices in one node
# if connection error, switch  DNS setting to 223.5.5.5, 119.29.29.29
pip install --upgrade jax jaxlib==0.1.72+cuda111 -f https://storage.googleapis.com/jax-releases/jax_releases.html

# alphafold official repo, with seperated featuring and modeling procedures.
cd /software
sudo git clone https://github.com/deepmind/alphafold.git
cd /software/alphafold/
alphafold_path=`pwd`
cd $CONDA_PREFIX/lib/python3.8/site-packages/ && patch -p0 < $alphafold_path/docker/openmm.patch

# upload self-managed run_feature_cpu.sh | run_feature_cpu.py | run_alphafold.sh | run_alphafold.py to $alphafold_path

# sms.py will use this for job-finished notifying
pip install aliyun-python-sdk-core

# get chemical props
wget -q -P $alphafold_path/alphafold/common/ https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt

# cuda 11.4
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-ubuntu2004.pin
sudo mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600
wget https://developer.download.nvidia.com/compute/cuda/11.4.2/local_installers/cuda-repo-ubuntu2004-11-4-local_11.4.2-470.57.02-1_amd64.deb
sudo dpkg -i cuda-repo-ubuntu2004-11-4-local_11.4.2-470.57.02-1_amd64.deb
sudo apt-key add /var/cuda-repo-ubuntu2004-11-4-local/7fa2af80.pub
sudo apt-get update
sudo apt-get -y install cuda

# database settings to avoid permission errors
# remember to rename each db file exactly the same as what it is in run_feature_cpu.sh and run_alphafold.sh
sudo chmod -R 755 /mnt/db

# runtime settings
# edit nproc to 16
cd $alphafold_path/alphafold/data/tools/
sudo vim hhblits.py hmmsearch.py jackhmmer.py


