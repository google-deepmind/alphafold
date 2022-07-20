# for A100 server jpas in ubuntu2004
# run the following as root

# user configuration
CONDA_PATH="/opt/anaconda3"
SOFTWARE_PATH='/software'
DB_PATH='/mnt/db/'

# basic tools
apt-get -y aria2

mkdir -p $SOFTWARE_PATH

wget https://mirrors.bfsu.edu.cn/anaconda/archive/Anaconda3-2021.05-Linux-x86_64.sh

# install anaconda3, and set conda path to /opt/anaconda3/

bash Anaconda3-2021.05-Linux-x86_64.sh -bfp ${CONDA_PATH}

# setup conda repository
echo "channels:
  - defaults
  - https://USERNAME:PASSWORD@conda.graylab.jhu.edu
  - conda-forge
show_channel_urls: true
default_channels:
  - https://mirrors.bfsu.edu.cn/anaconda/pkgs/main
  - https://mirrors.bfsu.edu.cn/anaconda/pkgs/r
  - https://mirrors.bfsu.edu.cn/anaconda/pkgs/msys2
custom_channels:
  conda-forge: https://mirrors.bfsu.edu.cn/anaconda/cloud
  msys2: https://mirrors.bfsu.edu.cn/anaconda/cloud
  bioconda: https://mirrors.bfsu.edu.cn/anaconda/cloud
  menpo: https://mirrors.bfsu.edu.cn/anaconda/cloud
  pytorch: https://mirrors.bfsu.edu.cn/anaconda/cloud
  simpleitk: https://mirrors.bfsu.edu.cn/anaconda/cloud
report_errors: false
" > ~/.condarc
# refresh the conda cache
conda clean -i

# create the conda environment for alphafold
conda create -y --name alphafold python==3.8

conda activate alphafold
conda install -y -c conda-forge openmm==7.5.1 cudnn==8.2.1.32 cudatoolkit==11.4.2 pdbfixer==1.7
conda install -y -c bioconda hmmer==3.3.2 hhsuite==3.3.0 kalign2==2.04

pip install  absl-py==0.13.0 biopython==1.79 chex==0.0.7 dm-haiku==0.0.4 dm-tree==0.1.6 immutabledict==2.0.0 jax==0.2.14 ml-collections==0.1.0 numpy==1.19.5 scipy==1.7.0 tensorflow-cpu==2.5.0

# fix CUBLAS_STATUS_EXECUTION_FAILED error for agis slurm with multiple GPU devices in one node
# if connection error, switch  DNS setting to 223.5.5.5, 119.29.29.29
pip install --upgrade jax jaxlib==0.1.72+cuda111 -f https://storage.googleapis.com/jax-releases/jax_releases.html

# alphafold official repo, with seperated featuring and modeling procedures.
cd $SOFTWARE_PATH
git clone https://github.com/YaoYinYing/alphafold.git
cd alphafold
alphafold_path=$SOFTWARE_PATH/alphafold/
cd $CONDA_PREFIX/lib/python3.8/site-packages/ && patch -p0 < $alphafold_path/docker/openmm.patch

# upload self-managed run_feature_cpu.sh | run_feature_cpu.py | run_alphafold.sh | run_alphafold.py to $alphafold_path

# sms.py will use this for job-finished notifying
pip install aliyun-python-sdk-core

# get chemical props
wget -q -P $alphafold_path/alphafold/common/ https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt

# cuda 11.4
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-ubuntu2004.pin
mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600
wget https://developer.download.nvidia.com/compute/cuda/11.4.2/local_installers/cuda-repo-ubuntu2004-11-4-local_11.4.2-470.57.02-1_amd64.deb
dpkg -i cuda-repo-ubuntu2004-11-4-local_11.4.2-470.57.02-1_amd64.deb
apt-key add /var/cuda-repo-ubuntu2004-11-4-local/7fa2af80.pub
apt-get update
apt-get -y install cuda


# database settings to avoid permission errors
# remember to rename each db file exactly the same as what it is in run_feature_cpu.sh and run_alphafold.sh


mkdir -p $DB_PATH/alphafold

chmod -R 755 /mnt/db

# runtime settings
# edit nproc to 16
cd $alphafold_path/alphafold/data/tools/
vim hhblits.py hmmsearch.py jackhmmer.py


