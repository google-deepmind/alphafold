#!/bin/zsh 

# Aim: install alphafold v2.3 dependencies

CUDA=12.2.2


# Set DEBIAN_FRONTEND to noninteractive to avoid interactive prompts
export DEBIAN_FRONTEND=noninteractive

# Preseed debconf to use the Europe/London timezone
echo "tzdata tzdata/Areas select Europe" | debconf-set-selections
echo "tzdata tzdata/Zones/Europe select London" | debconf-set-selections

# Install tzdata and any other necessary packages
apt-get update && apt-get install -y tzdata

sudo apt-get update 
sudo DEBIAN_FRONTEND=noninteractive
sudo apt-get install --no-install-recommends -y \
    build-essential \
    cmake \
    cuda-command-line-tools-$(cut -f1,2 -d- <<< ${CUDA//./-}) \
    git \
    hmmer \
    kalign \
    wget
sudo apt-get install tmux parallel nvtop aria2 -y
sudo rm -rf /var/lib/apt/lists/* 
sudo apt-get autoremove -y 
sudo apt-get clean


# Compile HHsuite from source.
WD=/tmp
sudo mkdir -p $WD/hh-suite
sudo git clone --branch v3.3.0 https://github.com/soedinglab/hh-suite.git $WD/hh-suite 
sudo mkdir $WD/hh-suite/build
sudo pushd $WD/hh-suite/build
sudo sudo cmake -DCMAKE_INSTALL_PREFIX=/opt/hhsuite ..
sudo sudo make -j 4 && sudo make install
sudo sudo ln -s /opt/hhsuite/bin/* /usr/bin
sudo popd
sudo rm -rf $WD/hh-suite


# Install conda packages 
source $HOME/.zshrc  
envname=alphafold
# create conda environment if it does not exist
if ! conda info --envs | grep -q $envname; then
    conda create -y -c conda-forge -n $envname python=3.11
fi
conda init && conda activate $envname | exit 1
conda install -y -c conda-forge -n $envname openmm pdbfixer
conda clean --all --force-pkgs-dirs --yes


# Add SETUID bit to the ldconfig binary so that non-root users can run it.
chmod u+s /sbin/ldconfig.real


# download stereo_chemical_props.txt 
BASE=$(dirname $(dirname $(realpath $0)))  # /workspaces/alphafold
wget -q -P $BASE/alphafold/common/ \
  https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt


# Install pip packages.
BASE=$(dirname $(dirname $(realpath $0)))  # /workspaces/alphafold
pip install -r $BASE/requirements.txt --no-cache-dir
pip install --upgrade --no-cache-dir jax==0.4.10 jaxlib==0.4.10+cuda12.cudnn88 \
    -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html