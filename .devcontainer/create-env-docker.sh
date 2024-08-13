#!/bin/zsh

# Aim: install alphafold v2.3 dependencies

set -e

CUDA=12.2.2

sudo locale-gen en_US.UTF-8

# Set DEBIAN_FRONTEND to noninteractive to avoid interactive prompts
export DEBIAN_FRONTEND=noninteractive

sudo apt-get update
export DEBIAN_FRONTEND=noninteractive
sudo apt-get install --no-install-recommends -y \
    build-essential \
    cmake \
    git \
    hmmer \
    kalign \
    wget \
    parallel \
    aria2


# Compile HHsuite from source.
WD=/tmp
sudo rm -rf $WD/hh-suite
sudo mkdir -p $WD/hh-suite
sudo git clone --branch v3.3.0 https://github.com/soedinglab/hh-suite.git $WD/hh-suite
sudo mkdir $WD/hh-suite/build
pushd $WD/hh-suite/build
sudo cmake -DCMAKE_INSTALL_PREFIX=/opt/hhsuite ..
sudo make -j 4 && sudo make install
sudo ln -s /opt/hhsuite/bin/* /usr/bin
popd
sudo rm -rf $WD/hh-suite


# Install conda packages
source /miniconda/etc/profile.d/conda.sh
envname=alphafold
if ! $(conda info --envs | grep -q $envname); then
    conda init
    conda create -y -c conda-forge -n $envname python=3.10  # do not use greater than 3.10, required by openmm
fi
conda activate $envname
conda install -y -c conda-forge -n $envname openmm=7.7.0 pdbfixer cudatoolkit


# Add SETUID bit to the ldconfig binary so that non-root users can run it.
sudo chmod u+s /sbin/ldconfig.real


# download stereo_chemical_props.txt
BASE=$(dirname $(dirname $(realpath $0)))  # /workspaces/alphafold
# only download if file not exist
[ ! -f $BASE/alphafold/common/stereo_chemical_props.txt ] && \
wget -q -P $BASE/alphafold/common/ \
  https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt


# Install pip packages.
BASE=$(dirname $(dirname $(realpath $0)))  # /workspaces/alphafold
pip install -r $BASE/requirements.txt --no-cache-dir
pip install --upgrade --no-cache-dir jax==0.4.10 jaxlib==0.4.10+cuda12.cudnn88 \
    -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

# downgrade ml_dtypes
pip install --no-cache-dir ml_dtypes==0.2.0

# install alphafold -e
pip install -e $BASE

# clean up
conda clean --all --force-pkgs-dirs --yes
pip cache purge
sudo apt autoremove -y
sudo rm -rf /var/lib/apt/lists/*
sudo apt-get clean