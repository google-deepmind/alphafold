# for A100 server jpas in ubuntu2004
# run the following as root

# user configuration

SOFTWARE_PATH=$(readlink -f /software)
DB_PATH='/mnt/db/'

# basic tools
if ! command -v aria2c &> /dev/null;then
  apt-get -y aria2
fi


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

# create the conda environment for alphafold from older versions
# supposed to be the name of the conda environment
conda create -n alphafold_2.3 --clone alphafold

conda activate alphafold_2.3

CONDA_PATH=$CONDA_PREFIX

# alphafold official repo, with seperated featuring and modeling procedures.
cd $SOFTWARE_PATH
git clone https://github.com/YaoYinYing/alphafold.git
cd alphafold
alphafold_path=$PWD
git checkout official-to-main

wget -q -P $alphafold_pat/alphafold/common/ \
  https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt

pip3 install --upgrade pip --no-cache-dir \
    && pip3 install -r $alphafold_path/requirements.txt --no-cache-dir \
    && pip3 install --upgrade --no-cache-dir \
      jax==0.3.25 \
      jaxlib==0.3.25+cuda11.cudnn805 \
      -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html


alphafold_path=$SOFTWARE_PATH/alphafold/
cd $CONDA_PREFIX/lib/python3.8/site-packages/ && patch -p0 < $alphafold_path/docker/openmm.patch

# sms.py will use this for job-finished notifying
pip install aliyun-python-sdk-core

# get chemical props
wget -q -P $alphafold_path/alphafold/common/ https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt


# database settings to avoid permission errors
# remember to rename each db file exactly the same as what it is in run_feature_cpu.sh and run_alphafold.sh


mkdir -p $DB_PATH/alphafold
pushd $DB_PATH/alphafold
# set the historical pretrained af parameters
awk '{
  split($0, arr, ",");
  if(arr[1]!="tag"){
    tag=arr[1];
    url=arr[2];
    short_name=arr[3];

    split(url,arr2,"/");
    filename=arr2[length(arr2)];
    print filename;

    # create directories
    system("mkdir "short_name);

    # downloading params
    system("pushd "short_name";if [[ -f "filename" && ! -f "filename".aria2c ]];then echo Find complete file "filename".; else aria2c -x 10 "url";fi; tar -xf "filename" ; rm -f "filename";parallel -k sha256sum {} ::: $(ls)  > "short_name".sha256; popd");

    }
  }'  $alphafold_path/pretrained_data_url.csv

bash $alphafold_path/scripts/update_pdb_mmcif.sh /mnt/db/
chmod -R 755 $DB_PATH/



