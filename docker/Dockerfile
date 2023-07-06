# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

ARG CUDA=11.1.1
FROM nvidia/cuda:${CUDA}-cudnn8-runtime-ubuntu18.04
# FROM directive resets ARGS, so we specify again (the value is retained if
# previously set).
ARG CUDA

# Use bash to support string substitution.
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
        git \
        wget \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get autoremove -y \
    && apt-get clean

RUN wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
RUN mv bin/micromamba /usr/local/bin/

WORKDIR /app

RUN git clone https://github.com/oteri/alphafold_MSA.git
WORKDIR alphafold_MSA
RUN micromamba create --prefix /app/env/condaenv  -y -c conda-forge openmm=7.7.0 cudatoolkit pdbfixer pip python=3.8
RUN micromamba run -p /app/env/condaenv pip3  install -r requirements.txt --no-cache-dir
RUN micromamba run -p /app/env/condaenv pip3  install --upgrade --no-cache-dir  \
  jax==0.3.25  \
  jaxlib==0.3.25+cuda11.cudnn805  \
  -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

RUN micromamba install  -p /app/env/condaenv  -y   -c nvidia cuda-nvcc=12.1
RUN cp /app/env/condaenv/bin/ptxas /usr/local/bin/

RUN micromamba clean --all --force-pkgs-dirs --yes

RUN wget -q -P /app/alphafold_MSA/alphafold/common/ \
  https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt

