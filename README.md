# AlphaFold

This package provides an implementation of the inference pipeline of AlphaFold
v2. Using MSA as only input.

AlphaFold-Multimer is not yet provided
## Installation and running your first prediction

You will need a machine running Linux, AlphaFold does not support other
operating systems. Full installation requires up to 10GB of disk space and a modern 
NVIDIA GPU (GPUs with more memory can predict larger protein structures).

Please follow these steps:

1.  Install [Docker](https://www.docker.com/).
    *   Install
        [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html)
        for GPU support.
    *   Setup running
        [Docker as a non-root user](https://docs.docker.com/engine/install/linux-postinstall/#manage-docker-as-a-non-root-user).

1.  Clone this repository and `cd` into it.

    ```bash
    git clone https://github.com/oteri/alphafold_serverless.git
    cd ./alphafold
    ```

1.  Download model parameters:

    *   Please use this commands to download the parameters

    ```bash
    mkdir params
    cd params
    wget "https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar"
    tar xf alphafold_params_2022-12-06.tar
    rm alphafold_params_2022-12-06.tar
    ```

1.  Check that AlphaFold will be able to use a GPU by running:

    ```bash
    docker run --rm --gpus all nvidia/cuda:11.0-base nvidia-smi
    ```

    The output of this command should show a list of your GPUs. If it doesn't,
    check if you followed all steps correctly when setting up the
    [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html)
    or take a look at the following
    [NVIDIA Docker issue](https://github.com/NVIDIA/nvidia-docker/issues/1447#issuecomment-801479573).


1.  Build the Docker image. A makefile is provided. So this is enough:

    ```bash
    make build
    ```

    If you encounter the following error:

    ```
    W: GPG error: https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64 InRelease: The following signatures couldn't be verified because the public key is not available: NO_PUBKEY A4B469963BF863CC
    E: The repository 'https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64 InRelease' is not signed.
    ```

    use the workaround described in
    https://github.com/deepmind/alphafold/issues/463#issuecomment-1124881779.

1.  Make sure that the output directory exists (the default is `/tmp/alphafold`)
    and that you have sufficient permissions to write into it.

1.  Use the Makefile to run AF2. **Attention**: for the moment the Makefile must be
edited to correctly run AF2
    ```bash
    make run
    ```
1.  Once the run is over, the output directory shall contain predicted
    structures of the target protein. Please check the documentation below for
    additional options and troubleshooting tips.

### Model parameters

While the AlphaFold code is licensed under the Apache 2.0 License, the AlphaFold
parameters and CASP15 prediction data are made available under the terms of the
CC BY 4.0 license. Please see the [Disclaimer](#license-and-disclaimer) below
for more detail.

The AlphaFold parameters are available from
https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar, and
are downloaded as part of the `scripts/download_all_data.sh` script. This script
will download parameters for:

*   5 models which were used during CASP14, and were extensively validated for
    structure prediction quality (see Jumper et al. 2021, Suppl. Methods 1.12
    for details).
*   5 pTM models, which were fine-tuned to produce pTM (predicted TM-score) and
    (PAE) predicted aligned error values alongside their structure predictions
    (see Jumper et al. 2021, Suppl. Methods 1.9.7 for details).
*   5 AlphaFold-Multimer models that produce pTM and PAE values alongside their
    structure predictions.

## Running AlphaFold

**The simplest way to run AlphaFold is using the provided Makefile.** 
For the moment only basic capabilities are implemented.

1.  Relaxation is disabled

### Running AlphaFold-Multimer
AlphaFold-Multimer is not yet supported


### AlphaFold output

The outputs will be saved in a subdirectory of the directory provided via the
`--output_dir` flag of `run_docker.py`.

The contents of each output directory is as follows:
```<target_name>/
    features.pkl
    selected_prediction.pdb
```
*   `features.pkl` – A `pickle` file containing the input feature NumPy arrays
    used by the models to produce the structures.
*   `selected_prediction.pdb` – A `pdb` file containing the computed structure.

The pLDDT confidence measure is stored in the B-factor field of the output PDB
files (although unlike a B-factor, higher pLDDT is better, so care must be taken
when using for tasks such as molecular replacement).

This code has been tested to match mean top-1 accuracy on a CASP14 test set with
pLDDT ranking over 5 model predictions (some CASP targets were run with earlier
versions of AlphaFold and some had manual interventions; see our forthcoming
publication for details). Some targets such as T1064 may also have high
individual run variance over random seeds.

## Community contributions

Colab notebooks provided by the community (please note that these notebooks may
vary from our full AlphaFold system and we did not validate their accuracy):

*   The
    [ColabFold AlphaFold2 notebook](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)
    by Martin Steinegger, Sergey Ovchinnikov and Milot Mirdita, which uses an
    API hosted at the Södinglab based on the MMseqs2 server
    [(Mirdita et al. 2019, Bioinformatics)](https://academic.oup.com/bioinformatics/article/35/16/2856/5280135)
    for the multiple sequence alignment creation.

## Acknowledgements

AlphaFold communicates with and/or references the following separate libraries
and packages:

*   [Abseil](https://github.com/abseil/abseil-py)
*   [Biopython](https://biopython.org)
*   [Chex](https://github.com/deepmind/chex)
*   [Colab](https://research.google.com/colaboratory/)
*   [Docker](https://www.docker.com)
*   [HH Suite](https://github.com/soedinglab/hh-suite)
*   [HMMER Suite](http://eddylab.org/software/hmmer)
*   [Haiku](https://github.com/deepmind/dm-haiku)
*   [Immutabledict](https://github.com/corenting/immutabledict)
*   [JAX](https://github.com/google/jax/)
*   [Kalign](https://msa.sbc.su.se/cgi-bin/msa.cgi)
*   [matplotlib](https://matplotlib.org/)
*   [ML Collections](https://github.com/google/ml_collections)
*   [NumPy](https://numpy.org)
*   [OpenMM](https://github.com/openmm/openmm)
*   [OpenStructure](https://openstructure.org)
*   [pandas](https://pandas.pydata.org/)
*   [pymol3d](https://github.com/avirshup/py3dmol)
*   [SciPy](https://scipy.org)
*   [Sonnet](https://github.com/deepmind/sonnet)
*   [TensorFlow](https://github.com/tensorflow/tensorflow)
*   [Tree](https://github.com/deepmind/tree)
*   [tqdm](https://github.com/tqdm/tqdm)

We thank all their contributors and maintainers!

