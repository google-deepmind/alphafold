![header](imgs/header.jpg)

# AlphaFold

This package provides an implementation of the inference pipeline of AlphaFold
v2.0. This is a completely new model that was entered in CASP14 and published in
Nature. For simplicity, we refer to this model as AlphaFold throughout the rest
of this document.

We also provide an implementation of AlphaFold-Multimer. This represents a work
in progress and AlphaFold-Multimer isn't expected to be as stable as our monomer
AlphaFold system.
[Read the guide](#updating-existing-alphafold-installation-to-include-alphafold-multimers)
for how to upgrade and update code.

Any publication that discloses findings arising from using this source code or the model parameters should [cite](#citing-this-work) the
[AlphaFold  paper](https://doi.org/10.1038/s41586-021-03819-2) and, if
applicable, the [AlphaFold-Multimer paper](https://www.biorxiv.org/content/10.1101/2021.10.04.463034v1).

Please also refer to the
[Supplementary Information](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-03819-2/MediaObjects/41586_2021_3819_MOESM1_ESM.pdf)
for a detailed description of the method.

**You can use a slightly simplified version of AlphaFold with
[this Colab
notebook](https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb)**
or community-supported versions (see below).

If you have any questions, please contact the AlphaFold team at
[alphafold@deepmind.com](mailto:alphafold@deepmind.com).

![CASP14 predictions](imgs/casp14_predictions.gif)

## First time setup

You will need a machine running Linux, AlphaFold does not support other
operating systems.

The following steps are required in order to run AlphaFold:

1.  Install [Docker](https://www.docker.com/).
    *   Install
        [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html)
        for GPU support.
    *   Setup running
        [Docker as a non-root user](https://docs.docker.com/engine/install/linux-postinstall/#manage-docker-as-a-non-root-user).
1.  Download genetic databases (see below).
1.  Download model parameters (see below).
1.  Check that AlphaFold will be able to use a GPU by running:

    ```bash
    docker run --rm --gpus all nvidia/cuda:11.0-base nvidia-smi
    ```

    The output of this command should show a list of your GPUs. If it doesn't,
    check if you followed all steps correctly when setting up the
    [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html)
    or take a look at the following
    [NVIDIA Docker issue](https://github.com/NVIDIA/nvidia-docker/issues/1447#issuecomment-801479573).

If you wish to run AlphaFold using Singularity (a common containerization platform on HPC systems) we recommend using some of the
third party Singularity setups as linked in
https://github.com/deepmind/alphafold/issues/10 or
https://github.com/deepmind/alphafold/issues/24.

### Genetic databases

This step requires `aria2c` to be installed on your machine.

AlphaFold needs multiple genetic (sequence) databases to run:

*   [BFD](https://bfd.mmseqs.com/),
*   [MGnify](https://www.ebi.ac.uk/metagenomics/),
*   [PDB70](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/),
*   [PDB](https://www.rcsb.org/) (structures in the mmCIF format),
*   [PDB seqres](https://www.rcsb.org/) – only for AlphaFold-Multimer,
*   [Uniclust30](https://uniclust.mmseqs.com/),
*   [UniProt](https://www.uniprot.org/uniprot/) – only for AlphaFold-Multimer,
*   [UniRef90](https://www.uniprot.org/help/uniref).

We provide a script `scripts/download_all_data.sh` that can be used to download
and set up all of these databases:

*   Default:

    ```bash
    scripts/download_all_data.sh <DOWNLOAD_DIR>
    ```

    will download the full databases.

*   With `reduced_dbs`:

    ```bash
    scripts/download_all_data.sh <DOWNLOAD_DIR> reduced_dbs
    ```

    will download a reduced version of the databases to be used with the
    `reduced_dbs` database preset.

:ledger: **Note: The download directory `<DOWNLOAD_DIR>` should _not_ be a
subdirectory in the AlphaFold repository directory.** If it is, the Docker build
will be slow as the large databases will be copied during the image creation.

We don't provide exactly the database versions used in CASP14 – see the [note on
reproducibility](#note-on-reproducibility). Some of the databases are mirrored
for speed, see [mirrored databases](#mirrored-databases).

:ledger: **Note: The total download size for the full databases is around 415 GB
and the total size when unzipped is 2.2 TB. Please make sure you have a large
enough hard drive space, bandwidth and time to download. We recommend using an
SSD for better genetic search performance.**

The `download_all_data.sh` script will also download the model parameter files.
Once the script has finished, you should have the following directory structure:

```
$DOWNLOAD_DIR/                             # Total: ~ 2.2 TB (download: 438 GB)
    bfd/                                   # ~ 1.7 TB (download: 271.6 GB)
        # 6 files.
    mgnify/                                # ~ 64 GB (download: 32.9 GB)
        mgy_clusters_2018_12.fa
    params/                                # ~ 3.5 GB (download: 3.5 GB)
        # 5 CASP14 models,
        # 5 pTM models,
        # 5 AlphaFold-Multimer models,
        # LICENSE,
        # = 16 files.
    pdb70/                                 # ~ 56 GB (download: 19.5 GB)
        # 9 files.
    pdb_mmcif/                             # ~ 206 GB (download: 46 GB)
        mmcif_files/
            # About 180,000 .cif files.
        obsolete.dat
    pdb_seqres/                            # ~ 0.2 GB (download: 0.2 GB)
        pdb_seqres.txt
    small_bfd/                             # ~ 17 GB (download: 9.6 GB)
        bfd-first_non_consensus_sequences.fasta
    uniclust30/                            # ~ 86 GB (download: 24.9 GB)
        uniclust30_2018_08/
            # 13 files.
    uniprot/                               # ~ 98.3 GB (download: 49 GB)
        uniprot.fasta
    uniref90/                              # ~ 58 GB (download: 29.7 GB)
        uniref90.fasta
```

`bfd/` is only downloaded if you download the full databases, and `small_bfd/`
is only downloaded if you download the reduced databases.

### Model parameters

While the AlphaFold code is licensed under the Apache 2.0 License, the AlphaFold
parameters are made available under the terms of the CC BY 4.0 license. Please
see the [Disclaimer](#license-and-disclaimer) below for more detail.

The AlphaFold parameters are available from
https://storage.googleapis.com/alphafold/alphafold_params_2022-03-02.tar, and
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

### Updating existing AlphaFold installation to include AlphaFold-Multimers

If you have AlphaFold v2.0.0 or v2.0.1 you can either reinstall AlphaFold fully
from scratch (remove everything and run the setup from scratch) or you can do an
incremental update that will be significantly faster but will require a bit more
work. Make sure you follow these steps in the exact order they are listed below:

1.  **Update the code.**
    *   Go to the directory with the cloned AlphaFold repository and run
        `git fetch origin main` to get all code updates.
1.  **Download the UniProt and PDB seqres databases.**
    *   Run `scripts/download_uniprot.sh <DOWNLOAD_DIR>`.
    *   Remove `<DOWNLOAD_DIR>/pdb_mmcif`. It is needed to have PDB SeqRes and
        PDB from exactly the same date. Failure to do this step will result in
        potential errors when searching for templates when running
        AlphaFold-Multimer.
    *   Run `scripts/download_pdb_mmcif.sh <DOWNLOAD_DIR>`.
    *   Run `scripts/download_pdb_seqres.sh <DOWNLOAD_DIR>`.
1.  **Update the model parameters.**
    *   Remove the old model parameters in `<DOWNLOAD_DIR>/params`.
    *   Download new model parameters using
        `scripts/download_alphafold_params.sh <DOWNLOAD_DIR>`.
1.  **Follow [Running AlphaFold](#running-alphafold).**

#### API changes between v2.0.0 and v2.1.0

We tried to keep the API as much backwards compatible as possible, but we had to
change the following:

*   The `RunModel.predict()` now needs a `random_seed` argument as MSA sampling
    happens inside the Multimer model.
*   The `preset` flag in `run_alphafold.py` and `run_docker.py` was split into
    `db_preset` and `model_preset`.
*   The models to use are not specified using `model_names` but rather using the
    `model_preset` flag. If you want to customize which models are used for each
    preset, you will have to modify the the `MODEL_PRESETS` dictionary in
    `alphafold/model/config.py`.
*   Setting the `data_dir` flag is now needed when using `run_docker.py`.

#### API changes between v2.1.0 and v2.2.0

The AlphaFold-Multimer model weights have been updated, these new models have
greatly reduced numbers of clashes on average and are slightly more accurate.

A flag `--num_multimer_predictions_per_model` has been added that controls how
many predictions will be made per model, by default the offline system will run
each model 5 times for a total of 25 predictions.

The `--is_prokaryote_list` flag has been removed along with the `is_prokaryote`
argument in `run_alphafold.predict_structure()`, eukaryotes and prokaryotes are
now paired in the same way.

To use the deprecated v2.1.0 AlphaFold-Multimer model weights:

1.  Change `SOURCE_URL` in `scripts/download_alphafold_params.sh` to
`https://storage.googleapis.com/alphafold/alphafold_params_2022-01-19.tar`,
and download the old parameters.
2.  Remove the `_v2` in the multimer `MODEL_PRESETS` in `config.py`.

## Running AlphaFold

**The simplest way to run AlphaFold is using the provided Docker script.** This
was tested on Google Cloud with a machine using the `nvidia-gpu-cloud-image`
with 12 vCPUs, 85 GB of RAM, a 100 GB boot disk, the databases on an additional
3 TB disk, and an A100 GPU.

1.  Clone this repository and `cd` into it.

    ```bash
    git clone https://github.com/deepmind/alphafold.git
    ```

1.  Build the Docker image:

    ```bash
    docker build -f docker/Dockerfile -t alphafold .
    ```

    If you encounter the following error:

    ```
    W: GPG error: https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64 InRelease: The following signatures couldn't be verified because the public key is not available: NO_PUBKEY A4B469963BF863CC
    E: The repository 'https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64 InRelease' is not signed.
    ```

    use the workaround described in
    https://github.com/deepmind/alphafold/issues/463#issuecomment-1124881779.

1.  Install the `run_docker.py` dependencies. Note: You may optionally wish to
    create a
    [Python Virtual Environment](https://docs.python.org/3/tutorial/venv.html)
    to prevent conflicts with your system's Python environment.

    ```bash
    pip3 install -r docker/requirements.txt
    ```

1.  Make sure that the output directory exists (the default is `/tmp/alphafold`)
    and that you have sufficient permissions to write into it. You can make sure
    that is the case by manually running `mkdir /tmp/alphafold` and
    `chmod 770 /tmp/alphafold`.

1.  Run `run_docker.py` pointing to a FASTA file containing the protein
    sequence(s) for which you wish to predict the structure. If you are
    predicting the structure of a protein that is already in PDB and you wish to
    avoid using it as a template, then `max_template_date` must be set to be
    before the release date of the structure. You must also provide the path to
    the directory containing the downloaded databases. For example, for the
    T1050 CASP14 target:

    ```bash
    python3 docker/run_docker.py \
      --fasta_paths=T1050.fasta \
      --max_template_date=2020-05-14 \
      --data_dir=$DOWNLOAD_DIR
    ```

    By default, Alphafold will attempt to use all visible GPU devices. To use a
    subset, specify a comma-separated list of GPU UUID(s) or index(es) using the
    `--gpu_devices` flag. See
    [GPU enumeration](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/user-guide.html#gpu-enumeration)
    for more details.

1.  You can control which AlphaFold model to run by adding the
    `--model_preset=` flag. We provide the following models:

    * **monomer**: This is the original model used at CASP14 with no ensembling.

    * **monomer\_casp14**: This is the original model used at CASP14 with
      `num_ensemble=8`, matching our CASP14 configuration. This is largely
      provided for reproducibility as it is 8x more computationally
      expensive for limited accuracy gain (+0.1 average GDT gain on CASP14
      domains).

    * **monomer\_ptm**: This is the original CASP14 model fine tuned with the
      pTM head, providing a pairwise confidence measure. It is slightly less
      accurate than the normal monomer model.

    * **multimer**: This is the [AlphaFold-Multimer](#citing-this-work) model.
      To use this model, provide a multi-sequence FASTA file. In addition, the
      UniProt database should have been downloaded.

1.  You can control MSA speed/quality tradeoff by adding
    `--db_preset=reduced_dbs` or `--db_preset=full_dbs` to the run command. We
    provide the following presets:

    *   **reduced\_dbs**: This preset is optimized for speed and lower hardware
        requirements. It runs with a reduced version of the BFD database.
        It requires 8 CPU cores (vCPUs), 8 GB of RAM, and 600 GB of disk space.

    *   **full\_dbs**: This runs with all genetic databases used at CASP14.

    Running the command above with the `monomer` model preset and the
    `reduced_dbs` data preset would look like this:

    ```bash
    python3 docker/run_docker.py \
      --fasta_paths=T1050.fasta \
      --max_template_date=2020-05-14 \
      --model_preset=monomer \
      --db_preset=reduced_dbs \
      --data_dir=$DOWNLOAD_DIR
    ```

### Running AlphaFold-Multimer

All steps are the same as when running the monomer system, but you will have to

*   provide an input fasta with multiple sequences,
*   set `--model_preset=multimer`,

An example that folds a protein complex `multimer.fasta`:

```bash
python3 docker/run_docker.py \
  --fasta_paths=multimer.fasta \
  --max_template_date=2020-05-14 \
  --model_preset=multimer \
  --data_dir=$DOWNLOAD_DIR
```

By default the multimer system will run 5 seeds per model (25 total predictions)
for a small drop in accuracy you may wish to run a single seed per model.  This
can be done via the `--num_multimer_predictions_per_model` flag, e.g. set it to
`--num_multimer_predictions_per_model=1` to run a single seed per model.

### Examples

Below are examples on how to use AlphaFold in different scenarios.

#### Folding a monomer

Say we have a monomer with the sequence `<SEQUENCE>`. The input fasta should be:

```fasta
>sequence_name
<SEQUENCE>
```

Then run the following command:

```bash
python3 docker/run_docker.py \
  --fasta_paths=monomer.fasta \
  --max_template_date=2021-11-01 \
  --model_preset=monomer \
  --data_dir=$DOWNLOAD_DIR
```

#### Folding a homomer

Say we have a homomer with 3 copies of the same sequence
`<SEQUENCE>`. The input fasta should be:

```fasta
>sequence_1
<SEQUENCE>
>sequence_2
<SEQUENCE>
>sequence_3
<SEQUENCE>
```

Then run the following command:

```bash
python3 docker/run_docker.py \
  --fasta_paths=homomer.fasta \
  --max_template_date=2021-11-01 \
  --model_preset=multimer \
  --data_dir=$DOWNLOAD_DIR
```

#### Folding a heteromer

Say we have an A2B3 heteromer, i.e. with 2 copies of
`<SEQUENCE A>` and 3 copies of `<SEQUENCE B>`. The input fasta should be:

```fasta
>sequence_1
<SEQUENCE A>
>sequence_2
<SEQUENCE A>
>sequence_3
<SEQUENCE B>
>sequence_4
<SEQUENCE B>
>sequence_5
<SEQUENCE B>
```

Then run the following command:

```bash
python3 docker/run_docker.py \
  --fasta_paths=heteromer.fasta \
  --max_template_date=2021-11-01 \
  --model_preset=multimer \
  --data_dir=$DOWNLOAD_DIR
```

#### Folding multiple monomers one after another

Say we have a two monomers, `monomer1.fasta` and `monomer2.fasta`.

We can fold both sequentially by using the following command:

```bash
python3 docker/run_docker.py \
  --fasta_paths=monomer1.fasta,monomer2.fasta \
  --max_template_date=2021-11-01 \
  --model_preset=monomer \
  --data_dir=$DOWNLOAD_DIR
```

#### Folding multiple multimers one after another

Say we have a two multimers, `multimer1.fasta` and `multimer2.fasta`.

We can fold both sequentially by using the following command:

```bash
python3 docker/run_docker.py \
  --fasta_paths=multimer1.fasta,multimer2.fasta \
  --max_template_date=2021-11-01 \
  --model_preset=multimer \
  --data_dir=$DOWNLOAD_DIR
```

### AlphaFold output

The outputs will be saved in a subdirectory of the directory provided via the
`--output_dir` flag of `run_docker.py` (defaults to `/tmp/alphafold/`). The
outputs include the computed MSAs, unrelaxed structures, relaxed structures,
ranked structures, raw model outputs, prediction metadata, and section timings.
The `--output_dir` directory will have the following structure:

```
<target_name>/
    features.pkl
    ranked_{0,1,2,3,4}.pdb
    ranking_debug.json
    relaxed_model_{1,2,3,4,5}.pdb
    result_model_{1,2,3,4,5}.pkl
    timings.json
    unrelaxed_model_{1,2,3,4,5}.pdb
    msas/
        bfd_uniclust_hits.a3m
        mgnify_hits.sto
        uniref90_hits.sto
```

The contents of each output file are as follows:

*   `features.pkl` – A `pickle` file containing the input feature NumPy arrays
    used by the models to produce the structures.
*   `unrelaxed_model_*.pdb` – A PDB format text file containing the predicted
    structure, exactly as outputted by the model.
*   `relaxed_model_*.pdb` – A PDB format text file containing the predicted
    structure, after performing an Amber relaxation procedure on the unrelaxed
    structure prediction (see Jumper et al. 2021, Suppl. Methods 1.8.6 for
    details).
*   `ranked_*.pdb` – A PDB format text file containing the relaxed predicted
    structures, after reordering by model confidence. Here `ranked_0.pdb` should
    contain the prediction with the highest confidence, and `ranked_4.pdb` the
    prediction with the lowest confidence. To rank model confidence, we use
    predicted LDDT (pLDDT) scores (see Jumper et al. 2021, Suppl. Methods 1.9.6
    for details).
*   `ranking_debug.json` – A JSON format text file containing the pLDDT values
    used to perform the model ranking, and a mapping back to the original model
    names.
*   `timings.json` – A JSON format text file containing the times taken to run
    each section of the AlphaFold pipeline.
*   `msas/` - A directory containing the files describing the various genetic
    tool hits that were used to construct the input MSA.
*   `result_model_*.pkl` – A `pickle` file containing a nested dictionary of the
    various NumPy arrays directly produced by the model. In addition to the
    output of the structure module, this includes auxiliary outputs such as:

    *   Distograms (`distogram/logits` contains a NumPy array of shape [N_res,
        N_res, N_bins] and `distogram/bin_edges` contains the definition of the
        bins).
    *   Per-residue pLDDT scores (`plddt` contains a NumPy array of shape
        [N_res] with the range of possible values from `0` to `100`, where `100`
        means most confident). This can serve to identify sequence regions
        predicted with high confidence or as an overall per-target confidence
        score when averaged across residues.
    *   Present only if using pTM models: predicted TM-score (`ptm` field
        contains a scalar). As a predictor of a global superposition metric,
        this score is designed to also assess whether the model is confident in
        the overall domain packing.
    *   Present only if using pTM models: predicted pairwise aligned errors
        (`predicted_aligned_error` contains a NumPy array of shape [N_res,
        N_res] with the range of possible values from `0` to
        `max_predicted_aligned_error`, where `0` means most confident). This can
        serve for a visualisation of domain packing confidence within the
        structure.

The pLDDT confidence measure is stored in the B-factor field of the output PDB
files (although unlike a B-factor, higher pLDDT is better, so care must be taken
when using for tasks such as molecular replacement).

This code has been tested to match mean top-1 accuracy on a CASP14 test set with
pLDDT ranking over 5 model predictions (some CASP targets were run with earlier
versions of AlphaFold and some had manual interventions; see our forthcoming
publication for details). Some targets such as T1064 may also have high
individual run variance over random seeds.

## Inferencing many proteins

The provided inference script is optimized for predicting the structure of a
single protein, and it will compile the neural network to be specialized to
exactly the size of the sequence, MSA, and templates. For large proteins, the
compile time is a negligible fraction of the runtime, but it may become more
significant for small proteins or if the multi-sequence alignments are already
precomputed. In the bulk inference case, it may make sense to use our
`make_fixed_size` function to pad the inputs to a uniform size, thereby reducing
the number of compilations required.

We do not provide a bulk inference script, but it should be straightforward to
develop on top of the `RunModel.predict` method with a parallel system for
precomputing multi-sequence alignments. Alternatively, this script can be run
repeatedly with only moderate overhead.

## Note on CASP14 reproducibility

AlphaFold's output for a small number of proteins has high inter-run variance,
and may be affected by changes in the input data. The CASP14 target T1064 is a
notable example; the large number of SARS-CoV-2-related sequences recently
deposited changes its MSA significantly. This variability is somewhat mitigated
by the model selection process; running 5 models and taking the most confident.

To reproduce the results of our CASP14 system as closely as possible you must
use the same database versions we used in CASP. These may not match the default
versions downloaded by our scripts.

For genetics:

*   UniRef90:
    [v2020_01](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2020_01/uniref/)
*   MGnify:
    [v2018_12](http://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2018_12/)
*   Uniclust30: [v2018_08](http://wwwuser.gwdg.de/~compbiol/uniclust/2018_08/)
*   BFD: [only version available](https://bfd.mmseqs.com/)

For templates:

*   PDB: (downloaded 2020-05-14)
*   PDB70: [2020-05-13](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/old-releases/pdb70_from_mmcif_200513.tar.gz)

An alternative for templates is to use the latest PDB and PDB70, but pass the
flag `--max_template_date=2020-05-14`, which restricts templates only to
structures that were available at the start of CASP14.

## Citing this work

If you use the code or data in this package, please cite:

```bibtex
@Article{AlphaFold2021,
  author  = {Jumper, John and Evans, Richard and Pritzel, Alexander and Green, Tim and Figurnov, Michael and Ronneberger, Olaf and Tunyasuvunakool, Kathryn and Bates, Russ and {\v{Z}}{\'\i}dek, Augustin and Potapenko, Anna and Bridgland, Alex and Meyer, Clemens and Kohl, Simon A A and Ballard, Andrew J and Cowie, Andrew and Romera-Paredes, Bernardino and Nikolov, Stanislav and Jain, Rishub and Adler, Jonas and Back, Trevor and Petersen, Stig and Reiman, David and Clancy, Ellen and Zielinski, Michal and Steinegger, Martin and Pacholska, Michalina and Berghammer, Tamas and Bodenstein, Sebastian and Silver, David and Vinyals, Oriol and Senior, Andrew W and Kavukcuoglu, Koray and Kohli, Pushmeet and Hassabis, Demis},
  journal = {Nature},
  title   = {Highly accurate protein structure prediction with {AlphaFold}},
  year    = {2021},
  volume  = {596},
  number  = {7873},
  pages   = {583--589},
  doi     = {10.1038/s41586-021-03819-2}
}
```

In addition, if you use the AlphaFold-Multimer mode, please cite:

```bibtex
@article {AlphaFold-Multimer2021,
  author       = {Evans, Richard and O{\textquoteright}Neill, Michael and Pritzel, Alexander and Antropova, Natasha and Senior, Andrew and Green, Tim and {\v{Z}}{\'\i}dek, Augustin and Bates, Russ and Blackwell, Sam and Yim, Jason and Ronneberger, Olaf and Bodenstein, Sebastian and Zielinski, Michal and Bridgland, Alex and Potapenko, Anna and Cowie, Andrew and Tunyasuvunakool, Kathryn and Jain, Rishub and Clancy, Ellen and Kohli, Pushmeet and Jumper, John and Hassabis, Demis},
  journal      = {bioRxiv}
  title        = {Protein complex prediction with AlphaFold-Multimer},
  year         = {2021},
  elocation-id = {2021.10.04.463034},
  doi          = {10.1101/2021.10.04.463034},
  URL          = {https://www.biorxiv.org/content/early/2021/10/04/2021.10.04.463034},
  eprint       = {https://www.biorxiv.org/content/early/2021/10/04/2021.10.04.463034.full.pdf},
}
```

## Community contributions

Colab notebooks provided by the community (please note that these notebooks may
vary from our full AlphaFold system and we did not validate their accuracy):

*   The [ColabFold AlphaFold2 notebook](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)
    by Martin Steinegger, Sergey Ovchinnikov and Milot Mirdita, which uses an
    API hosted at the Södinglab based on the MMseqs2 server [(Mirdita et al.
    2019, Bioinformatics)](https://academic.oup.com/bioinformatics/article/35/16/2856/5280135)
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

## Get in Touch

If you have any questions not covered in this overview, please contact the
AlphaFold team at [alphafold@deepmind.com](mailto:alphafold@deepmind.com).

We would love to hear your feedback and understand how AlphaFold has been
useful in your research. Share your stories with us at
[alphafold@deepmind.com](mailto:alphafold@deepmind.com).

## License and Disclaimer

This is not an officially supported Google product.

Copyright 2021 DeepMind Technologies Limited.

### AlphaFold Code License

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this file except in compliance with the License. You may obtain a copy of the
License at https://www.apache.org/licenses/LICENSE-2.0.

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

### Model Parameters License

The AlphaFold parameters are made available under the terms of the Creative
Commons Attribution 4.0 International (CC BY 4.0) license. You can find details
at: https://creativecommons.org/licenses/by/4.0/legalcode

### Third-party software

Use of the third-party software, libraries or code referred to in the
[Acknowledgements](#acknowledgements) section above may be governed by separate
terms and conditions or license provisions. Your use of the third-party
software, libraries or code is subject to any such terms and you should check
that you can comply with any applicable restrictions or terms and conditions
before use.

### Mirrored Databases

The following databases have been mirrored by DeepMind, and are available with reference to the following:

*   [BFD](https://bfd.mmseqs.com/) (unmodified), by Steinegger M. and Söding J., available under a [Creative Commons Attribution-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-sa/4.0/).

*   [BFD](https://bfd.mmseqs.com/) (modified), by Steinegger M. and Söding J., modified by DeepMind, available under a [Creative Commons Attribution-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-sa/4.0/). See the Methods section of the [AlphaFold proteome paper](https://www.nature.com/articles/s41586-021-03828-1) for details.

*   [Uniclust30: v2018_08](http://wwwuser.gwdg.de/~compbiol/uniclust/2018_08/) (unmodified), by Mirdita M. et al., available under a [Creative Commons Attribution-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-sa/4.0/).

*   [MGnify: v2018_12](http://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/current_release/README.txt) (unmodified), by Mitchell AL et al., available free of all copyright restrictions and made fully and freely available for both non-commercial and commercial use under [CC0 1.0 Universal (CC0 1.0) Public Domain Dedication](https://creativecommons.org/publicdomain/zero/1.0/).
