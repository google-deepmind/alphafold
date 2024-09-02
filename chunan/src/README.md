# Commands

Folder structure

```plaintext
src/
├── host/
│   ├── build-af2-image.sh
│   └── run-af2m-container.sh
├── container/
│   ├── run-af2m-msa.sh
│   └── run-af2m-struct.sh
├── template/
│   ├── run_msa.sh
│   └── run_struct_pred.sh
└── README.md (Current file)
```

- `host`: has scripts run on the host machine.
- `container`: has scripts run inside the container (built into image).
- `template`: has template scripts for running MSA and structure prediction.

## ./host/build-af2-image.sh

This script builds the docker image for the AF2M container.

Execute the following command on the host machine to build the image:

```sh
defaultImageName="$USER/alphafold2.3:base"
./host/build-af2-image.sh $defaultImageName
```

Change the `defaultImageName` to the desired name of the image.

## ./host/run-af2m-container.sh

This script launches a container from the AlphaFold v2.3 image.

By default, the network access is disabled.
If you need to access the internet from the container, you need to
add the `--network host` option to the `docker run` command in the script.

Usage: (1) compute MSA; (2) Predict structure

```sh
wdDir="out"  # host working dir e.g. wd/{seqs,out} where seqs contains fasta files, out contains results
data="/mnt/data/alphafold"  # host alphafold data dir
containerName="af2mrun-msa"  # container name
./host/run-af2m-container.sh -o $wdDir -d $data -c $containerName
./host/run-af2m-container.sh -o $wdDir -d $data -c $containerName
```

- `wdDir` is mapped to `/home/vscode/out` in the container.

```sh
wdDir="out"  # host output dir
data="/mnt/data/alphafold"  # host alphafold data dir
containerName="af2mrun-struct"  # container name
./host/run-af2m-container.sh -o $wdDir -d $data -c $containerName
```

## ./container/run-af2m-msa.sh

This script is used to run the MSA prediction inside the container.

> !!! DO NOT EXECUTE THIS ON THE HOST MACHINE !!!

## ./container/run-af2m-struct.sh

This script is used to run the structure prediction inside the container.

> !!! DO NOT EXECUTE THIS ON THE HOST MACHINE !!!

## Run MSA and Structure Prediction

Execute the following commands to run MSA and structure prediction
using the containers launched by `host/run-af2m-container.sh`.

```sh
containerName="af2mrun-msa"  # container name

# example exec for computing MSAs
docker exec $containerName \
    zsh /home/vscode/alphafold/chunan/src/run-af2m-msa.sh \
    --fasta /home/vscode/out/agName/abName/abag.fasta \
    --outdir /home/vscode/out/agName/abName/ \
    --data /mnt/data/alphafold

# example exec for predicting structures
docker exec $containerName \
    zsh /home/vscode/alphafold/chunan/src/run-af2m-struct.sh \
    --fasta /home/vscode/out/agName/abName/abag.fasta \
    --outdir /home/vscode/out/agName/abName/ \
    --data /mnt/data/alphafold

```

NOTE:

The arguments for both steps are identical. The only difference is the script name.
This is because predicting structure script makes use of the precomputed MSAs.
And the location of precomputed MSAs are determined by the input FASTA file name and output directory.
