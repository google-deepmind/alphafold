# Commands

```
cmds
├── host
│   ├── build-af2-image.sh
│   └── run-af2m-container.sh
├── README.md
├── run-af2m-msa.sh
├── run-af2m-struct.sh
└── template
    ├── run_msa.sh
    └── run_struct_pred.sh
```

## host/build-af2-image.sh

This script builds the docker image for the AF2M container.

Execute the following command on the host machine to build the image:

```sh
defaultImageName="$USER/alphafold2.3:base"
./host/build-af2-image.sh $defaultImageName
```
change the `defaultImageName` to the desired name of the image.

## host/run-af2m-container.sh

This script launches a container from the AlphaFold v2.3 image.

By default, the network access is disabled.
If you need to access the internet from the container, you need to
add the `--network host` option to the `docker run` command in the script.

Usage: (1) compute MSA; (2) Predict structure

```sh
outdir="out"  # host output dir
data="/mnt/data/alphafold"  # host alphafold data dir
containerName="af2mrun-msa"  # container name
./cmds/run-af2m-container.sh -o $outdir -d $data -c $containerName
```

```sh
outdir="out"  # host output dir
data="/mnt/data/alphafold"  # host alphafold data dir
containerName="af2mrun-struct"  # container name
./cmds/run-af2m-container.sh -o $outdir -d $data -c $containerName
```

## cmds/run-af2m-msa.sh

This script is used to run the MSA prediction inside the container.

!!! DO NOT EXECUTE THIS ON THE HOST MACHINE !!!


## cmds/run-af2m-struct.sh

This script is used to run the structure prediction inside the container.

!!! DO NOT EXECUTE THIS ON THE HOST MACHINE !!!

## Run MSA and Structure Prediction

Execute the following commands to run MSA and structure prediction
using the containers launched by `host/run-af2m-container.sh`.

```sh
containerName="af2mrun-msa"  # container name

# example exec for computing MSAs
docker exec $containerName \
    zsh /home/vscode/alphafold/cmds/run-af2m-msa.sh \
    --fasta /home/vscode/out/agName/abName/abag.fasta \
    --outdir /home/vscode/out/agName/abName/ \
    --data /mnt/data/alphafold

# example exec for predicting structures
docker exec $containerName \
    zsh /home/vscode/alphafold/cmds/run-af2m-struct.sh \
    --fasta /home/vscode/out/agName/abName/abag.fasta \
    --outdir /home/vscode/out/agName/abName/ \
    --data /mnt/data/alphafold

```

NOTE:

The arguments for both steps are identical. The only difference is the script name.
This is because predicting structure script makes use of the precomputed MSAs.
And the location of precomputed MSAs are determined by the input FASTA file name and output directory.