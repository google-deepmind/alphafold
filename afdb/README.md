# AlphaFold Protein Structure Database

## Introduction

The AlphaFold UniProt release (214M predictions) is hosted on
[Google Cloud Public Datasets](https://console.cloud.google.com/marketplace/product/bigquery-public-data/deepmind-alphafold),
and is available to download at no cost under a
[CC-BY-4.0 licence](http://creativecommons.org/licenses/by/4.0/legalcode). The
dataset is in a Cloud Storage bucket, and metadata is available on BigQuery. A
Google Cloud account is required for the download, but the data can be freely
used under the terms of the
[CC-BY 4.0 Licence](http://creativecommons.org/licenses/by/4.0/legalcode).

This document provides an overview of how to access and download the dataset for
different use cases. Please refer to the [AlphaFold database FAQ](https://www.alphafold.com/faq)
for further information on what proteins are in the database and a changelog of
releases.

:ledger: **Note: The full dataset is difficult to manipulate without significant
computational resources (the size of the dataset is 23 TiB, 3 * 214M files).**

There are also alternatives to downloading the full dataset:

1.  Download a premade subset (covering important species / Swiss-Prot) via our
    [download page](https://alphafold.ebi.ac.uk/download).
2.  Download a custom subset of the data. See below.

If you need to download the full dataset then please see the "Bulk download"
section. See "Creating a Google Cloud Account" below for more information on how
to avoid any surprise costs when using Google Cloud Public Datasets.

## Licence

Data is available for academic and commercial use, under a
[CC-BY-4.0 licence](http://creativecommons.org/licenses/by/4.0/legalcode).

EMBL-EBI expects attribution (e.g. in publications, services or products) for
any of its online services, databases or software in accordance with good
scientific practice.

If you make use of an AlphaFold prediction, please cite the following papers:

*   [Jumper, J et al. Highly accurate protein structure prediction with
    AlphaFold. Nature
    (2021).](https://www.nature.com/articles/s41586-021-03819-2)
*   [Varadi, M et al. AlphaFold Protein Structure Database: massively expanding
    the structural coverage of protein-sequence space with high-accuracy models.
    Nucleic Acids Research
    (2021).](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkab1061/6430488)

AlphaFold Data Copyright (2022) DeepMind Technologies Limited.

## Disclaimer

The AlphaFold Data and other information provided on this site is for
theoretical modelling only, caution should be exercised in its use. It is
provided 'as-is' without any warranty of any kind, whether expressed or implied.
For clarity, no warranty is given that use of the information shall not infringe
the rights of any third party. The information is not intended to be a
substitute for professional medical advice, diagnosis, or treatment, and does
not constitute medical or other professional advice.

## Format

Dataset file names start with a protein identifier of the form `AF-[a UniProt
accession]-F[a fragment number]`.

Three files are provided for each entry:

*   **model_v4.cif** – contains the atomic coordinates for the predicted protein
    structure, along with some metadata. Useful references for this file format
    are the [ModelCIF](https://github.com/ihmwg/ModelCIF) and
    [PDBx/mmCIF](https://mmcif.wwpdb.org) project sites.
*   **confidence_v4.json** – contains a confidence metric output by AlphaFold
    called pLDDT. This provides a number for each residue, indicating how
    confident AlphaFold is in the *local* surrounding structure. pLDDT ranges
    from 0 to 100, where 100 is most confident. This is also contained in the
    CIF file.
*   **predicted_aligned_error_v4.json** – contains a confidence metric output by
    AlphaFold called PAE. This provides a number for every pair of residues,
    which is lower when AlphaFold is more confident in the relative position of
    the two residues. PAE is more suitable than pLDDT for judging confidence in 
    relative domain placements.
    [See here](https://alphafold.ebi.ac.uk/faq#faq-7) for a description of the
    format.

Predictions grouped by NCBI taxonomy ID are available as
`proteomes/proteome-tax_id-[TAX ID]-[SHARD ID]_v4.tar` within the same
bucket.

There are also two extra files stored in the bucket:

*   `accession_ids.csv` – This file contains a list of all the UniProt
    accessions that have predictions in AlphaFold DB. The file is in CSV format
    and includes the following columns, separated by a comma:
    *   UniProt accession, e.g. A8H2R3
    *   First residue index (UniProt numbering), e.g. 1
    *   Last residue index (UniProt numbering), e.g. 199
    *   AlphaFold DB identifier, e.g. AF-A8H2R3-F1
    *   Latest version, e.g. 4
*   `sequences.fasta` – This file contains sequences for all proteins in the
    current database version in FASTA format. The identifier rows start with
    ">AFDB", followed by the AlphaFold DB identifier and the name of the
    protein. The sequence rows contain the corresponding amino acid sequences.
    Each sequence is on a single line, i.e. there is no wrapping.

## Creating a Google Cloud Account

Downloading from the Google Cloud Public Datasets (rather than from AFDB or 3D
Beacons) requires a Google Cloud account. See the
[Google Cloud get started](https://cloud.google.com/docs/get-started) page, and
explore the [free tier account usage limits](https://cloud.google.com/free).

**IMPORTANT: After the trial period has finished (90 days), to continue access,
you are required to upgrade to a billing account. While your free tier access
(including access to the Public Datasets storage bucket) continues, usage beyond
the free tier will incur costs – please familiarise yourself with the pricing
for the services that you use to avoid any surprises.**

1.  Go to
    [https://cloud.google.com/datasets](https://cloud.google.com/datasets).
2.  Create an account:
    1.  Click "get started for free" in the top right corner.
    2.  Agree to all terms of service.
    3.  Follow the setup instructions. Note that a payment method is required,
        but this will not be used unless you enable billing.
    4.  Access to the Google Cloud Public Datasets storage bucket is always at
        no cost and you will have access to the
        [free tier.](https://cloud.google.com/free/docs/gcp-free-tier#free-tier-usage-limits)
3.  Set up a project:
    1.  In the top left corner, click the navigation menu (three horizontal bars
        icon).
    2.  Select: "Cloud overview" -> "Dashboard".
    3.  In the top left corner there is a project menu bar (likely says "My
        First Project"). Select this and a "Select a Project" box will appear.
    4.  To keep using this project, click "Cancel" at the bottom of the box.
    5.  To create a new project, click "New Project" at the top of the box:
        1.  Select a project name.
        2.  For location, if your organization has a Cloud account then select
            this, otherwise leave as is.
4.  Install `gsutil`:
    1.  Follow these
        [instructions](https://cloud.google.com/storage/docs/gsutil_install).

## Accessing the dataset

The data is available from:

*   GCS data bucket:
    [gs://public-datasets-deepmind-alphafold-v4](https://console.cloud.google.com/storage/browser/public-datasets-deepmind-alphafold-v4)

## Bulk download

We don't recommend downloading the full dataset unless required for processing
with local computational resources, for example in an academic high performance
computing centre.

We estimate that a 1 Gbps internet connection will allow download of the full
database in roughly 2.5 days.

While we don’t know the exact nature of your computational infrastructure, below
are some suggested approaches for downloading the dataset. Please reach out to
[alphafold@deepmind.com](mailto:alphafold@deepmind.com) if you have any
questions.

The recommended way of downloading the whole database is by downloading
1,015,797 sharded proteome tar files using the command below. This is
significantly faster than downloading all of the individual files because of
large constant per-file latency.

```bash
gsutil -m cp -r gs://public-datasets-deepmind-alphafold-v4/proteomes/ .
```

You will then have to un-tar all of the proteomes and un-gzip all of the
individual files. Note that after un-taring, there will be about 644M files, so
make sure your filesystem can handle this.

### Storage Transfer Service

Some users might find the
[Storage Transfer Service](https://cloud.google.com/storage-transfer-service) a
convenient way to set up the transfer between this bucket and another bucket, or
another cloud service. *Using this service may incur costs*. Please check the
[pricing page](https://cloud.google.com/storage-transfer/pricing) for more
detail, particularly for transfers to other cloud services.

## Downloading subsets of the data

### AlphaFold Database search

For simple queries, for example by protein name, gene name or UniProt accession
you can use the main search bar on
[alphafold.ebi.ac.uk](https://alphafold.ebi.ac.uk).

### 3D Beacons

[3D-Beacons](https://3d-beacons.org) is an international collaboration of
protein structure data providers to create a federated network with unified data
access mechanisms. The 3D-Beacons platform allows users to retrieve coordinate
files and metadata of experimentally determined and theoretical protein models
from data providers such as AlphaFold DB.

More information about how to access AlphaFold predictions using 3D-Beacons is
available at
[3D-Beacons documentation](https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/docs).

### Other premade species subsets

Downloads for some model organism proteomes, global health proteomes and
Swiss-Prot are available on the
[AFDB website](https://alphafold.ebi.ac.uk/download). These are generated from
[reference proteomes](https://www.uniprot.org/help/reference_proteome). If you
want other species, or *all* proteins for a particular species, please continue
reading.

We provide 1,015,797 sharded tar files for all species in
[gs://public-datasets-deepmind-alphafold-v4/proteomes/](https://console.cloud.google.com/storage/browser/public-datasets-deepmind-alphafold-v4/proteomes/).
We shard each proteome so that each shard contains at most 10,000 proteins
(which corresponds to 30,000 files per shard, since there are 3 files per
protein). To download a proteome of your choice, you have to do the following
steps:

1.  Find the [NCBI taxonomy ID](https://www.ncbi.nlm.nih.gov/taxonomy)
    (`[TAX_ID]`) of the species in question.
2.  Run `gsutil -m cp
    gs://public-datasets-deepmind-alphafold-v4/proteomes/proteome-tax_id-[TAX
    ID]-*_v4.tar .` to download all shards for this proteome.
3.  Un-tar all of the downloaded files and un-gzip all of the individual files.

### File manifests

Pre-made lists of files (manifests) are available at
[gs://public-datasets-deepmind-alphafold-v4/manifests](https://console.cloud.google.com/storage/browser/public-datasets-deepmind-alphafold-v4/manifests/).
Note that these filenames do not include the bucket prefix, but this can be
added once the files have been downloaded to your filesystem.

One can also define their own list of files, for example created by BigQuery
(see below). `gsutil` can be used to download these files with

```bash
cat [manifest file] | gsutil -m cp -I .
```

This will be much slower than downloading the tar files (grouped by species)
because each file has an associated overhead.

### BigQuery

**IMPORTANT: The
[free tier](https://cloud.google.com/bigquery/pricing#free-tier) of Google Cloud
comes with [BigQuery Sandbox](https://cloud.google.com/bigquery/docs/sandbox)
with 1 TB of free processed query data each month. Repeated queries within a
month could exceed this limit and if you have
[upgraded to a paid Cloud Billing account](https://cloud.google.com/free/docs/gcp-free-tier#how-to-upgrade)
you may be charged.**

**This should be sufficient for running a number of queries on the metadata
table, though the usage depends on the size of the columns queried and selected.
Please look at the
[BigQuery pricing page](https://cloud.google.com/bigquery/pricing) for more
information.**

**This is the user's responsibility so please ensure you keep track of your
billing settings and resource usage in the console.**

BigQuery provides a serverless and highly scalable analytics tool enabling SQL
queries over large datasets. The metadata for the UniProt dataset takes up ​​113
GiB and so can be challenging to process and analyse locally. The table name is:

*   BigQuery metadata table:
    [bigquery-public-data.deepmind_alphafold.metadata](https://console.cloud.google.com/bigquery?project=bigquery-public-data&ws=!1m5!1m4!4m3!1sbigquery-public-data!2sdeepmind_alphafold!3smetadata)

With BigQuery SQL you can do complex queries, e.g. find all high accuracy
predictions for a particular species, or even join on to other datasets, e.g. to
an experimental dataset by the `uniprotSequence`, or to the NCBI taxonomy by
`taxId`.

If you would find additional information in the metadata useful please file a
GitHub issue.

#### Setup

Follow the
[BigQuery Sandbox set up guide](https://cloud.google.com/bigquery/docs/sandbox).

#### Exploring the metadata

The column names and associated data types available can be found using the
following query.

```sql
SELECT column_name, data_type FROM bigquery-public-data.deepmind_alphafold.INFORMATION_SCHEMA.COLUMNS
WHERE table_name = 'metadata'
```

**Column name**        | **Data type**   | **Description**
---------------------- | --------------- | ---------------
allVersions            | `ARRAY<INT64>`  | An array of AFDB versions this prediction has had
entryId                | `STRING`        | The AFDB entry ID, e.g. "AF-Q1HGU3-F1"
fractionPlddtConfident | `FLOAT64`       | Fraction of the residues in the prediction with pLDDT between 70 and 90
fractionPlddtLow       | `FLOAT64`       | Fraction of the residues in the prediction with pLDDT between 50 and 70
fractionPlddtVeryHigh  | `FLOAT64`       | Fraction of the residues in the prediction with pLDDT greater than 90
fractionPlddtVeryLow   | `FLOAT64`       | Fraction of the residues in the prediction with pLDDT less than 50
gene                   | `STRING`        | The name of the gene if known, e.g. "COII"
geneSynonyms           | `ARRAY<STRING>` | Additional synonyms for the gene
globalMetricValue      | `FLOAT64`       | The mean pLDDT of this prediction
isReferenceProteome    | `BOOL`          | Is this protein part of the reference proteome?
isReviewed             | `BOOL`          | Has this protein been reviewed, i.e. is it part of SwissProt?
latestVersion          | `INT64`         | The latest AFDB version for this prediction
modelCreatedDate       | `DATE`          | The date of creation for this entry, e.g. "2022-06-01"
organismCommonNames    | `ARRAY<STRING>` | List of common organism names
organismScientificName | `STRING`        | The scientific name of the organism
organismSynonyms       | `ARRAY<STRING>` | List of synonyms for the organism
proteinFullNames       | `ARRAY<STRING>` | Full names of the protein
proteinShortNames      | `ARRAY<STRING>` | Short names of the protein
sequenceChecksum       | `STRING`        | [CRC64 hash](https://www.uniprot.org/help/checksum) of the sequence. Can be used for cheaper lookups.
sequenceVersionDate    | `DATE`          | Date when the sequence data was last modified in UniProt
taxId                  | `INT64`         | NCBI taxonomy id of the originating species
uniprotAccession       | `STRING`        | Uniprot accession ID
uniprotDescription     | `STRING`        | The name recommended by the UniProt consortium
uniprotEnd             | `INT64`         | Number of the last residue in the entry relative to the UniProt entry. This is equal to the length of the protein unless we are dealing with protein fragments.
uniprotId              | `STRING`        | The Uniprot EntryName field
uniprotSequence        | `STRING`        | Amino acid sequence for this prediction
uniprotStart           | `INT64`         | Number of the first residue in the entry relative to the UniProt entry. This is 1 unless we are dealing with protein fragments.

#### Producing summary statistics

The following query gives the mean of the prediction confidence fractions per
species.

```sql
SELECT
 organismScientificName AS name,
 SUM(fractionPlddtVeryLow) / COUNT(fractionPlddtVeryLow) AS mean_plddt_very_low,
 SUM(fractionPlddtLow) / COUNT(fractionPlddtLow) AS mean_plddt_low,
 SUM(fractionPlddtConfident) / COUNT(fractionPlddtConfident) AS mean_plddt_confident,
 SUM(fractionPlddtVeryHigh) / COUNT(fractionPlddtVeryHigh) AS mean_plddt_very_high,
 COUNT(organismScientificName) AS num_predictions
FROM bigquery-public-data.deepmind_alphafold.metadata
GROUP by name
ORDER BY num_predictions DESC;
```

#### Producing lists of files

We expect that the most important use for the metadata will be to create subsets
of proteins according to various criteria, so that users can choose to only copy
a subset of the 214M proteins that exist in the dataset. An example query is
given below:

```sql
with file_rows AS (
  with file_cols AS (
    SELECT
      CONCAT(entryID, '-model_v4.cif') as m,
      CONCAT(entryID, '-predicted_aligned_error_v4.json') as p
    FROM bigquery-public-data.deepmind_alphafold.metadata
    WHERE organismScientificName = "Homo sapiens"
      AND (fractionPlddtVeryHigh + fractionPlddtConfident) > 0.5
  )
  SELECT * FROM file_cols UNPIVOT (files for filetype in (m, p))
)
SELECT CONCAT('gs://public-datasets-deepmind-alphafold-v4/', files) as files
from file_rows
```

In this case, the list has been filtered to only include proteins from *Homo
sapiens* for which over half the residues are confident or better (>70 pLDDT).

This creates a table with one column "files", where each row is the cloud
location of one of the two file types that has been provided for each protein.
There is an additional `confidence_v4.json` file which contains the
per-residue pLDDT. This information is already in the CIF file but may be
preferred if only this information is required.

This allows users to bulk download the exact proteins they need, without having
to download the entire dataset. Other columns may also be used to select subsets
of proteins, and we point the user to the
[BigQuery documentation](https://cloud.google.com/bigquery/docs) to understand
other ways to filter for their desired protein lists. Likewise, the
documentation should be followed to download these file subsets locally, as the
most appropriate approach will depend on the filesize. Note that it may be
easier to download large files using [Colab](https://colab.research.google.com/)
(e.g. pandas to_csv).

#### Previous versions
Previous versions of AFDB will remain available at
[gs://public-datasets-deepmind-alphafold](https://console.cloud.google.com/storage/browser/public-datasets-deepmind-alphafold)
to enable reproducible research. We recommend using the latest version (v4).