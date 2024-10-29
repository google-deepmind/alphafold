# JSON file format for AlphaFold Server jobs

You can
[download an example JSON file here](https://github.com/google-deepmind/alphafold/blob/main/server/example.json);
here we describe the contents of this example JSON file.

This JSON file consists of a list of dictionaries (even in the case of a single
dictionary, a single-element list must be used), with each dictionary containing
a job description. Therefore, you can specify multiple jobs in one JSON file.

Each job description contains a job name, a list of PRNG seeds (which can be an
empty list for automated random seed assignment), and a list of entities
(molecules) to be modeled.

AlphaFold Server JSON files are especially useful for automation of repetitive
modeling jobs (e.g. to screen interactions of one protein with a small number of
others). The easiest way to construct an initial JSON file is to run a modeling
job via AlphaFold Server GUI and use it as a template. AlphaFold Server will
produce a zip file containing modeling results. Inside the zip file you will
find a JSON file named `<job_name>_job_request.json` containing the job inputs.
These files offer a convenient starting point for generating new jobs as they
are easily editable in standard text editors or in programming environments like
Google Colab notebooks.

Note that comments are not allowed in JSON files.

## Job name, seeds and sequences

*   `name` is a string with the job name. This is how the job will appear as in
    the job history table.
*   `modelSeeds` is a list of strings of uint32 seed values (e.g.
    `["1593933729", "4273"]`). Seeds are used to run the modeling. We recommend
    providing an empty list, in which case a single random seed will be used.
    This is the recommended option.
*   `sequences` is a list of dictionaries that carry descriptions of the
    entities (molecules) for modeling.

```json
{
  "name": "Test Fold Job Number One",
  "modelSeeds": [],
  "sequences": [...]
}
```

## Entity types

Valid entity types mirror those available in the AlphaFold Server web interface:

*   `proteinChain` – used for proteins
*   `dnaSequence` – used for DNA (single strand)
*   `rnaSequence` – used for RNA (single strand)
*   `ligand` – used for allowed ligands
*   `ion` – used for allowed ions

### Protein chains

`sequence` is a string containing protein sequence; the same limitations as in
the UI are in place, e.g. only letters corresponding to amino acids are allowed,
as defined by IUPAC. Only 20 standard amino acid type are supported.

`count` is the number of copies of this protein chain (integer).

`glycans` is an optional list of dictionaries that carries descriptions of the
protein glycosylation.

*   `residues` is a string defining glycan. Please refer to the
    [FAQ](https://alphafoldserver.com/faq) for the format description and
    allowed glycans.
*   `position` is a position of the amino acid to which the glycan is attached
    (integer, 1-based indexing).

`modifications` is an optional list of dictionaries that carries descriptions of
the post-translational modifications.

*   `ptmType` is a string containing the
    [CCD code](https://www.wwpdb.org/data/ccd) of the modification; the same
    codes are allowed as in the UI.
*   `position` is a position of the modified amino acid (integer).
*   Allowed modifications: `CCD_SEP`, `CCD_TPO`, `CCD_PTR`, `CCD_NEP`,
    `CCD_HIP`, `CCD_ALY`, `CCD_MLY`, `CCD_M3L`, `CCD_MLZ`, `CCD_2MR`, `CCD_AGM`,
    `CCD_MCS`, `CCD_HYP`, `CCD_HY3`, `CCD_LYZ`, `CCD_AHB`, `CCD_P1L`, `CCD_SNN`,
    `CCD_SNC`, `CCD_TRF`, `CCD_KCR`, `CCD_CIR`, `CCD_YHA`

```json
{
  "proteinChain": {
    "sequence": "PREACHINGS",

    "glycans": [
      {
        "residues": "NAG(NAG)(BMA)",
        "position": 8
      },
      {
        "residues": "BMA",
        "position": 10
      }
    ],

    "modifications": [
      {
        "ptmType": "CCD_HY3",
        "ptmPosition": 1
      },
      {
        "ptmType": "CCD_P1L",
        "ptmPosition": 5
      }
    ],

    "count": 1
  }
},
{
  "proteinChain": {
    "sequence": "REACHER",
    "count": 1
  }
}
```

### DNA chains

Please note that the `dnaSequence` type refers to single stranded DNA. If you
wish to model double stranded DNA, please add a second `"dnaSequence`", carrying
the sequence of the reverse complement strand.

`sequence` is a string containing a DNA sequence; the same limitations as in the
UI are in place, i.e. only letters A, T, G, C are allowed.

`count` is a number of copies of this DNA chain (integer).

`modifications` is an optional list of dictionaries that carries descriptions of
the DNA chemical modifications.

*   `modificationType` is a string containing
    [CCD code](https://www.wwpdb.org/data/ccd) of modification; the same codes
    are allowed as in the UI.
*   `basePosition` is a position of the modified nucleotide (integer).
*   Allowed modifications: `CCD_5CM`, `CCD_C34`, `CCD_5HC`, `CCD_6OG`,
    `CCD_6MA`, `CCD_1CC`, `CCD_8OG`, `CCD_5FC`, `CCD_3DR`

```json
{
  "dnaSequence": {
    "sequence": "GATTACA",

    "modifications": [
      {
        "modificationType": "CCD_6OG",
        "basePosition": 1
      },
      {
        "modificationType": "CCD_6MA",
        "basePosition": 2
      }
    ],

    "count": 1
  }
},
{
  "dnaSequence": {
    "sequence": "TGTAATC",
    "count": 1
  }
}
```

### RNA chains

`sequence` is a string containing RNA sequence (single strand); the same
limitations as in the UI are in place, e.g. only letters A, U, G, C are allowed.

`count` is a number of copies of this RNA chain (integer).

`modifications` is an optional list of dictionaries that carries descriptions of
the RNA chemical modifications.

*   `modificationType` is a string containing
    [CCD code](https://www.wwpdb.org/data/ccd) of modification; the same codes
    are allowed as in the UI.
*   `basePosition` is a position of the modified nucleotide (integer).
*   Allowed modifications: `CCD_PSU`, `CCD_5MC`, `CCD_OMC`, `CCD_4OC`,
    `CCD_5MU`, `CCD_OMU`, `CCD_UR3`, `CCD_A2M`, `CCD_MA6`, `CCD_6MZ`, `CCD_2MG`,
    `CCD_OMG`, `CCD_7MG`, `CCD_RSQ`

```json
{
  "rnaSequence": {
    "sequence": "GUAC",

    "modifications": [
      {
        "modificationType": "CCD_2MG",
        "basePosition": 1
      },
      {
        "modificationType": "CCD_5MC",
        "basePosition": 4
      }
    ],

    "count": 1
  }
}
```

### Ligands

`ligand` is a string containing the [CCD code](https://www.wwpdb.org/data/ccd)
of the ligand; the same codes are allowed as in the UI.

`count` is the number of copies of this ligand (integer).

Allowed ligands: `CCD_ADP`, `CCD_ATP`, `CCD_AMP`, `CCD_GTP`, `CCD_GDP`,
`CCD_FAD`, `CCD_NAD`, `CCD_NAP`, `CCD_NDP`, `CCD_HEM`, `CCD_HEC`, `CCD_PLM`,
`CCD_OLA`, `CCD_MYR`, `CCD_CIT`, `CCD_CLA`, `CCD_CHL`, `CCD_BCL`, `CCD_BCB`

```json
{
  "ligand": {
    "ligand": "CCD_ATP",
    "count": 1
  }
},
{
  "ligand": {
    "ligand": "CCD_HEM",
    "count": 2
  }
}
```

### Ions

`ion` is a string containing [CCD code](https://www.wwpdb.org/data/ccd) of the
ion; the same codes are allowed as in the UI. The ion charge is implicitly
specified by the CCD code.

`count` is a number of copies of this ion (integer).

Allowed ions: `MG`, `ZN`, `CL`, `CA`, `NA`, `MN`, `K`, `FE`, `CU`, `CO`

```json
{
  "ion": {
    "ion": "MG",
    "count": 2
  }
},
{
  "ion": {
    "ion": "NA",
    "count": 3
  }
}
```

# Additional modeling jobs

You may specify multiple jobs in one JSON file. This is an example of a simple
job request for one protein chain and two copies of the palindromic DNA
sequence:

```json
{
  "name": "Test Fold Job Number Two",
  "modelSeeds": [],
  "sequences": [
    {
      "proteinChain": {
        "sequence": "TEACHINGS",
        "count": 1
      }
    },
    {
      "dnaSequence": {
        "sequence": "TAGCTA",
        "count": 2
      }
    }
  ]
}
```
