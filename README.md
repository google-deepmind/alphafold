![header](imgs/header.jpg)

# AlphaFold

This package provides an implementation of the inference pipeline of AlphaFold
v2.0. This is a completely new model that was entered in CASP14 and published in
Nature. For simplicity, we refer to this model as AlphaFold throughout the rest
of this document.

Any publication that discloses findings arising from using this source code or
the model parameters should [cite](#citing-this-work) the
[AlphaFold paper](https://doi.org/10.1038/s41586-021-03819-2). Please also refer
to the
[Supplementary Information](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-03819-2/MediaObjects/41586_2021_3819_MOESM1_ESM.pdf)
for a detailed description of the method.

**You can use a slightly simplified version of AlphaFold with
[this Colab
notebook](https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb)**
or community-supported versions (see below).

![CASP14 predictions](imgs/casp14_predictions.gif)

# Modefied By WTTAT

Support AWS Batch.
Further use info on the way, no time now.

晚安67373.
See you next time ~

## Note on reproducibility

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
*   [pymol3d](https://github.com/avirshup/py3dmol)
*   [SciPy](https://scipy.org)
*   [Sonnet](https://github.com/deepmind/sonnet)
*   [TensorFlow](https://github.com/tensorflow/tensorflow)
*   [Tree](https://github.com/deepmind/tree)
*   [tqdm](https://github.com/tqdm/tqdm)

We thank all their contributors and maintainers!

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

The AlphaFold parameters are made available for non-commercial use only, under
the terms of the Creative Commons Attribution-NonCommercial 4.0 International
(CC BY-NC 4.0) license. You can find details at:
https://creativecommons.org/licenses/by-nc/4.0/legalcode

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
