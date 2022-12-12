# AlphaFold v2.3.0

This technical note describes updates in the code and model weights that were
made to produce AlphaFold v2.3.0 including updated training data.

We have fine-tuned new AlphaFold-Multimer weights using identical model
architecture but a new training cutoff of 2021-09-30. Previously released
versions of AlphaFold and AlphaFold-Multimer were trained using PDB structures
with a release date before 2018-04-30, a cutoff date chosen to coincide with the
start of the 2018 CASP13 assessment. The new training cutoff represents ~30%
more data to train AlphaFold and more importantly includes much more data on
large protein complexes. The new training cutoff includes 4× the number of
electron microscopy structures and in aggregate twice the number of large
structures (more than 2,000 residues)[^1]. Due to the significant increase in
the number of large structures, we are also able to increase the size of
training crops (subsets of the structure used to train AlphaFold) from 384 to
640 residues. These new AlphaFold-Multimer models are expected to be
substantially more accurate on large protein complexes even though we use the
same model architecture and training methodology as our previously released
AlphaFold-Multimer paper.

These models were initially developed in response to a request from the CASP
organizers to better understand baselines for the progress of structure
prediction in CASP15, and because of the significant increase in accuracy for
large targets, we are making them available as the default multimer models.
Since they were developed as baselines, we have emphasized minimal changes to
our previous AlphaFold-Multimer system while accommodating larger complexes.
In particular, we increase the number of chains used at training time from 8 to
20 and increase the maximum number of MSA sequences from 1,152 to 2,048 for 3 of
the 5 AlphaFold-Multimer models.

For the CASP15 baseline, we also used somewhat more expensive inference settings
that have been found externally to improve AlphaFold accuracy. We increase the
number of seeds per model to 20[^2] and increase the maximum number of
recyclings to 20 with early stopping[^3]. Increasing the number of seeds to 20
is recommended for very large or difficult targets but is not the default due to
increased computational time.

Overall, we expect these new models to be the preferred models whenever the
stoichiometry of the complex is known, including known monomeric structures. In
cases where the stoichiometry is unknown, such as in genome-scale prediction, it
is likely that single chain AlphaFold will be more accurate on average unless
the chain has several thousand residues.

The predicted structures used for the CASP15 baselines are available
[here](https://github.com/deepmind/alphafold/blob/main/docs/casp15_predictions.zip).


[^1]: wwPDB Consortium. "Protein Data Bank: the single global archive for 3D
  macromolecular structure data." Nucleic Acids Res. 47, D520–D528 (2018).

[^2]: Johansson-Åkhe, Isak, and Björn Wallner. "Improving peptide-protein
  docking with AlphaFold-Multimer using forced sampling." Frontiers in
  bioinformatics 2 (2022): 959160-959160.

[^3]: Gao, Mu, et al. "AF2Complex predicts direct physical interactions in
  multimeric proteins with deep learning." Nature communications 13.1 (2022):
  1-13.
