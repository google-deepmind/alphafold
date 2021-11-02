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

"""Pairing logic for multimer data pipeline."""

import collections
import functools
import re
import string
from typing import Any, Dict, Iterable, List, Sequence

from alphafold.common import residue_constants
from alphafold.data import pipeline
import numpy as np
import pandas as pd
import scipy.linalg

ALPHA_ACCESSION_ID_MAP = {x: y for y, x in enumerate(string.ascii_uppercase)}
ALPHANUM_ACCESSION_ID_MAP = {
    chr: num for num, chr in enumerate(string.ascii_uppercase + string.digits)
}  # A-Z,0-9
NUM_ACCESSION_ID_MAP = {str(x): x for x in range(10)}                # 0-9

MSA_GAP_IDX = residue_constants.restypes_with_x_and_gap.index('-')
SEQUENCE_GAP_CUTOFF = 0.5
SEQUENCE_SIMILARITY_CUTOFF = 0.9

MSA_PAD_VALUES = {'msa_all_seq': MSA_GAP_IDX,
                  'msa_mask_all_seq': 1,
                  'deletion_matrix_all_seq': 0,
                  'deletion_matrix_int_all_seq': 0,
                  'msa': MSA_GAP_IDX,
                  'msa_mask': 1,
                  'deletion_matrix': 0,
                  'deletion_matrix_int': 0}

MSA_FEATURES = ('msa', 'msa_mask', 'deletion_matrix', 'deletion_matrix_int')
SEQ_FEATURES = ('residue_index', 'aatype', 'all_atom_positions',
                'all_atom_mask', 'seq_mask', 'between_segment_residues',
                'has_alt_locations', 'has_hetatoms', 'asym_id', 'entity_id',
                'sym_id', 'entity_mask', 'deletion_mean',
                'prediction_atom_mask',
                'literature_positions', 'atom_indices_to_group_indices',
                'rigid_group_default_frame')
TEMPLATE_FEATURES = ('template_aatype', 'template_all_atom_positions',
                     'template_all_atom_mask')
CHAIN_FEATURES = ('num_alignments', 'seq_length')


domain_name_pattern = re.compile(
    r'''^(?P<pdb>[a-z\d]{4})
    \{(?P<bioassembly>[\d+(\+\d+)?])\}
    (?P<chain>[a-zA-Z\d]+)
    \{(?P<transform_index>\d+)\}$
    ''', re.VERBOSE)


def create_paired_features(
    chains: Iterable[pipeline.FeatureDict],
    prokaryotic: bool,
    ) ->  List[pipeline.FeatureDict]:
  """Returns the original chains with paired NUM_SEQ features.

  Args:
    chains:  A list of feature dictionaries for each chain.
    prokaryotic: Whether the target complex is from a prokaryotic organism.
      Used to determine the distance metric for pairing.

  Returns:
    A list of feature dictionaries with sequence features including only
    rows to be paired.
  """
  chains = list(chains)
  chain_keys = chains[0].keys()

  if len(chains) < 2:
    return chains
  else:
    updated_chains = []
    paired_chains_to_paired_row_indices = pair_sequences(
        chains, prokaryotic)
    paired_rows = reorder_paired_rows(
        paired_chains_to_paired_row_indices)

    for chain_num, chain in enumerate(chains):
      new_chain = {k: v for k, v in chain.items() if '_all_seq' not in k}
      for feature_name in chain_keys:
        if feature_name.endswith('_all_seq'):
          feats_padded = pad_features(chain[feature_name], feature_name)
          new_chain[feature_name] = feats_padded[paired_rows[:, chain_num]]
      new_chain['num_alignments_all_seq'] = np.asarray(
          len(paired_rows[:, chain_num]))
      updated_chains.append(new_chain)
    return updated_chains


def pad_features(feature: np.ndarray, feature_name: str) -> np.ndarray:
  """Add a 'padding' row at the end of the features list.

  The padding row will be selected as a 'paired' row in the case of partial
  alignment - for the chain that doesn't have paired alignment.

  Args:
    feature: The feature to be padded.
    feature_name: The name of the feature to be padded.

  Returns:
    The feature with an additional padding row.
  """
  assert feature.dtype != np.dtype(np.string_)
  if feature_name in ('msa_all_seq', 'msa_mask_all_seq',
                      'deletion_matrix_all_seq', 'deletion_matrix_int_all_seq'):
    num_res = feature.shape[1]
    padding = MSA_PAD_VALUES[feature_name] * np.ones([1, num_res],
                                                     feature.dtype)
  elif feature_name in ('msa_uniprot_accession_identifiers_all_seq',
                        'msa_species_identifiers_all_seq'):
    padding = [b'']
  else:
    return feature
  feats_padded = np.concatenate([feature, padding], axis=0)
  return feats_padded


def _make_msa_df(chain_features: pipeline.FeatureDict) -> pd.DataFrame:
  """Makes dataframe with msa features needed for msa pairing."""
  chain_msa = chain_features['msa_all_seq']
  query_seq = chain_msa[0]
  per_seq_similarity = np.sum(
      query_seq[None] == chain_msa, axis=-1) / float(len(query_seq))
  per_seq_gap = np.sum(chain_msa == 21, axis=-1) / float(len(query_seq))
  msa_df = pd.DataFrame({
      'msa_species_identifiers':
          chain_features['msa_species_identifiers_all_seq'],
      'msa_uniprot_accession_identifiers':
          chain_features['msa_uniprot_accession_identifiers_all_seq'],
      'msa_row':
          np.arange(len(
              chain_features['msa_uniprot_accession_identifiers_all_seq'])),
      'msa_similarity': per_seq_similarity,
      'gap': per_seq_gap
  })
  return msa_df


def _create_species_dict(msa_df: pd.DataFrame) -> Dict[bytes, pd.DataFrame]:
  """Creates mapping from species to msa dataframe of that species."""
  species_lookup = {}
  for species, species_df in msa_df.groupby('msa_species_identifiers'):
    species_lookup[species] = species_df
  return species_lookup


@functools.lru_cache(maxsize=65536)
def encode_accession(accession_id: str) -> int:
  """Map accession codes to the serial order in which they were assigned."""
  alpha = ALPHA_ACCESSION_ID_MAP        # A-Z
  alphanum = ALPHANUM_ACCESSION_ID_MAP  # A-Z,0-9
  num = NUM_ACCESSION_ID_MAP            # 0-9

  coding = 0

  # This is based on the uniprot accession id format
  # https://www.uniprot.org/help/accession_numbers
  if accession_id[0] in {'O', 'P', 'Q'}:
    bases = (alpha, num, alphanum, alphanum, alphanum, num)
  elif len(accession_id) == 6:
    bases = (alpha, num, alpha, alphanum, alphanum, num)
  elif len(accession_id) == 10:
    bases = (alpha, num, alpha, alphanum, alphanum, num, alpha, alphanum,
             alphanum, num)

  product = 1
  for place, base in zip(reversed(accession_id), reversed(bases)):
    coding += base[place] * product
    product *= len(base)

  return coding


def _calc_id_diff(id_a: bytes, id_b: bytes) -> int:
  return abs(encode_accession(id_a.decode()) - encode_accession(id_b.decode()))


def _find_all_accession_matches(accession_id_lists: List[List[bytes]],
                                diff_cutoff: int = 20
                                ) -> List[List[Any]]:
  """Finds accession id matches across the chains based on their difference."""
  all_accession_tuples = []
  current_tuple = []
  tokens_used_in_answer = set()

  def _matches_all_in_current_tuple(inp: bytes, diff_cutoff: int) -> bool:
    return all((_calc_id_diff(s, inp) < diff_cutoff for s in current_tuple))

  def _all_tokens_not_used_before() -> bool:
    return all((s not in tokens_used_in_answer for s in current_tuple))

  def dfs(level, accession_id, diff_cutoff=diff_cutoff) -> None:
    if level == len(accession_id_lists) - 1:
      if _all_tokens_not_used_before():
        all_accession_tuples.append(list(current_tuple))
        for s in current_tuple:
          tokens_used_in_answer.add(s)
      return

    if level == -1:
      new_list = accession_id_lists[level+1]
    else:
      new_list = [(_calc_id_diff(accession_id, s), s) for
                  s in accession_id_lists[level+1]]
      new_list = sorted(new_list)
      new_list = [s for d, s in new_list]

    for s in new_list:
      if (_matches_all_in_current_tuple(s, diff_cutoff) and
          s not in tokens_used_in_answer):
        current_tuple.append(s)
        dfs(level + 1, s)
        current_tuple.pop()
  dfs(-1, '')
  return all_accession_tuples


def _accession_row(msa_df: pd.DataFrame, accession_id: bytes) -> pd.Series:
  matched_df = msa_df[msa_df.msa_uniprot_accession_identifiers == accession_id]
  return matched_df.iloc[0]


def _match_rows_by_genetic_distance(
    this_species_msa_dfs: List[pd.DataFrame],
    cutoff: int = 20) -> List[List[int]]:
  """Finds MSA sequence pairings across chains within a genetic distance cutoff.

  The genetic distance between two sequences is approximated by taking the
  difference in their UniProt accession ids.

  Args:
    this_species_msa_dfs: a list of dataframes containing MSA features for
      sequences for a specific species. If species is missing for a chain, the
      dataframe is set to None.
    cutoff: the genetic distance cutoff.

  Returns:
    A list of lists, each containing M indices corresponding to paired MSA rows,
    where M is the number of chains.
  """
  num_examples = len(this_species_msa_dfs)  # N

  accession_id_lists = []  # M
  match_index_to_chain_index = {}
  for chain_index, species_df in enumerate(this_species_msa_dfs):
    if species_df is not None:
      accession_id_lists.append(
          list(species_df.msa_uniprot_accession_identifiers.values))
      # Keep track of which of the this_species_msa_dfs are not None.
      match_index_to_chain_index[len(accession_id_lists) - 1] = chain_index

  all_accession_id_matches = _find_all_accession_matches(
      accession_id_lists, cutoff)  # [k, M]

  all_paired_msa_rows = []  # [k, N]
  for accession_id_match in all_accession_id_matches:
    paired_msa_rows = []
    for match_index, accession_id in enumerate(accession_id_match):
      # Map back to chain index.
      chain_index = match_index_to_chain_index[match_index]
      seq_series = _accession_row(
          this_species_msa_dfs[chain_index], accession_id)

      if (seq_series.msa_similarity > SEQUENCE_SIMILARITY_CUTOFF or
          seq_series.gap > SEQUENCE_GAP_CUTOFF):
        continue
      else:
        paired_msa_rows.append(seq_series.msa_row)
    # If a sequence is skipped based on sequence similarity to the respective
    # target sequence or a gap cuttoff, the lengths of accession_id_match and
    # paired_msa_rows will be different. Skip this match.
    if len(paired_msa_rows) == len(accession_id_match):
      paired_and_non_paired_msa_rows = np.array([-1] * num_examples)
      matched_chain_indices = list(match_index_to_chain_index.values())
      paired_and_non_paired_msa_rows[matched_chain_indices] = paired_msa_rows
      all_paired_msa_rows.append(list(paired_and_non_paired_msa_rows))
  return all_paired_msa_rows


def _match_rows_by_sequence_similarity(this_species_msa_dfs: List[pd.DataFrame]
                                       ) -> List[List[int]]:
  """Finds MSA sequence pairings across chains based on sequence similarity.

  Each chain's MSA sequences are first sorted by their sequence similarity to
  their respective target sequence. The sequences are then paired, starting
  from the sequences most similar to their target sequence.

  Args:
    this_species_msa_dfs: a list of dataframes containing MSA features for
      sequences for a specific species.

  Returns:
   A list of lists, each containing M indices corresponding to paired MSA rows,
   where M is the number of chains.
  """
  all_paired_msa_rows = []

  num_seqs = [len(species_df) for species_df in this_species_msa_dfs
              if species_df is not None]
  take_num_seqs = np.min(num_seqs)

  sort_by_similarity = (
      lambda x: x.sort_values('msa_similarity', axis=0, ascending=False))

  for species_df in this_species_msa_dfs:
    if species_df is not None:
      species_df_sorted = sort_by_similarity(species_df)
      msa_rows = species_df_sorted.msa_row.iloc[:take_num_seqs].values
    else:
      msa_rows = [-1] * take_num_seqs  # take the last 'padding' row
    all_paired_msa_rows.append(msa_rows)
  all_paired_msa_rows = list(np.array(all_paired_msa_rows).transpose())
  return all_paired_msa_rows


def pair_sequences(examples: List[pipeline.FeatureDict],
                   prokaryotic: bool) -> Dict[int, np.ndarray]:
  """Returns indices for paired MSA sequences across chains."""

  num_examples = len(examples)

  all_chain_species_dict = []
  common_species = set()
  for chain_features in examples:
    msa_df = _make_msa_df(chain_features)
    species_dict = _create_species_dict(msa_df)
    all_chain_species_dict.append(species_dict)
    common_species.update(set(species_dict))

  common_species = sorted(common_species)
  common_species.remove(b'')  # Remove target sequence species.

  all_paired_msa_rows = [np.zeros(len(examples), int)]
  all_paired_msa_rows_dict = {k: [] for k in range(num_examples)}
  all_paired_msa_rows_dict[num_examples] = [np.zeros(len(examples), int)]

  for species in common_species:
    if not species:
      continue
    this_species_msa_dfs = []
    species_dfs_present = 0
    for species_dict in all_chain_species_dict:
      if species in species_dict:
        this_species_msa_dfs.append(species_dict[species])
        species_dfs_present += 1
      else:
        this_species_msa_dfs.append(None)

    # Skip species that are present in only one chain.
    if species_dfs_present <= 1:
      continue

    if np.any(
        np.array([len(species_df) for species_df in
                  this_species_msa_dfs if
                  isinstance(species_df, pd.DataFrame)]) > 600):
      continue

    # In prokaryotes (and some eukaryotes), interacting genes are often
    # co-located on the chromosome into operons. Because of that we can assume
    # that if two proteins' intergenic distance is less than a threshold, they
    # two proteins will form an an interacting pair.
    # In most eukaryotes, a single protein's MSA can contain many paralogs.
    # Two genes may interact even if they are not close by genomic distance.
    # In case of eukaryotes, some methods pair MSA sequences using sequence
    # similarity method.
    # See Jinbo Xu's work:
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6030867/#B28.
    if prokaryotic:
      paired_msa_rows = _match_rows_by_genetic_distance(this_species_msa_dfs)

      if not paired_msa_rows:
        continue
    else:
      paired_msa_rows = _match_rows_by_sequence_similarity(this_species_msa_dfs)
    all_paired_msa_rows.extend(paired_msa_rows)
    all_paired_msa_rows_dict[species_dfs_present].extend(paired_msa_rows)
  all_paired_msa_rows_dict = {
      num_examples: np.array(paired_msa_rows) for
      num_examples, paired_msa_rows in all_paired_msa_rows_dict.items()
  }
  return all_paired_msa_rows_dict


def reorder_paired_rows(all_paired_msa_rows_dict: Dict[int, np.ndarray]
                        ) -> np.ndarray:
  """Creates a list of indices of paired MSA rows across chains.

  Args:
    all_paired_msa_rows_dict: a mapping from the number of paired chains to the
      paired indices.

  Returns:
    a list of lists, each containing indices of paired MSA rows across chains.
    The paired-index lists are ordered by:
      1) the number of chains in the paired alignment, i.e, all-chain pairings
         will come first.
      2) e-values
  """
  all_paired_msa_rows = []

  for num_pairings in sorted(all_paired_msa_rows_dict, reverse=True):
    paired_rows = all_paired_msa_rows_dict[num_pairings]
    paired_rows_product = abs(np.array([np.prod(rows) for rows in paired_rows]))
    paired_rows_sort_index = np.argsort(paired_rows_product)
    all_paired_msa_rows.extend(paired_rows[paired_rows_sort_index])

  return np.array(all_paired_msa_rows)


def block_diag(*arrs: np.ndarray, pad_value: float = 0.0) -> np.ndarray:
  """Like scipy.linalg.block_diag but with an optional padding value."""
  ones_arrs = [np.ones_like(x) for x in arrs]
  off_diag_mask = 1.0 - scipy.linalg.block_diag(*ones_arrs)
  diag = scipy.linalg.block_diag(*arrs)
  diag += (off_diag_mask * pad_value).astype(diag.dtype)
  return diag


def _correct_post_merged_feats(
    np_example: pipeline.FeatureDict,
    np_chains_list: Sequence[pipeline.FeatureDict],
    pair_msa_sequences: bool) -> pipeline.FeatureDict:
  """Adds features that need to be computed/recomputed post merging."""

  np_example['seq_length'] = np.asarray(np_example['aatype'].shape[0],
                                        dtype=np.int32)
  np_example['num_alignments'] = np.asarray(np_example['msa'].shape[0],
                                            dtype=np.int32)

  if not pair_msa_sequences:
    # Generate a bias that is 1 for the first row of every block in the
    # block diagonal MSA - i.e. make sure the cluster stack always includes
    # the query sequences for each chain (since the first row is the query
    # sequence).
    cluster_bias_masks = []
    for chain in np_chains_list:
      mask = np.zeros(chain['msa'].shape[0])
      mask[0] = 1
      cluster_bias_masks.append(mask)
    np_example['cluster_bias_mask'] = np.concatenate(cluster_bias_masks)

    # Initialize Bert mask with masked out off diagonals.
    msa_masks = [np.ones(x['msa'].shape, dtype=np.float32)
                 for x in np_chains_list]

    np_example['bert_mask'] = block_diag(
        *msa_masks, pad_value=0)
  else:
    np_example['cluster_bias_mask'] = np.zeros(np_example['msa'].shape[0])
    np_example['cluster_bias_mask'][0] = 1

    # Initialize Bert mask with masked out off diagonals.
    msa_masks = [np.ones(x['msa'].shape, dtype=np.float32) for
                 x in np_chains_list]
    msa_masks_all_seq = [np.ones(x['msa_all_seq'].shape, dtype=np.float32) for
                         x in np_chains_list]

    msa_mask_block_diag = block_diag(
        *msa_masks, pad_value=0)
    msa_mask_all_seq = np.concatenate(msa_masks_all_seq, axis=1)
    np_example['bert_mask'] = np.concatenate(
        [msa_mask_all_seq, msa_mask_block_diag], axis=0)
  return np_example


def _pad_templates(chains: Sequence[pipeline.FeatureDict],
                   max_templates: int) -> Sequence[pipeline.FeatureDict]:
  """For each chain pad the number of templates to a fixed size.

  Args:
    chains: A list of protein chains.
    max_templates: Each chain will be padded to have this many templates.

  Returns:
    The list of chains, updated to have template features padded to
    max_templates.
  """
  for chain in chains:
    for k, v in chain.items():
      if k in TEMPLATE_FEATURES:
        padding = np.zeros_like(v.shape)
        padding[0] = max_templates - v.shape[0]
        padding = [(0, p) for p in padding]
        chain[k] = np.pad(v, padding, mode='constant')
  return chains


def _merge_features_from_multiple_chains(
    chains: Sequence[pipeline.FeatureDict],
    pair_msa_sequences: bool) -> pipeline.FeatureDict:
  """Merge features from multiple chains.

  Args:
    chains: A list of feature dictionaries that we want to merge.
    pair_msa_sequences: Whether to concatenate MSA features along the
      num_res dimension (if True), or to block diagonalize them (if False).

  Returns:
    A feature dictionary for the merged example.
  """
  merged_example = {}
  for feature_name in chains[0]:
    feats = [x[feature_name] for x in chains]
    feature_name_split = feature_name.split('_all_seq')[0]
    if feature_name_split in MSA_FEATURES:
      if pair_msa_sequences or '_all_seq' in feature_name:
        merged_example[feature_name] = np.concatenate(feats, axis=1)
      else:
        merged_example[feature_name] = block_diag(
            *feats, pad_value=MSA_PAD_VALUES[feature_name])
    elif feature_name_split in SEQ_FEATURES:
      merged_example[feature_name] = np.concatenate(feats, axis=0)
    elif feature_name_split in TEMPLATE_FEATURES:
      merged_example[feature_name] = np.concatenate(feats, axis=1)
    elif feature_name_split in CHAIN_FEATURES:
      merged_example[feature_name] = np.sum(x for x in feats).astype(np.int32)
    else:
      merged_example[feature_name] = feats[0]
  return merged_example


def _merge_homomers_dense_msa(
    chains: Iterable[pipeline.FeatureDict]) -> Sequence[pipeline.FeatureDict]:
  """Merge all identical chains, making the resulting MSA dense.

  Args:
    chains: An iterable of features for each chain.

  Returns:
    A list of feature dictionaries.  All features with the same entity_id
    will be merged - MSA features will be concatenated along the num_res
    dimension - making them dense.
  """
  entity_chains = collections.defaultdict(list)
  for chain in chains:
    entity_id = chain['entity_id'][0]
    entity_chains[entity_id].append(chain)

  grouped_chains = []
  for entity_id in sorted(entity_chains):
    chains = entity_chains[entity_id]
    grouped_chains.append(chains)
  chains = [
      _merge_features_from_multiple_chains(chains, pair_msa_sequences=True)
      for chains in grouped_chains]
  return chains


def _concatenate_paired_and_unpaired_features(
    example: pipeline.FeatureDict) -> pipeline.FeatureDict:
  """Merges paired and block-diagonalised features."""
  features = MSA_FEATURES
  for feature_name in features:
    if feature_name in example:
      feat = example[feature_name]
      feat_all_seq = example[feature_name + '_all_seq']
      merged_feat = np.concatenate([feat_all_seq, feat], axis=0)
      example[feature_name] = merged_feat
  example['num_alignments'] = np.array(example['msa'].shape[0],
                                       dtype=np.int32)
  return example


def merge_chain_features(np_chains_list: List[pipeline.FeatureDict],
                         pair_msa_sequences: bool,
                         max_templates: int) -> pipeline.FeatureDict:
  """Merges features for multiple chains to single FeatureDict.

  Args:
    np_chains_list: List of FeatureDicts for each chain.
    pair_msa_sequences: Whether to merge paired MSAs.
    max_templates: The maximum number of templates to include.

  Returns:
    Single FeatureDict for entire complex.
  """
  np_chains_list = _pad_templates(
      np_chains_list, max_templates=max_templates)
  np_chains_list = _merge_homomers_dense_msa(np_chains_list)
  # Unpaired MSA features will be always block-diagonalised; paired MSA
  # features will be concatenated.
  np_example = _merge_features_from_multiple_chains(
      np_chains_list, pair_msa_sequences=False)
  if pair_msa_sequences:
    np_example = _concatenate_paired_and_unpaired_features(np_example)
  np_example = _correct_post_merged_feats(
      np_example=np_example,
      np_chains_list=np_chains_list,
      pair_msa_sequences=pair_msa_sequences)

  return np_example


def deduplicate_unpaired_sequences(
    np_chains: List[pipeline.FeatureDict]) -> List[pipeline.FeatureDict]:
  """Removes unpaired sequences which duplicate a paired sequence."""

  feature_names = np_chains[0].keys()
  msa_features = MSA_FEATURES

  for chain in np_chains:
    sequence_set = set(tuple(s) for s in chain['msa_all_seq'])
    keep_rows = []
    # Go through unpaired MSA seqs and remove any rows that correspond to the
    # sequences that are already present in the paired MSA.
    for row_num, seq in enumerate(chain['msa']):
      if tuple(seq) not in sequence_set:
        keep_rows.append(row_num)
    for feature_name in feature_names:
      if feature_name in msa_features:
        if keep_rows:
          chain[feature_name] = chain[feature_name][keep_rows]
        else:
          new_shape = list(chain[feature_name].shape)
          new_shape[0] = 0
          chain[feature_name] = np.zeros(new_shape,
                                         dtype=chain[feature_name].dtype)
    chain['num_alignments'] = np.array(chain['msa'].shape[0], dtype=np.int32)
  return np_chains
