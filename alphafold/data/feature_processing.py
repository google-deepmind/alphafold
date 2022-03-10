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

"""Feature processing logic for multimer data pipeline."""

from typing import Iterable, MutableMapping, List

from alphafold.common import residue_constants
from alphafold.data import msa_pairing
from alphafold.data import pipeline
import numpy as np

REQUIRED_FEATURES = frozenset({
    'aatype', 'all_atom_mask', 'all_atom_positions', 'all_chains_entity_ids',
    'all_crops_all_chains_mask', 'all_crops_all_chains_positions',
    'all_crops_all_chains_residue_ids', 'assembly_num_chains', 'asym_id',
    'bert_mask', 'cluster_bias_mask', 'deletion_matrix', 'deletion_mean',
    'entity_id', 'entity_mask', 'mem_peak', 'msa', 'msa_mask', 'num_alignments',
    'num_templates', 'queue_size', 'residue_index', 'resolution',
    'seq_length', 'seq_mask', 'sym_id', 'template_aatype',
    'template_all_atom_mask', 'template_all_atom_positions'
})

MAX_TEMPLATES = 4
MSA_CROP_SIZE = 2048


def _is_homomer_or_monomer(chains: Iterable[pipeline.FeatureDict]) -> bool:
  """Checks if a list of chains represents a homomer/monomer example."""
  # Note that an entity_id of 0 indicates padding.
  num_unique_chains = len(np.unique(np.concatenate(
      [np.unique(chain['entity_id'][chain['entity_id'] > 0]) for
       chain in chains])))
  return num_unique_chains == 1


def pair_and_merge(
    all_chain_features: MutableMapping[str, pipeline.FeatureDict]
    ) -> pipeline.FeatureDict:
  """Runs processing on features to augment, pair and merge.

  Args:
    all_chain_features: A MutableMap of dictionaries of features for each chain.

  Returns:
    A dictionary of features.
  """

  process_unmerged_features(all_chain_features)

  np_chains_list = list(all_chain_features.values())

  pair_msa_sequences = not _is_homomer_or_monomer(np_chains_list)

  if pair_msa_sequences:
    np_chains_list = msa_pairing.create_paired_features(
        chains=np_chains_list)
    np_chains_list = msa_pairing.deduplicate_unpaired_sequences(np_chains_list)
  np_chains_list = crop_chains(
      np_chains_list,
      msa_crop_size=MSA_CROP_SIZE,
      pair_msa_sequences=pair_msa_sequences,
      max_templates=MAX_TEMPLATES)
  np_example = msa_pairing.merge_chain_features(
      np_chains_list=np_chains_list, pair_msa_sequences=pair_msa_sequences,
      max_templates=MAX_TEMPLATES)
  np_example = process_final(np_example)
  return np_example


def crop_chains(
    chains_list: List[pipeline.FeatureDict],
    msa_crop_size: int,
    pair_msa_sequences: bool,
    max_templates: int) -> List[pipeline.FeatureDict]:
  """Crops the MSAs for a set of chains.

  Args:
    chains_list: A list of chains to be cropped.
    msa_crop_size: The total number of sequences to crop from the MSA.
    pair_msa_sequences: Whether we are operating in sequence-pairing mode.
    max_templates: The maximum templates to use per chain.

  Returns:
    The chains cropped.
  """

  # Apply the cropping.
  cropped_chains = []
  for chain in chains_list:
    cropped_chain = _crop_single_chain(
        chain,
        msa_crop_size=msa_crop_size,
        pair_msa_sequences=pair_msa_sequences,
        max_templates=max_templates)
    cropped_chains.append(cropped_chain)

  return cropped_chains


def _crop_single_chain(chain: pipeline.FeatureDict,
                       msa_crop_size: int,
                       pair_msa_sequences: bool,
                       max_templates: int) -> pipeline.FeatureDict:
  """Crops msa sequences to `msa_crop_size`."""
  msa_size = chain['num_alignments']

  if pair_msa_sequences:
    msa_size_all_seq = chain['num_alignments_all_seq']
    msa_crop_size_all_seq = np.minimum(msa_size_all_seq, msa_crop_size // 2)

    # We reduce the number of un-paired sequences, by the number of times a
    # sequence from this chain's MSA is included in the paired MSA.  This keeps
    # the MSA size for each chain roughly constant.
    msa_all_seq = chain['msa_all_seq'][:msa_crop_size_all_seq, :]
    num_non_gapped_pairs = np.sum(
        np.any(msa_all_seq != msa_pairing.MSA_GAP_IDX, axis=1))
    num_non_gapped_pairs = np.minimum(num_non_gapped_pairs,
                                      msa_crop_size_all_seq)

    # Restrict the unpaired crop size so that paired+unpaired sequences do not
    # exceed msa_seqs_per_chain for each chain.
    max_msa_crop_size = np.maximum(msa_crop_size - num_non_gapped_pairs, 0)
    msa_crop_size = np.minimum(msa_size, max_msa_crop_size)
  else:
    msa_crop_size = np.minimum(msa_size, msa_crop_size)

  include_templates = 'template_aatype' in chain and max_templates
  if include_templates:
    num_templates = chain['template_aatype'].shape[0]
    templates_crop_size = np.minimum(num_templates, max_templates)

  for k in chain:
    k_split = k.split('_all_seq')[0]
    if k_split in msa_pairing.TEMPLATE_FEATURES:
      chain[k] = chain[k][:templates_crop_size, :]
    elif k_split in msa_pairing.MSA_FEATURES:
      if '_all_seq' in k and pair_msa_sequences:
        chain[k] = chain[k][:msa_crop_size_all_seq, :]
      else:
        chain[k] = chain[k][:msa_crop_size, :]

  chain['num_alignments'] = np.asarray(msa_crop_size, dtype=np.int32)
  if include_templates:
    chain['num_templates'] = np.asarray(templates_crop_size, dtype=np.int32)
  if pair_msa_sequences:
    chain['num_alignments_all_seq'] = np.asarray(
        msa_crop_size_all_seq, dtype=np.int32)
  return chain


def process_final(np_example: pipeline.FeatureDict) -> pipeline.FeatureDict:
  """Final processing steps in data pipeline, after merging and pairing."""
  np_example = _correct_msa_restypes(np_example)
  np_example = _make_seq_mask(np_example)
  np_example = _make_msa_mask(np_example)
  np_example = _filter_features(np_example)
  return np_example


def _correct_msa_restypes(np_example):
  """Correct MSA restype to have the same order as residue_constants."""
  new_order_list = residue_constants.MAP_HHBLITS_AATYPE_TO_OUR_AATYPE
  np_example['msa'] = np.take(new_order_list, np_example['msa'], axis=0)
  np_example['msa'] = np_example['msa'].astype(np.int32)
  return np_example


def _make_seq_mask(np_example):
  np_example['seq_mask'] = (np_example['entity_id'] > 0).astype(np.float32)
  return np_example


def _make_msa_mask(np_example):
  """Mask features are all ones, but will later be zero-padded."""

  np_example['msa_mask'] = np.ones_like(np_example['msa'], dtype=np.float32)

  seq_mask = (np_example['entity_id'] > 0).astype(np.float32)
  np_example['msa_mask'] *= seq_mask[None]

  return np_example


def _filter_features(np_example: pipeline.FeatureDict) -> pipeline.FeatureDict:
  """Filters features of example to only those requested."""
  return {k: v for (k, v) in np_example.items() if k in REQUIRED_FEATURES}


def process_unmerged_features(
    all_chain_features: MutableMapping[str, pipeline.FeatureDict]):
  """Postprocessing stage for per-chain features before merging."""
  num_chains = len(all_chain_features)
  for chain_features in all_chain_features.values():
    # Convert deletion matrices to float.
    chain_features['deletion_matrix'] = np.asarray(
        chain_features.pop('deletion_matrix_int'), dtype=np.float32)
    if 'deletion_matrix_int_all_seq' in chain_features:
      chain_features['deletion_matrix_all_seq'] = np.asarray(
          chain_features.pop('deletion_matrix_int_all_seq'), dtype=np.float32)

    chain_features['deletion_mean'] = np.mean(
        chain_features['deletion_matrix'], axis=0)

    # Add all_atom_mask and dummy all_atom_positions based on aatype.
    all_atom_mask = residue_constants.STANDARD_ATOM_MASK[
        chain_features['aatype']]
    chain_features['all_atom_mask'] = all_atom_mask
    chain_features['all_atom_positions'] = np.zeros(
        list(all_atom_mask.shape) + [3])

    # Add assembly_num_chains.
    chain_features['assembly_num_chains'] = np.asarray(num_chains)

  # Add entity_mask.
  for chain_features in all_chain_features.values():
    chain_features['entity_mask'] = (
        chain_features['entity_id'] != 0).astype(np.int32)
