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

"""Functions for building the features for the AlphaFold multimer model."""

import collections
import contextlib
import copy
import dataclasses
import json
import os
import tempfile
from typing import Mapping, MutableMapping, Sequence

from absl import logging
from alphafold.common import protein
from alphafold.common import residue_constants
from alphafold.data import feature_processing
from alphafold.data import msa_pairing
from alphafold.data import parsers
from alphafold.data import pipeline
from alphafold.data.tools import jackhmmer
import numpy as np

# Internal import (7716).


@dataclasses.dataclass(frozen=True)
class _FastaChain:
  sequence: str
  description: str


def _make_chain_id_map(*,
                       sequences: Sequence[str],
                       descriptions: Sequence[str],
                       ) -> Mapping[str, _FastaChain]:
  """Makes a mapping from PDB-format chain ID to sequence and description."""
  if len(sequences) != len(descriptions):
    raise ValueError('sequences and descriptions must have equal length. '
                     f'Got {len(sequences)} != {len(descriptions)}.')
  if len(sequences) > protein.PDB_MAX_CHAINS:
    raise ValueError('Cannot process more chains than the PDB format supports. '
                     f'Got {len(sequences)} chains.')
  chain_id_map = {}
  for chain_id, sequence, description in zip(
      protein.PDB_CHAIN_IDS, sequences, descriptions):
    chain_id_map[chain_id] = _FastaChain(
        sequence=sequence, description=description)
  return chain_id_map


@contextlib.contextmanager
def temp_fasta_file(fasta_str: str):
  with tempfile.NamedTemporaryFile('w', suffix='.fasta') as fasta_file:
    fasta_file.write(fasta_str)
    fasta_file.seek(0)
    yield fasta_file.name


def convert_monomer_features(
    monomer_features: pipeline.FeatureDict,
    chain_id: str) -> pipeline.FeatureDict:
  """Reshapes and modifies monomer features for multimer models."""
  converted = {}
  converted['auth_chain_id'] = np.asarray(chain_id, dtype=np.object_)
  unnecessary_leading_dim_feats = {
      'sequence', 'domain_name', 'num_alignments', 'seq_length'}
  for feature_name, feature in monomer_features.items():
    if feature_name in unnecessary_leading_dim_feats:
      # asarray ensures it's a np.ndarray.
      feature = np.asarray(feature[0], dtype=feature.dtype)
    elif feature_name == 'aatype':
      # The multimer model performs the one-hot operation itself.
      feature = np.argmax(feature, axis=-1).astype(np.int32)
    elif feature_name == 'template_aatype':
      feature = np.argmax(feature, axis=-1).astype(np.int32)
      new_order_list = residue_constants.MAP_HHBLITS_AATYPE_TO_OUR_AATYPE
      feature = np.take(new_order_list, feature.astype(np.int32), axis=0)
    elif feature_name == 'template_all_atom_masks':
      feature_name = 'template_all_atom_mask'
    converted[feature_name] = feature
  return converted


def int_id_to_str_id(num: int) -> str:
  """Encodes a number as a string, using reverse spreadsheet style naming.

  Args:
    num: A positive integer.

  Returns:
    A string that encodes the positive integer using reverse spreadsheet style,
    naming e.g. 1 = A, 2 = B, ..., 27 = AA, 28 = BA, 29 = CA, ... This is the
    usual way to encode chain IDs in mmCIF files.
  """
  if num <= 0:
    raise ValueError(f'Only positive integers allowed, got {num}.')

  num = num - 1  # 1-based indexing.
  output = []
  while num >= 0:
    output.append(chr(num % 26 + ord('A')))
    num = num // 26 - 1
  return ''.join(output)


def add_assembly_features(
    all_chain_features: MutableMapping[str, pipeline.FeatureDict],
    ) -> MutableMapping[str, pipeline.FeatureDict]:
  """Add features to distinguish between chains.

  Args:
    all_chain_features: A dictionary which maps chain_id to a dictionary of
      features for each chain.

  Returns:
    all_chain_features: A dictionary which maps strings of the form
      `<seq_id>_<sym_id>` to the corresponding chain features. E.g. two
      chains from a homodimer would have keys A_1 and A_2. Two chains from a
      heterodimer would have keys A_1 and B_1.
  """
  # Group the chains by sequence
  seq_to_entity_id = {}
  grouped_chains = collections.defaultdict(list)
  for chain_id, chain_features in all_chain_features.items():
    seq = str(chain_features['sequence'])
    if seq not in seq_to_entity_id:
      seq_to_entity_id[seq] = len(seq_to_entity_id) + 1
    grouped_chains[seq_to_entity_id[seq]].append(chain_features)

  new_all_chain_features = {}
  chain_id = 1
  for entity_id, group_chain_features in grouped_chains.items():
    for sym_id, chain_features in enumerate(group_chain_features, start=1):
      new_all_chain_features[
          f'{int_id_to_str_id(entity_id)}_{sym_id}'] = chain_features
      seq_length = chain_features['seq_length']
      chain_features['asym_id'] = chain_id * np.ones(seq_length)
      chain_features['sym_id'] = sym_id * np.ones(seq_length)
      chain_features['entity_id'] = entity_id * np.ones(seq_length)
      chain_id += 1

  return new_all_chain_features


def pad_msa(np_example, min_num_seq):
  np_example = dict(np_example)
  num_seq = np_example['msa'].shape[0]
  if num_seq < min_num_seq:
    for feat in ('msa', 'deletion_matrix', 'bert_mask', 'msa_mask'):
      np_example[feat] = np.pad(
          np_example[feat], ((0, min_num_seq - num_seq), (0, 0)))
    np_example['cluster_bias_mask'] = np.pad(
        np_example['cluster_bias_mask'], ((0, min_num_seq - num_seq),))
  return np_example


class DataPipeline:
  """Runs the alignment tools and assembles the input features."""

  def __init__(self,
               monomer_data_pipeline: pipeline.DataPipeline,
               jackhmmer_binary_path: str,
               uniprot_database_path: str,
               max_uniprot_hits: int = 50000,
               use_precomputed_msas: bool = False):
    """Initializes the data pipeline.

    Args:
      monomer_data_pipeline: An instance of pipeline.DataPipeline - that runs
        the data pipeline for the monomer AlphaFold system.
      jackhmmer_binary_path: Location of the jackhmmer binary.
      uniprot_database_path: Location of the unclustered uniprot sequences, that
        will be searched with jackhmmer and used for MSA pairing.
      max_uniprot_hits: The maximum number of hits to return from uniprot.
      use_precomputed_msas: Whether to use pre-existing MSAs; see run_alphafold.
    """
    self._monomer_data_pipeline = monomer_data_pipeline
    self._uniprot_msa_runner = jackhmmer.Jackhmmer(
        binary_path=jackhmmer_binary_path,
        database_path=uniprot_database_path)
    self._max_uniprot_hits = max_uniprot_hits
    self.use_precomputed_msas = use_precomputed_msas

  def _process_single_chain(
      self,
      chain_id: str,
      sequence: str,
      description: str,
      msa_output_dir: str,
      is_homomer_or_monomer: bool) -> pipeline.FeatureDict:
    """Runs the monomer pipeline on a single chain."""
    chain_fasta_str = f'>chain_{chain_id}\n{sequence}\n'
    chain_msa_output_dir = os.path.join(msa_output_dir, chain_id)
    if not os.path.exists(chain_msa_output_dir):
      os.makedirs(chain_msa_output_dir)
    with temp_fasta_file(chain_fasta_str) as chain_fasta_path:
      logging.info('Running monomer pipeline on chain %s: %s',
                   chain_id, description)
      chain_features = self._monomer_data_pipeline.process(
          input_fasta_path=chain_fasta_path,
          msa_output_dir=chain_msa_output_dir)

      # We only construct the pairing features if there are 2 or more unique
      # sequences.
      if not is_homomer_or_monomer:
        all_seq_msa_features = self._all_seq_msa_features(chain_fasta_path,
                                                          chain_msa_output_dir)
        chain_features.update(all_seq_msa_features)
    return chain_features

  def _all_seq_msa_features(self, input_fasta_path, msa_output_dir):
    """Get MSA features for unclustered uniprot, for pairing."""
    out_path = os.path.join(msa_output_dir, 'uniprot_hits.sto')
    result = pipeline.run_msa_tool(
        self._uniprot_msa_runner, input_fasta_path, out_path, 'sto',
        self.use_precomputed_msas)
    msa = parsers.parse_stockholm(result['sto'])
    msa = msa.truncate(max_seqs=self._max_uniprot_hits)
    all_seq_features = pipeline.make_msa_features([msa])
    valid_feats = msa_pairing.MSA_FEATURES + (
        'msa_species_identifiers',
    )
    feats = {f'{k}_all_seq': v for k, v in all_seq_features.items()
             if k in valid_feats}
    return feats

  def process(self,
              input_fasta_path: str,
              msa_output_dir: str) -> pipeline.FeatureDict:
    """Runs alignment tools on the input sequences and creates features."""
    with open(input_fasta_path) as f:
      input_fasta_str = f.read()
    input_seqs, input_descs = parsers.parse_fasta(input_fasta_str)

    chain_id_map = _make_chain_id_map(sequences=input_seqs,
                                      descriptions=input_descs)
    chain_id_map_path = os.path.join(msa_output_dir, 'chain_id_map.json')
    with open(chain_id_map_path, 'w') as f:
      chain_id_map_dict = {chain_id: dataclasses.asdict(fasta_chain)
                           for chain_id, fasta_chain in chain_id_map.items()}
      json.dump(chain_id_map_dict, f, indent=4, sort_keys=True)

    all_chain_features = {}
    sequence_features = {}
    is_homomer_or_monomer = len(set(input_seqs)) == 1
    for chain_id, fasta_chain in chain_id_map.items():
      if fasta_chain.sequence in sequence_features:
        all_chain_features[chain_id] = copy.deepcopy(
            sequence_features[fasta_chain.sequence])
        continue
      chain_features = self._process_single_chain(
          chain_id=chain_id,
          sequence=fasta_chain.sequence,
          description=fasta_chain.description,
          msa_output_dir=msa_output_dir,
          is_homomer_or_monomer=is_homomer_or_monomer)

      chain_features = convert_monomer_features(chain_features,
                                                chain_id=chain_id)
      all_chain_features[chain_id] = chain_features
      sequence_features[fasta_chain.sequence] = chain_features

    all_chain_features = add_assembly_features(all_chain_features)

    np_example = feature_processing.pair_and_merge(
        all_chain_features=all_chain_features)

    # Pad MSA to avoid zero-sized extra_msa.
    np_example = pad_msa(np_example, 512)

    return np_example
