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

"""Functions for building the input features for the AlphaFold model."""

import os
from typing import Mapping, Sequence

import numpy as np

# Internal import (7716).

from alphafold.common import residue_constants
from alphafold.data import parsers
from alphafold.data import templates
from alphafold.data.tools import hhblits
from alphafold.data.tools import hhsearch
from alphafold.data.tools import jackhmmer

FeatureDict = Mapping[str, np.ndarray]


def make_sequence_features(
    sequence: str, description: str, num_res: int) -> FeatureDict:
  """Constructs a feature dict of sequence features."""
  features = {}
  features['aatype'] = residue_constants.sequence_to_onehot(
      sequence=sequence,
      mapping=residue_constants.restype_order_with_x,
      map_unknown_to_x=True)
  features['between_segment_residues'] = np.zeros((num_res,), dtype=np.int32)
  features['domain_name'] = np.array([description.encode('utf-8')],
                                     dtype=np.object_)
  features['residue_index'] = np.array(range(num_res), dtype=np.int32)
  features['seq_length'] = np.array([num_res] * num_res, dtype=np.int32)
  features['sequence'] = np.array([sequence.encode('utf-8')], dtype=np.object_)
  return features


def make_msa_features(
    msas: Sequence[Sequence[str]],
    deletion_matrices: Sequence[parsers.DeletionMatrix]) -> FeatureDict:
  """Constructs a feature dict of MSA features."""
  if not msas:
    raise ValueError('At least one MSA must be provided.')

  int_msa = []
  deletion_matrix = []
  seen_sequences = set()
  for msa_index, msa in enumerate(msas):
    if not msa:
      raise ValueError(f'MSA {msa_index} must contain at least one sequence.')
    for sequence_index, sequence in enumerate(msa):
      if sequence in seen_sequences:
        continue
      seen_sequences.add(sequence)
      int_msa.append(
          [residue_constants.HHBLITS_AA_TO_ID[res] for res in sequence])
      deletion_matrix.append(deletion_matrices[msa_index][sequence_index])

  num_res = len(msas[0][0])
  num_alignments = len(int_msa)
  features = {}
  features['deletion_matrix_int'] = np.array(deletion_matrix, dtype=np.int32)
  features['msa'] = np.array(int_msa, dtype=np.int32)
  features['num_alignments'] = np.array(
      [num_alignments] * num_res, dtype=np.int32)
  return features


class DataPipeline:
  """Runs the alignment tools and assembles the input features."""

  def __init__(self,
               jackhmmer_binary_path: str,
               hhblits_binary_path: str,
               hhsearch_binary_path: str,
               uniref90_database_path: str,
               mgnify_database_path: str,
               bfd_database_path: str,
               uniclust30_database_path: str,
               pdb70_database_path: str,
               template_featurizer: templates.TemplateHitFeaturizer,
               mgnify_max_hits: int = 501,
               uniref_max_hits: int = 10000):
    """Constructs a feature dict for a given FASTA file."""
    self.jackhmmer_uniref90_runner = jackhmmer.Jackhmmer(
        binary_path=jackhmmer_binary_path,
        database_path=uniref90_database_path)
    self.hhblits_bfd_uniclust_runner = hhblits.HHBlits(
        binary_path=hhblits_binary_path,
        databases=[bfd_database_path, uniclust30_database_path])
    self.jackhmmer_mgnify_runner = jackhmmer.Jackhmmer(
        binary_path=jackhmmer_binary_path,
        database_path=mgnify_database_path)
    self.hhsearch_pdb70_runner = hhsearch.HHSearch(
        binary_path=hhsearch_binary_path,
        databases=[pdb70_database_path])
    self.template_featurizer = template_featurizer
    self.mgnify_max_hits = mgnify_max_hits
    self.uniref_max_hits = uniref_max_hits

  def process(self, input_fasta_path: str, msa_output_dir: str) -> FeatureDict:
    """Runs alignment tools on the input sequence and creates features."""
    with open(input_fasta_path) as f:
      input_fasta_str = f.read()
    input_seqs, input_descs = parsers.parse_fasta(input_fasta_str)
    if len(input_seqs) != 1:
      raise ValueError(
          f'More than one input sequence found in {input_fasta_path}.')
    input_sequence = input_seqs[0]
    input_description = input_descs[0]
    num_res = len(input_sequence)

    jackhmmer_uniref90_result = self.jackhmmer_uniref90_runner.query(
        input_fasta_path)
    jackhmmer_mgnify_result = self.jackhmmer_mgnify_runner.query(
        input_fasta_path)

    uniref90_msa_as_a3m = parsers.convert_stockholm_to_a3m(
        jackhmmer_uniref90_result['sto'], max_sequences=self.uniref_max_hits)
    hhsearch_result = self.hhsearch_pdb70_runner.query(uniref90_msa_as_a3m)

    uniref90_out_path = os.path.join(msa_output_dir, 'uniref90_hits.sto')
    with open(uniref90_out_path, 'w') as f:
      f.write(jackhmmer_uniref90_result['sto'])

    mgnify_out_path = os.path.join(msa_output_dir, 'mgnify_hits.sto')
    with open(mgnify_out_path, 'w') as f:
      f.write(jackhmmer_mgnify_result['sto'])

    uniref90_msa, uniref90_deletion_matrix = parsers.parse_stockholm(
        jackhmmer_uniref90_result['sto'])
    mgnify_msa, mgnify_deletion_matrix = parsers.parse_stockholm(
        jackhmmer_mgnify_result['sto'])
    hhsearch_hits = parsers.parse_hhr(hhsearch_result)
    mgnify_msa = mgnify_msa[:self.mgnify_max_hits]
    mgnify_deletion_matrix = mgnify_deletion_matrix[:self.mgnify_max_hits]

    hhblits_bfd_uniclust_result = self.hhblits_bfd_uniclust_runner.query(
        input_fasta_path)

    bfd_out_path = os.path.join(msa_output_dir, 'bfd_uniclust_hits.a3m')
    with open(bfd_out_path, 'w') as f:
      f.write(hhblits_bfd_uniclust_result['a3m'])

    bfd_msa, bfd_deletion_matrix = parsers.parse_a3m(
        hhblits_bfd_uniclust_result['a3m'])

    templates_result = self.template_featurizer.get_templates(
        query_sequence=input_sequence,
        query_pdb_code=None,
        query_release_date=None,
        hhr_hits=hhsearch_hits)

    sequence_features = make_sequence_features(
        sequence=input_sequence,
        description=input_description,
        num_res=num_res)

    msa_features = make_msa_features(
        msas=(uniref90_msa, bfd_msa, mgnify_msa),
        deletion_matrices=(uniref90_deletion_matrix,
                           bfd_deletion_matrix,
                           mgnify_deletion_matrix))

    return {**sequence_features, **msa_features, **templates_result.features}
