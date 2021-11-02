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

"""Helper methods for the AlphaFold Colab notebook."""
import enum
import json
from typing import Any, Mapping, Optional, Sequence, Tuple

from alphafold.common import residue_constants
from alphafold.data import parsers
from matplotlib import pyplot as plt
import numpy as np


@enum.unique
class ModelType(enum.Enum):
  MONOMER = 0
  MULTIMER = 1


def clean_and_validate_sequence(
    input_sequence: str, min_length: int, max_length: int) -> str:
  """Checks that the input sequence is ok and returns a clean version of it."""
  # Remove all whitespaces, tabs and end lines; upper-case.
  clean_sequence = input_sequence.translate(
      str.maketrans('', '', ' \n\t')).upper()
  aatypes = set(residue_constants.restypes)  # 20 standard aatypes.
  if not set(clean_sequence).issubset(aatypes):
    raise ValueError(
        f'Input sequence contains non-amino acid letters: '
        f'{set(clean_sequence) - aatypes}. AlphaFold only supports 20 standard '
        'amino acids as inputs.')
  if len(clean_sequence) < min_length:
    raise ValueError(
        f'Input sequence is too short: {len(clean_sequence)} amino acids, '
        f'while the minimum is {min_length}')
  if len(clean_sequence) > max_length:
    raise ValueError(
        f'Input sequence is too long: {len(clean_sequence)} amino acids, while '
        f'the maximum is {max_length}. You may be able to run it with the full '
        f'AlphaFold system depending on your resources (system memory, '
        f'GPU memory).')
  return clean_sequence


def validate_input(
    input_sequences: Sequence[str],
    min_length: int,
    max_length: int,
    max_multimer_length: int) -> Tuple[Sequence[str], ModelType]:
  """Validates and cleans input sequences and determines which model to use."""
  sequences = []

  for input_sequence in input_sequences:
    if input_sequence.strip():
      input_sequence = clean_and_validate_sequence(
          input_sequence=input_sequence,
          min_length=min_length,
          max_length=max_length)
      sequences.append(input_sequence)

  if len(sequences) == 1:
    print('Using the single-chain model.')
    return sequences, ModelType.MONOMER

  elif len(sequences) > 1:
    total_multimer_length = sum([len(seq) for seq in sequences])
    if total_multimer_length > max_multimer_length:
      raise ValueError(f'The total length of multimer sequences is too long: '
                       f'{total_multimer_length}, while the maximum is '
                       f'{max_multimer_length}. Please use the full AlphaFold '
                       f'system for long multimers.')
    elif total_multimer_length > 1536:
      print('WARNING: The accuracy of the system has not been fully validated '
            'above 1536 residues, and you may experience long running times or '
            f'run out of memory for your complex with {total_multimer_length} '
            'residues.')
    print(f'Using the multimer model with {len(sequences)} sequences.')
    return sequences, ModelType.MULTIMER

  else:
    raise ValueError('No input amino acid sequence provided, please provide at '
                     'least one sequence.')


def merge_chunked_msa(
    results: Sequence[Mapping[str, Any]],
    max_hits: Optional[int] = None
    ) -> parsers.Msa:
  """Merges chunked database hits together into hits for the full database."""
  unsorted_results = []
  for chunk_index, chunk in enumerate(results):
    msa = parsers.parse_stockholm(chunk['sto'])
    e_values_dict = parsers.parse_e_values_from_tblout(chunk['tbl'])
    # Jackhmmer lists sequences as <sequence name>/<residue from>-<residue to>.
    e_values = [e_values_dict[t.partition('/')[0]] for t in msa.descriptions]
    chunk_results = zip(
        msa.sequences, msa.deletion_matrix, msa.descriptions, e_values)
    if chunk_index != 0:
      next(chunk_results)  # Only take query (first hit) from the first chunk.
    unsorted_results.extend(chunk_results)

  sorted_by_evalue = sorted(unsorted_results, key=lambda x: x[-1])
  merged_sequences, merged_deletion_matrix, merged_descriptions, _ = zip(
      *sorted_by_evalue)
  merged_msa = parsers.Msa(sequences=merged_sequences,
                           deletion_matrix=merged_deletion_matrix,
                           descriptions=merged_descriptions)
  if max_hits is not None:
    merged_msa = merged_msa.truncate(max_seqs=max_hits)

  return merged_msa


def show_msa_info(
    single_chain_msas: Sequence[parsers.Msa],
    sequence_index: int):
  """Prints info and shows a plot of the deduplicated single chain MSA."""
  full_single_chain_msa = []
  for single_chain_msa in single_chain_msas:
    full_single_chain_msa.extend(single_chain_msa.sequences)

  # Deduplicate but preserve order (hence can't use set).
  deduped_full_single_chain_msa = list(dict.fromkeys(full_single_chain_msa))
  total_msa_size = len(deduped_full_single_chain_msa)
  print(f'\n{total_msa_size} unique sequences found in total for sequence '
        f'{sequence_index}\n')

  aa_map = {res: i for i, res in enumerate('ABCDEFGHIJKLMNOPQRSTUVWXYZ-')}
  msa_arr = np.array(
      [[aa_map[aa] for aa in seq] for seq in deduped_full_single_chain_msa])

  plt.figure(figsize=(12, 3))
  plt.title(f'Per-Residue Count of Non-Gap Amino Acids in the MSA for Sequence '
            f'{sequence_index}')
  plt.plot(np.sum(msa_arr != aa_map['-'], axis=0), color='black')
  plt.ylabel('Non-Gap Count')
  plt.yticks(range(0, total_msa_size + 1, max(1, int(total_msa_size / 3))))
  plt.show()


def empty_placeholder_template_features(
    num_templates: int, num_res: int) -> Mapping[str, np.ndarray]:
  return {
      'template_aatype': np.zeros(
          (num_templates, num_res,
           len(residue_constants.restypes_with_x_and_gap)), dtype=np.float32),
      'template_all_atom_masks': np.zeros(
          (num_templates, num_res, residue_constants.atom_type_num),
          dtype=np.float32),
      'template_all_atom_positions': np.zeros(
          (num_templates, num_res, residue_constants.atom_type_num, 3),
          dtype=np.float32),
      'template_domain_names': np.zeros([num_templates], dtype=np.object),
      'template_sequence': np.zeros([num_templates], dtype=np.object),
      'template_sum_probs': np.zeros([num_templates], dtype=np.float32),
  }


def get_pae_json(pae: np.ndarray, max_pae: float) -> str:
  """Returns the PAE in the same format as is used in the AFDB."""
  rounded_errors = np.round(pae.astype(np.float64), decimals=1)
  indices = np.indices((len(rounded_errors), len(rounded_errors))) + 1
  indices_1 = indices[0].flatten().tolist()
  indices_2 = indices[1].flatten().tolist()
  return json.dumps(
      [{'residue1': indices_1,
        'residue2': indices_2,
        'distance': rounded_errors.flatten().tolist(),
        'max_predicted_aligned_error': max_pae}],
      indent=None, separators=(',', ':'))
