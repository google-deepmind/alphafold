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

"""Protein data type."""

import collections
import dataclasses
import functools
import io
from typing import Any, Dict, List, Mapping, Optional, Tuple
from alphafold.common import mmcif_metadata
from alphafold.common import residue_constants
from Bio.PDB import MMCIFParser
from Bio.PDB import PDBParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.Structure import Structure
import numpy as np

FeatureDict = Mapping[str, np.ndarray]
ModelOutput = Mapping[str, Any]  # Is a nested dict.

# Complete sequence of chain IDs supported by the PDB format.
PDB_CHAIN_IDS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
PDB_MAX_CHAINS = len(PDB_CHAIN_IDS)  # := 62.

# Data to fill the _chem_comp table when writing mmCIFs.
_CHEM_COMP: Mapping[str, Tuple[Tuple[str, str], ...]] = {
    'L-peptide linking': (
        ('ALA', 'ALANINE'),
        ('ARG', 'ARGININE'),
        ('ASN', 'ASPARAGINE'),
        ('ASP', 'ASPARTIC ACID'),
        ('CYS', 'CYSTEINE'),
        ('GLN', 'GLUTAMINE'),
        ('GLU', 'GLUTAMIC ACID'),
        ('HIS', 'HISTIDINE'),
        ('ILE', 'ISOLEUCINE'),
        ('LEU', 'LEUCINE'),
        ('LYS', 'LYSINE'),
        ('MET', 'METHIONINE'),
        ('PHE', 'PHENYLALANINE'),
        ('PRO', 'PROLINE'),
        ('SER', 'SERINE'),
        ('THR', 'THREONINE'),
        ('TRP', 'TRYPTOPHAN'),
        ('TYR', 'TYROSINE'),
        ('VAL', 'VALINE'),
    ),
    'peptide linking': (('GLY', 'GLYCINE'),),
}


@dataclasses.dataclass(frozen=True)
class Protein:
  """Protein structure representation."""

  # Cartesian coordinates of atoms in angstroms. The atom types correspond to
  # residue_constants.atom_types, i.e. the first three are N, CA, CB.
  atom_positions: np.ndarray  # [num_res, num_atom_type, 3]

  # Amino-acid type for each residue represented as an integer between 0 and
  # 20, where 20 is 'X'.
  aatype: np.ndarray  # [num_res]

  # Binary float mask to indicate presence of a particular atom. 1.0 if an atom
  # is present and 0.0 if not. This should be used for loss masking.
  atom_mask: np.ndarray  # [num_res, num_atom_type]

  # Residue index as used in PDB. It is not necessarily continuous or 0-indexed.
  residue_index: np.ndarray  # [num_res]

  # 0-indexed number corresponding to the chain in the protein that this residue
  # belongs to.
  chain_index: np.ndarray  # [num_res]

  # B-factors, or temperature factors, of each residue (in sq. angstroms units),
  # representing the displacement of the residue from its ground truth mean
  # value.
  b_factors: np.ndarray  # [num_res, num_atom_type]

  def __post_init__(self):
    if len(np.unique(self.chain_index)) > PDB_MAX_CHAINS:
      raise ValueError(
          f'Cannot build an instance with more than {PDB_MAX_CHAINS} chains '
          'because these cannot be written to PDB format.')


def _from_bio_structure(
    structure: Structure, chain_id: Optional[str] = None
) -> Protein:
  """Takes a Biopython structure and creates a `Protein` instance.

  WARNING: All non-standard residue types will be converted into UNK. All
    non-standard atoms will be ignored.

  Args:
    structure: Structure from the Biopython library.
    chain_id: If chain_id is specified (e.g. A), then only that chain is parsed.
      Otherwise all chains are parsed.

  Returns:
    A new `Protein` created from the structure contents.

  Raises:
    ValueError: If the number of models included in the structure is not 1.
    ValueError: If insertion code is detected at a residue.
  """
  models = list(structure.get_models())
  if len(models) != 1:
    raise ValueError(
        'Only single model PDBs/mmCIFs are supported. Found'
        f' {len(models)} models.'
    )
  model = models[0]

  atom_positions = []
  aatype = []
  atom_mask = []
  residue_index = []
  chain_ids = []
  b_factors = []

  for chain in model:
    if chain_id is not None and chain.id != chain_id:
      continue
    for res in chain:
      if res.id[2] != ' ':
        raise ValueError(
            f'PDB/mmCIF contains an insertion code at chain {chain.id} and'
            f' residue index {res.id[1]}. These are not supported.'
        )
      res_shortname = residue_constants.restype_3to1.get(res.resname, 'X')
      restype_idx = residue_constants.restype_order.get(
          res_shortname, residue_constants.restype_num)
      pos = np.zeros((residue_constants.atom_type_num, 3))
      mask = np.zeros((residue_constants.atom_type_num,))
      res_b_factors = np.zeros((residue_constants.atom_type_num,))
      for atom in res:
        if atom.name not in residue_constants.atom_types:
          continue
        pos[residue_constants.atom_order[atom.name]] = atom.coord
        mask[residue_constants.atom_order[atom.name]] = 1.
        res_b_factors[residue_constants.atom_order[atom.name]] = atom.bfactor
      if np.sum(mask) < 0.5:
        # If no known atom positions are reported for the residue then skip it.
        continue
      aatype.append(restype_idx)
      atom_positions.append(pos)
      atom_mask.append(mask)
      residue_index.append(res.id[1])
      chain_ids.append(chain.id)
      b_factors.append(res_b_factors)

  # Chain IDs are usually characters so map these to ints.
  unique_chain_ids = np.unique(chain_ids)
  chain_id_mapping = {cid: n for n, cid in enumerate(unique_chain_ids)}
  chain_index = np.array([chain_id_mapping[cid] for cid in chain_ids])

  return Protein(
      atom_positions=np.array(atom_positions),
      atom_mask=np.array(atom_mask),
      aatype=np.array(aatype),
      residue_index=np.array(residue_index),
      chain_index=chain_index,
      b_factors=np.array(b_factors))


def from_pdb_string(pdb_str: str, chain_id: Optional[str] = None) -> Protein:
  """Takes a PDB string and constructs a `Protein` object.

  WARNING: All non-standard residue types will be converted into UNK. All
    non-standard atoms will be ignored.

  Args:
    pdb_str: The contents of the pdb file
    chain_id: If chain_id is specified (e.g. A), then only that chain is parsed.
      Otherwise all chains are parsed.

  Returns:
    A new `Protein` parsed from the pdb contents.
  """
  with io.StringIO(pdb_str) as pdb_fh:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(id='none', file=pdb_fh)
    return _from_bio_structure(structure, chain_id)


def from_mmcif_string(
    mmcif_str: str, chain_id: Optional[str] = None
) -> Protein:
  """Takes a mmCIF string and constructs a `Protein` object.

  WARNING: All non-standard residue types will be converted into UNK. All
    non-standard atoms will be ignored.

  Args:
    mmcif_str: The contents of the mmCIF file
    chain_id: If chain_id is specified (e.g. A), then only that chain is parsed.
      Otherwise all chains are parsed.

  Returns:
    A new `Protein` parsed from the mmCIF contents.
  """
  with io.StringIO(mmcif_str) as mmcif_fh:
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(structure_id='none', filename=mmcif_fh)
    return _from_bio_structure(structure, chain_id)


def _chain_end(atom_index, end_resname, chain_name, residue_index) -> str:
  chain_end = 'TER'
  return (f'{chain_end:<6}{atom_index:>5}      {end_resname:>3} '
          f'{chain_name:>1}{residue_index:>4}')


def to_pdb(prot: Protein) -> str:
  """Converts a `Protein` instance to a PDB string.

  Args:
    prot: The protein to convert to PDB.

  Returns:
    PDB string.
  """
  restypes = residue_constants.restypes + ['X']
  res_1to3 = lambda r: residue_constants.restype_1to3.get(restypes[r], 'UNK')
  atom_types = residue_constants.atom_types

  pdb_lines = []

  atom_mask = prot.atom_mask
  aatype = prot.aatype
  atom_positions = prot.atom_positions
  residue_index = prot.residue_index.astype(np.int32)
  chain_index = prot.chain_index.astype(np.int32)
  b_factors = prot.b_factors

  if np.any(aatype > residue_constants.restype_num):
    raise ValueError('Invalid aatypes.')

  # Construct a mapping from chain integer indices to chain ID strings.
  chain_ids = {}
  for i in np.unique(chain_index):  # np.unique gives sorted output.
    if i >= PDB_MAX_CHAINS:
      raise ValueError(
          f'The PDB format supports at most {PDB_MAX_CHAINS} chains.')
    chain_ids[i] = PDB_CHAIN_IDS[i]

  pdb_lines.append('MODEL     1')
  atom_index = 1
  last_chain_index = chain_index[0]
  # Add all atom sites.
  for i in range(aatype.shape[0]):
    # Close the previous chain if in a multichain PDB.
    if last_chain_index != chain_index[i]:
      pdb_lines.append(_chain_end(
          atom_index, res_1to3(aatype[i - 1]), chain_ids[chain_index[i - 1]],
          residue_index[i - 1]))
      last_chain_index = chain_index[i]
      atom_index += 1  # Atom index increases at the TER symbol.

    res_name_3 = res_1to3(aatype[i])
    for atom_name, pos, mask, b_factor in zip(
        atom_types, atom_positions[i], atom_mask[i], b_factors[i]):
      if mask < 0.5:
        continue

      record_type = 'ATOM'
      name = atom_name if len(atom_name) == 4 else f' {atom_name}'
      alt_loc = ''
      insertion_code = ''
      occupancy = 1.00
      element = atom_name[0]  # Protein supports only C, N, O, S, this works.
      charge = ''
      # PDB is a columnar format, every space matters here!
      atom_line = (f'{record_type:<6}{atom_index:>5} {name:<4}{alt_loc:>1}'
                   f'{res_name_3:>3} {chain_ids[chain_index[i]]:>1}'
                   f'{residue_index[i]:>4}{insertion_code:>1}   '
                   f'{pos[0]:>8.3f}{pos[1]:>8.3f}{pos[2]:>8.3f}'
                   f'{occupancy:>6.2f}{b_factor:>6.2f}          '
                   f'{element:>2}{charge:>2}')
      pdb_lines.append(atom_line)
      atom_index += 1

  # Close the final chain.
  pdb_lines.append(_chain_end(atom_index, res_1to3(aatype[-1]),
                              chain_ids[chain_index[-1]], residue_index[-1]))
  pdb_lines.append('ENDMDL')
  pdb_lines.append('END')

  # Pad all lines to 80 characters.
  pdb_lines = [line.ljust(80) for line in pdb_lines]
  return '\n'.join(pdb_lines) + '\n'  # Add terminating newline.


def ideal_atom_mask(prot: Protein) -> np.ndarray:
  """Computes an ideal atom mask.

  `Protein.atom_mask` typically is defined according to the atoms that are
  reported in the PDB. This function computes a mask according to heavy atoms
  that should be present in the given sequence of amino acids.

  Args:
    prot: `Protein` whose fields are `numpy.ndarray` objects.

  Returns:
    An ideal atom mask.
  """
  return residue_constants.STANDARD_ATOM_MASK[prot.aatype]


def from_prediction(
    features: FeatureDict,
    result: ModelOutput,
    b_factors: Optional[np.ndarray] = None,
    remove_leading_feature_dimension: bool = True) -> Protein:
  """Assembles a protein from a prediction.

  Args:
    features: Dictionary holding model inputs.
    result: Dictionary holding model outputs.
    b_factors: (Optional) B-factors to use for the protein.
    remove_leading_feature_dimension: Whether to remove the leading dimension
      of the `features` values.

  Returns:
    A protein instance.
  """
  fold_output = result['structure_module']

  def _maybe_remove_leading_dim(arr: np.ndarray) -> np.ndarray:
    return arr[0] if remove_leading_feature_dimension else arr

  if 'asym_id' in features:
    chain_index = _maybe_remove_leading_dim(features['asym_id'])
  else:
    chain_index = np.zeros_like(_maybe_remove_leading_dim(features['aatype']))

  if b_factors is None:
    b_factors = np.zeros_like(fold_output['final_atom_mask'])

  return Protein(
      aatype=_maybe_remove_leading_dim(features['aatype']),
      atom_positions=fold_output['final_atom_positions'],
      atom_mask=fold_output['final_atom_mask'],
      residue_index=_maybe_remove_leading_dim(features['residue_index']) + 1,
      chain_index=chain_index,
      b_factors=b_factors)


def to_mmcif(
    prot: Protein,
    file_id: str,
    model_type: str,
) -> str:
  """Converts a `Protein` instance to an mmCIF string.

  WARNING 1: The _entity_poly_seq is filled with unknown (UNK) residues for any
    missing residue indices in the range from min(1, min(residue_index)) to
    max(residue_index). E.g. for a protein object with positions for residues
    2 (MET), 3 (LYS), 6 (GLY), this method would set the _entity_poly_seq to:
    1 UNK
    2 MET
    3 LYS
    4 UNK
    5 UNK
    6 GLY
    This is done to preserve the residue numbering.

  WARNING 2: Converting ground truth mmCIF file to Protein and then back to
    mmCIF using this method will convert all non-standard residue types to UNK.
    If you need this behaviour, you need to store more mmCIF metadata in the
    Protein object (e.g. all fields except for the _atom_site loop).

  WARNING 3: Converting ground truth mmCIF file to Protein and then back to
    mmCIF using this method will not retain the original chain indices.

  WARNING 4: In case of multiple identical chains, they are assigned different
    `_atom_site.label_entity_id` values.

  Args:
    prot: A protein to convert to mmCIF string.
    file_id: The file ID (usually the PDB ID) to be used in the mmCIF.
    model_type: 'Multimer' or 'Monomer'.

  Returns:
    A valid mmCIF string.

  Raises:
    ValueError: If aminoacid types array contains entries with too many protein
    types.
  """
  atom_mask = prot.atom_mask
  aatype = prot.aatype
  atom_positions = prot.atom_positions
  residue_index = prot.residue_index.astype(np.int32)
  chain_index = prot.chain_index.astype(np.int32)
  b_factors = prot.b_factors

  # Construct a mapping from chain integer indices to chain ID strings.
  chain_ids = {}
  # We count unknown residues as protein residues.
  for entity_id in np.unique(chain_index):  # np.unique gives sorted output.
    chain_ids[entity_id] = _int_id_to_str_id(entity_id + 1)

  mmcif_dict = collections.defaultdict(list)

  mmcif_dict['data_'] = file_id.upper()
  mmcif_dict['_entry.id'] = file_id.upper()

  label_asym_id_to_entity_id = {}
  # Entity and chain information.
  for entity_id, chain_id in chain_ids.items():
    # Add all chain information to the _struct_asym table.
    label_asym_id_to_entity_id[str(chain_id)] = str(entity_id)
    mmcif_dict['_struct_asym.id'].append(chain_id)
    mmcif_dict['_struct_asym.entity_id'].append(str(entity_id))
    # Add information about the entity to the _entity_poly table.
    mmcif_dict['_entity_poly.entity_id'].append(str(entity_id))
    mmcif_dict['_entity_poly.type'].append(residue_constants.PROTEIN_CHAIN)
    mmcif_dict['_entity_poly.pdbx_strand_id'].append(chain_id)
    # Generate the _entity table.
    mmcif_dict['_entity.id'].append(str(entity_id))
    mmcif_dict['_entity.type'].append(residue_constants.POLYMER_CHAIN)

  # Add the residues to the _entity_poly_seq table.
  for entity_id, (res_ids, aas) in _get_entity_poly_seq(
      aatype, residue_index, chain_index
  ).items():
    for res_id, aa in zip(res_ids, aas):
      mmcif_dict['_entity_poly_seq.entity_id'].append(str(entity_id))
      mmcif_dict['_entity_poly_seq.num'].append(str(res_id))
      mmcif_dict['_entity_poly_seq.mon_id'].append(
          residue_constants.resnames[aa]
      )

  # Populate the chem comp table.
  for chem_type, chem_comp in _CHEM_COMP.items():
    for chem_id, chem_name in chem_comp:
      mmcif_dict['_chem_comp.id'].append(chem_id)
      mmcif_dict['_chem_comp.type'].append(chem_type)
      mmcif_dict['_chem_comp.name'].append(chem_name)

  # Add all atom sites.
  atom_index = 1
  for i in range(aatype.shape[0]):
    res_name_3 = residue_constants.resnames[aatype[i]]
    if aatype[i] <= len(residue_constants.restypes):
      atom_names = residue_constants.atom_types
    else:
      raise ValueError(
          'Amino acid types array contains entries with too many protein types.'
      )
    for atom_name, pos, mask, b_factor in zip(
        atom_names, atom_positions[i], atom_mask[i], b_factors[i]
    ):
      if mask < 0.5:
        continue
      type_symbol = residue_constants.atom_id_to_type(atom_name)

      mmcif_dict['_atom_site.group_PDB'].append('ATOM')
      mmcif_dict['_atom_site.id'].append(str(atom_index))
      mmcif_dict['_atom_site.type_symbol'].append(type_symbol)
      mmcif_dict['_atom_site.label_atom_id'].append(atom_name)
      mmcif_dict['_atom_site.label_alt_id'].append('.')
      mmcif_dict['_atom_site.label_comp_id'].append(res_name_3)
      mmcif_dict['_atom_site.label_asym_id'].append(chain_ids[chain_index[i]])
      mmcif_dict['_atom_site.label_entity_id'].append(
          label_asym_id_to_entity_id[chain_ids[chain_index[i]]]
      )
      mmcif_dict['_atom_site.label_seq_id'].append(str(residue_index[i]))
      mmcif_dict['_atom_site.pdbx_PDB_ins_code'].append('.')
      mmcif_dict['_atom_site.Cartn_x'].append(f'{pos[0]:.3f}')
      mmcif_dict['_atom_site.Cartn_y'].append(f'{pos[1]:.3f}')
      mmcif_dict['_atom_site.Cartn_z'].append(f'{pos[2]:.3f}')
      mmcif_dict['_atom_site.occupancy'].append('1.00')
      mmcif_dict['_atom_site.B_iso_or_equiv'].append(f'{b_factor:.2f}')
      mmcif_dict['_atom_site.auth_seq_id'].append(str(residue_index[i]))
      mmcif_dict['_atom_site.auth_asym_id'].append(chain_ids[chain_index[i]])
      mmcif_dict['_atom_site.pdbx_PDB_model_num'].append('1')

      atom_index += 1

  metadata_dict = mmcif_metadata.add_metadata_to_mmcif(mmcif_dict, model_type)
  mmcif_dict.update(metadata_dict)

  return _create_mmcif_string(mmcif_dict)


@functools.lru_cache(maxsize=256)
def _int_id_to_str_id(num: int) -> str:
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


def _get_entity_poly_seq(
    aatypes: np.ndarray, residue_indices: np.ndarray, chain_indices: np.ndarray
) -> Dict[int, Tuple[List[int], List[int]]]:
  """Constructs gapless residue index and aatype lists for each chain.

  Args:
    aatypes: A numpy array with aatypes.
    residue_indices: A numpy array with residue indices.
    chain_indices: A numpy array with chain indices.

  Returns:
    A dictionary mapping chain indices to a tuple with list of residue indices
    and a list of aatypes. Missing residues are filled with UNK residue type.
  """
  if (
      aatypes.shape[0] != residue_indices.shape[0]
      or aatypes.shape[0] != chain_indices.shape[0]
  ):
    raise ValueError(
        'aatypes, residue_indices, chain_indices must have the same length.'
    )

  # Group the present residues by chain index.
  present = collections.defaultdict(list)
  for chain_id, res_id, aa in zip(chain_indices, residue_indices, aatypes):
    present[chain_id].append((res_id, aa))

  # Add any missing residues (from 1 to the first residue and for any gaps).
  entity_poly_seq = {}
  for chain_id, present_residues in present.items():
    present_residue_indices = set([x[0] for x in present_residues])
    min_res_id = min(present_residue_indices)  # Could be negative.
    max_res_id = max(present_residue_indices)

    new_residue_indices = []
    new_aatypes = []
    present_index = 0
    for i in range(min(1, min_res_id), max_res_id + 1):
      new_residue_indices.append(i)
      if i in present_residue_indices:
        new_aatypes.append(present_residues[present_index][1])
        present_index += 1
      else:
        new_aatypes.append(20)  # Unknown amino acid type.
    entity_poly_seq[chain_id] = (new_residue_indices, new_aatypes)
  return entity_poly_seq


def _create_mmcif_string(mmcif_dict: Dict[str, Any]) -> str:
  """Converts mmCIF dictionary into mmCIF string."""
  mmcifio = MMCIFIO()
  mmcifio.set_dict(mmcif_dict)

  with io.StringIO() as file_handle:
    mmcifio.save(file_handle)
    return file_handle.getvalue()
