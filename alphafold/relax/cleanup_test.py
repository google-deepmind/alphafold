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

"""Tests for relax.cleanup."""
import io

from absl.testing import absltest
from alphafold.relax import cleanup
from openmm.app.internal import pdbstructure


def _pdb_to_structure(pdb_str):
  handle = io.StringIO(pdb_str)
  return pdbstructure.PdbStructure(handle)


def _lines_to_structure(pdb_lines):
  return _pdb_to_structure('\n'.join(pdb_lines))


class CleanupTest(absltest.TestCase):

  def test_missing_residues(self):
    pdb_lines = ['SEQRES   1 C    3  CYS GLY LEU',
                 'ATOM      1  N   CYS C   1     -12.262  20.115  60.959  1.00 '
                 '19.08           N',
                 'ATOM      2  CA  CYS C   1     -11.065  20.934  60.773  1.00 '
                 '17.23           C',
                 'ATOM      3  C   CYS C   1     -10.002  20.742  61.844  1.00 '
                 '15.38           C',
                 'ATOM      4  O   CYS C   1     -10.284  20.225  62.929  1.00 '
                 '16.04           O',
                 'ATOM      5  N   LEU C   3      -7.688  18.700  62.045  1.00 '
                 '14.75           N',
                 'ATOM      6  CA  LEU C   3      -7.256  17.320  62.234  1.00 '
                 '16.81           C',
                 'ATOM      7  C   LEU C   3      -6.380  16.864  61.070  1.00 '
                 '16.95           C',
                 'ATOM      8  O   LEU C   3      -6.551  17.332  59.947  1.00 '
                 '16.97           O']
    input_handle = io.StringIO('\n'.join(pdb_lines))
    alterations = {}
    result = cleanup.fix_pdb(input_handle, alterations)
    structure = _pdb_to_structure(result)
    residue_names = [r.get_name() for r in structure.iter_residues()]
    self.assertCountEqual(residue_names, ['CYS', 'GLY', 'LEU'])
    self.assertCountEqual(alterations['missing_residues'].values(), [['GLY']])

  def test_missing_atoms(self):
    pdb_lines = ['SEQRES   1 A    1  PRO',
                 'ATOM      1  CA  PRO A   1       1.000   1.000   1.000  1.00 '
                 ' 0.00           C']
    input_handle = io.StringIO('\n'.join(pdb_lines))
    alterations = {}
    result = cleanup.fix_pdb(input_handle, alterations)
    structure = _pdb_to_structure(result)
    atom_names = [a.get_name() for a in structure.iter_atoms()]
    self.assertCountEqual(atom_names, ['N', 'CD', 'HD2', 'HD3', 'CG', 'HG2',
                                       'HG3', 'CB', 'HB2', 'HB3', 'CA', 'HA',
                                       'C', 'O', 'H2', 'H3', 'OXT'])
    missing_atoms_by_residue = list(alterations['missing_heavy_atoms'].values())
    self.assertLen(missing_atoms_by_residue, 1)
    atoms_added = [a.name for a in missing_atoms_by_residue[0]]
    self.assertCountEqual(atoms_added, ['N', 'CD', 'CG', 'CB', 'C', 'O'])
    missing_terminals_by_residue = alterations['missing_terminals']
    self.assertLen(missing_terminals_by_residue, 1)
    has_missing_terminal = [r.name for r in missing_terminals_by_residue.keys()]
    self.assertCountEqual(has_missing_terminal, ['PRO'])
    self.assertCountEqual([t for t in missing_terminals_by_residue.values()],
                          [['OXT']])

  def test_remove_heterogens(self):
    pdb_lines = ['SEQRES   1 A    1  GLY',
                 'ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00 '
                 ' 0.00           C',
                 'ATOM      2   O  HOH A   2       0.000   0.000   0.000  1.00 '
                 ' 0.00           O']
    input_handle = io.StringIO('\n'.join(pdb_lines))
    alterations = {}
    result = cleanup.fix_pdb(input_handle, alterations)
    structure = _pdb_to_structure(result)
    self.assertCountEqual([res.get_name() for res in structure.iter_residues()],
                          ['GLY'])
    self.assertEqual(alterations['removed_heterogens'], set(['HOH']))

  def test_fix_nonstandard_residues(self):
    pdb_lines = ['SEQRES   1 A    1  DAL',
                 'ATOM      1  CA  DAL A   1       0.000   0.000   0.000  1.00 '
                 ' 0.00           C']
    input_handle = io.StringIO('\n'.join(pdb_lines))
    alterations = {}
    result = cleanup.fix_pdb(input_handle, alterations)
    structure = _pdb_to_structure(result)
    residue_names = [res.get_name() for res in structure.iter_residues()]
    self.assertCountEqual(residue_names, ['ALA'])
    self.assertLen(alterations['nonstandard_residues'], 1)
    original_res, new_name = alterations['nonstandard_residues'][0]
    self.assertEqual(original_res.id, '1')
    self.assertEqual(new_name, 'ALA')

  def test_replace_met_se(self):
    pdb_lines = ['SEQRES   1 A    1  MET',
                 'ATOM      1  SD  MET A   1       0.000   0.000   0.000  1.00 '
                 ' 0.00          Se']
    structure = _lines_to_structure(pdb_lines)
    alterations = {}
    cleanup._replace_met_se(structure, alterations)
    sd = [a for a in structure.iter_atoms() if a.get_name() == 'SD']
    self.assertLen(sd, 1)
    self.assertEqual(sd[0].element_symbol, 'S')
    self.assertCountEqual(alterations['Se_in_MET'], [sd[0].residue_number])

  def test_remove_chains_of_length_one(self):
    pdb_lines = ['SEQRES   1 A    1  GLY',
                 'ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00 '
                 ' 0.00           C']
    structure = _lines_to_structure(pdb_lines)
    alterations = {}
    cleanup._remove_chains_of_length_one(structure, alterations)
    chains = list(structure.iter_chains())
    self.assertEmpty(chains)
    self.assertCountEqual(alterations['removed_chains'].values(), [['A']])


if __name__ == '__main__':
  absltest.main()
