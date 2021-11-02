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

"""Tests for protein."""

import os

from absl.testing import absltest
from absl.testing import parameterized
from alphafold.common import protein
from alphafold.common import residue_constants
import numpy as np
# Internal import (7716).

TEST_DATA_DIR = 'alphafold/common/testdata/'


class ProteinTest(parameterized.TestCase):

  def _check_shapes(self, prot, num_res):
    """Check that the processed shapes are correct."""
    num_atoms = residue_constants.atom_type_num
    self.assertEqual((num_res, num_atoms, 3), prot.atom_positions.shape)
    self.assertEqual((num_res,), prot.aatype.shape)
    self.assertEqual((num_res, num_atoms), prot.atom_mask.shape)
    self.assertEqual((num_res,), prot.residue_index.shape)
    self.assertEqual((num_res,), prot.chain_index.shape)
    self.assertEqual((num_res, num_atoms), prot.b_factors.shape)

  @parameterized.named_parameters(
      dict(testcase_name='chain_A',
           pdb_file='2rbg.pdb', chain_id='A', num_res=282, num_chains=1),
      dict(testcase_name='chain_B',
           pdb_file='2rbg.pdb', chain_id='B', num_res=282, num_chains=1),
      dict(testcase_name='multichain',
           pdb_file='2rbg.pdb', chain_id=None, num_res=564, num_chains=2))
  def test_from_pdb_str(self, pdb_file, chain_id, num_res, num_chains):
    pdb_file = os.path.join(absltest.get_default_test_srcdir(), TEST_DATA_DIR,
                            pdb_file)
    with open(pdb_file) as f:
      pdb_string = f.read()
    prot = protein.from_pdb_string(pdb_string, chain_id)
    self._check_shapes(prot, num_res)
    self.assertGreaterEqual(prot.aatype.min(), 0)
    # Allow equal since unknown restypes have index equal to restype_num.
    self.assertLessEqual(prot.aatype.max(), residue_constants.restype_num)
    self.assertLen(np.unique(prot.chain_index), num_chains)

  def test_to_pdb(self):
    with open(
        os.path.join(absltest.get_default_test_srcdir(), TEST_DATA_DIR,
                     '2rbg.pdb')) as f:
      pdb_string = f.read()
    prot = protein.from_pdb_string(pdb_string)
    pdb_string_reconstr = protein.to_pdb(prot)

    for line in pdb_string_reconstr.splitlines():
      self.assertLen(line, 80)

    prot_reconstr = protein.from_pdb_string(pdb_string_reconstr)

    np.testing.assert_array_equal(prot_reconstr.aatype, prot.aatype)
    np.testing.assert_array_almost_equal(
        prot_reconstr.atom_positions, prot.atom_positions)
    np.testing.assert_array_almost_equal(
        prot_reconstr.atom_mask, prot.atom_mask)
    np.testing.assert_array_equal(
        prot_reconstr.residue_index, prot.residue_index)
    np.testing.assert_array_equal(
        prot_reconstr.chain_index, prot.chain_index)
    np.testing.assert_array_almost_equal(
        prot_reconstr.b_factors, prot.b_factors)

  def test_ideal_atom_mask(self):
    with open(
        os.path.join(absltest.get_default_test_srcdir(), TEST_DATA_DIR,
                     '2rbg.pdb')) as f:
      pdb_string = f.read()
    prot = protein.from_pdb_string(pdb_string)
    ideal_mask = protein.ideal_atom_mask(prot)
    non_ideal_residues = set([102] + list(range(127, 286)))
    for i, (res, atom_mask) in enumerate(
        zip(prot.residue_index, prot.atom_mask)):
      if res in non_ideal_residues:
        self.assertFalse(np.all(atom_mask == ideal_mask[i]), msg=f'{res}')
      else:
        self.assertTrue(np.all(atom_mask == ideal_mask[i]), msg=f'{res}')

  def test_too_many_chains(self):
    num_res = protein.PDB_MAX_CHAINS + 1
    num_atom_type = residue_constants.atom_type_num
    with self.assertRaises(ValueError):
      _ = protein.Protein(
          atom_positions=np.random.random([num_res, num_atom_type, 3]),
          aatype=np.random.randint(0, 21, [num_res]),
          atom_mask=np.random.randint(0, 2, [num_res]).astype(np.float32),
          residue_index=np.arange(1, num_res+1),
          chain_index=np.arange(num_res),
          b_factors=np.random.uniform(1, 100, [num_res]))


if __name__ == '__main__':
  absltest.main()
