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

"""Test that residue_constants generates correct values."""

from absl.testing import absltest
from absl.testing import parameterized
from alphafold.common import residue_constants
import numpy as np


class ResidueConstantsTest(parameterized.TestCase):

  @parameterized.parameters(
      ('ALA', 0),
      ('CYS', 1),
      ('HIS', 2),
      ('MET', 3),
      ('LYS', 4),
      ('ARG', 4),
  )
  def testChiAnglesAtoms(self, residue_name, chi_num):
    chi_angles_atoms = residue_constants.chi_angles_atoms[residue_name]
    self.assertLen(chi_angles_atoms, chi_num)
    for chi_angle_atoms in chi_angles_atoms:
      self.assertLen(chi_angle_atoms, 4)

  def testChiGroupsForAtom(self):
    for k, chi_groups in residue_constants.chi_groups_for_atom.items():
      res_name, atom_name = k
      for chi_group_i, atom_i in chi_groups:
        self.assertEqual(
            atom_name,
            residue_constants.chi_angles_atoms[res_name][chi_group_i][atom_i])

  @parameterized.parameters(
      ('ALA', 5), ('ARG', 11), ('ASN', 8), ('ASP', 8), ('CYS', 6), ('GLN', 9),
      ('GLU', 9), ('GLY', 4), ('HIS', 10), ('ILE', 8), ('LEU', 8), ('LYS', 9),
      ('MET', 8), ('PHE', 11), ('PRO', 7), ('SER', 6), ('THR', 7), ('TRP', 14),
      ('TYR', 12), ('VAL', 7)
  )
  def testResidueAtoms(self, atom_name, num_residue_atoms):
    residue_atoms = residue_constants.residue_atoms[atom_name]
    self.assertLen(residue_atoms, num_residue_atoms)

  def testStandardAtomMask(self):
    with self.subTest('Check shape'):
      self.assertEqual(residue_constants.STANDARD_ATOM_MASK.shape, (21, 37,))

    with self.subTest('Check values'):
      str_to_row = lambda s: [c == '1' for c in s]  # More clear/concise.
      np.testing.assert_array_equal(
          residue_constants.STANDARD_ATOM_MASK,
          np.array([
              # NB This was defined by c+p but looks sane.
              str_to_row('11111                                '),  # ALA
              str_to_row('111111     1           1     11 1    '),  # ARG
              str_to_row('111111         11                    '),  # ASP
              str_to_row('111111          11                   '),  # ASN
              str_to_row('11111     1                          '),  # CYS
              str_to_row('111111     1             11          '),  # GLU
              str_to_row('111111     1              11         '),  # GLN
              str_to_row('111 1                                '),  # GLY
              str_to_row('111111       11     1    1           '),  # HIS
              str_to_row('11111 11    1                        '),  # ILE
              str_to_row('111111      11                       '),  # LEU
              str_to_row('111111     1       1               1 '),  # LYS
              str_to_row('111111            11                 '),  # MET
              str_to_row('111111      11      11          1    '),  # PHE
              str_to_row('111111     1                         '),  # PRO
              str_to_row('11111   1                            '),  # SER
              str_to_row('11111  1 1                           '),  # THR
              str_to_row('111111      11       11 1   1    11  '),  # TRP
              str_to_row('111111      11      11         11    '),  # TYR
              str_to_row('11111 11                             '),  # VAL
              str_to_row('                                     '),  # UNK
          ]))

    with self.subTest('Check row totals'):
      # Check each row has the right number of atoms.
      for row, restype in enumerate(residue_constants.restypes):  # A, R, ...
        long_restype = residue_constants.restype_1to3[restype]  # ALA, ARG, ...
        atoms_names = residue_constants.residue_atoms[
            long_restype]  # ['C', 'CA', 'CB', 'N', 'O'], ...
        self.assertLen(atoms_names,
                       residue_constants.STANDARD_ATOM_MASK[row, :].sum(),
                       long_restype)

  def testAtomTypes(self):
    self.assertEqual(residue_constants.atom_type_num, 37)

    self.assertEqual(residue_constants.atom_types[0], 'N')
    self.assertEqual(residue_constants.atom_types[1], 'CA')
    self.assertEqual(residue_constants.atom_types[2], 'C')
    self.assertEqual(residue_constants.atom_types[3], 'CB')
    self.assertEqual(residue_constants.atom_types[4], 'O')

    self.assertEqual(residue_constants.atom_order['N'], 0)
    self.assertEqual(residue_constants.atom_order['CA'], 1)
    self.assertEqual(residue_constants.atom_order['C'], 2)
    self.assertEqual(residue_constants.atom_order['CB'], 3)
    self.assertEqual(residue_constants.atom_order['O'], 4)
    self.assertEqual(residue_constants.atom_type_num, 37)

  def testRestypes(self):
    three_letter_restypes = [
        residue_constants.restype_1to3[r] for r  in residue_constants.restypes]
    for restype, exp_restype in zip(
        three_letter_restypes, sorted(residue_constants.restype_1to3.values())):
      self.assertEqual(restype, exp_restype)
    self.assertEqual(residue_constants.restype_num, 20)

  def testSequenceToOneHotHHBlits(self):
    one_hot = residue_constants.sequence_to_onehot(
        'ABCDEFGHIJKLMNOPQRSTUVWXYZ-', residue_constants.HHBLITS_AA_TO_ID)
    exp_one_hot = np.array(
        [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
         [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
         [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]])
    np.testing.assert_array_equal(one_hot, exp_one_hot)

  def testSequenceToOneHotStandard(self):
    one_hot = residue_constants.sequence_to_onehot(
        'ARNDCQEGHILKMFPSTWYV', residue_constants.restype_order)
    np.testing.assert_array_equal(one_hot, np.eye(20))

  def testSequenceToOneHotUnknownMapping(self):
    seq = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    expected_out = np.zeros([26, 21])
    for row, position in enumerate(
        [0, 20, 4, 3, 6, 13, 7, 8, 9, 20, 11, 10, 12, 2, 20, 14, 5, 1, 15, 16,
         20, 19, 17, 20, 18, 20]):
      expected_out[row, position] = 1
    aa_types = residue_constants.sequence_to_onehot(
        sequence=seq,
        mapping=residue_constants.restype_order_with_x,
        map_unknown_to_x=True)
    self.assertTrue((aa_types == expected_out).all())

  @parameterized.named_parameters(
      ('lowercase', 'aaa'),  # Insertions in A3M.
      ('gaps', '---'),  # Gaps in A3M.
      ('dots', '...'),  # Gaps in A3M.
      ('metadata', '>TEST'),  # FASTA metadata line.
  )
  def testSequenceToOneHotUnknownMappingError(self, seq):
    with self.assertRaises(ValueError):
      residue_constants.sequence_to_onehot(
          sequence=seq,
          mapping=residue_constants.restype_order_with_x,
          map_unknown_to_x=True)


if __name__ == '__main__':
  absltest.main()
