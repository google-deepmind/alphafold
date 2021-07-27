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

"""Tests for utils."""

import os

from absl.testing import absltest
from alphafold.common import protein
from alphafold.relax import utils
import numpy as np
# Internal import (7716).


class UtilsTest(absltest.TestCase):

  def test_overwrite_b_factors(self):
    testdir = os.path.join(
        absltest.get_default_test_srcdir(),
        'alphafold/relax/testdata/'
        'multiple_disulfides_target.pdb')
    with open(testdir) as f:
      test_pdb = f.read()
    n_residues = 191
    bfactors = np.stack([np.arange(0, n_residues)] * 37, axis=-1)

    output_pdb = utils.overwrite_b_factors(test_pdb, bfactors)

    # Check that the atom lines are unchanged apart from the B-factors.
    atom_lines_original = [l for l in test_pdb.split('\n') if l[:4] == ('ATOM')]
    atom_lines_new = [l for l in output_pdb.split('\n') if l[:4] == ('ATOM')]
    for line_original, line_new in zip(atom_lines_original, atom_lines_new):
      self.assertEqual(line_original[:60].strip(), line_new[:60].strip())
      self.assertEqual(line_original[66:].strip(), line_new[66:].strip())

    # Check B-factors are correctly set for all atoms present.
    as_protein = protein.from_pdb_string(output_pdb)
    np.testing.assert_almost_equal(
        np.where(as_protein.atom_mask > 0, as_protein.b_factors, 0),
        np.where(as_protein.atom_mask > 0, bfactors, 0))


if __name__ == '__main__':
  absltest.main()
