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

"""Tests for relax."""
import os

from absl.testing import absltest
from alphafold.common import protein
from alphafold.relax import relax
import numpy as np
# Internal import (7716).


class RunAmberRelaxTest(absltest.TestCase):

  def setUp(self):
    super().setUp()
    self.test_dir = os.path.join(
        absltest.get_default_test_srcdir(),
        'alphafold/relax/testdata/')
    self.test_config = {
        'max_iterations': 1,
        'tolerance': 2.39,
        'stiffness': 10.0,
        'exclude_residues': [],
        'max_outer_iterations': 1,
        'use_gpu': False}

  def test_process(self):
    amber_relax = relax.AmberRelaxation(**self.test_config)

    with open(os.path.join(self.test_dir, 'model_output.pdb')) as f:
      test_prot = protein.from_pdb_string(f.read())
    pdb_min, debug_info, num_violations = amber_relax.process(prot=test_prot)

    self.assertCountEqual(debug_info.keys(),
                          set({'initial_energy', 'final_energy',
                               'attempts', 'rmsd'}))
    self.assertLess(debug_info['final_energy'], debug_info['initial_energy'])
    self.assertGreater(debug_info['rmsd'], 0)

    prot_min = protein.from_pdb_string(pdb_min)
    # Most protein properties should be unchanged.
    np.testing.assert_almost_equal(test_prot.aatype, prot_min.aatype)
    np.testing.assert_almost_equal(test_prot.residue_index,
                                   prot_min.residue_index)
    # Atom mask and bfactors identical except for terminal OXT of last residue.
    np.testing.assert_almost_equal(test_prot.atom_mask[:-1, :],
                                   prot_min.atom_mask[:-1, :])
    np.testing.assert_almost_equal(test_prot.b_factors[:-1, :],
                                   prot_min.b_factors[:-1, :])
    np.testing.assert_almost_equal(test_prot.atom_mask[:, :-1],
                                   prot_min.atom_mask[:, :-1])
    np.testing.assert_almost_equal(test_prot.b_factors[:, :-1],
                                   prot_min.b_factors[:, :-1])
    # There are no residues with violations.
    np.testing.assert_equal(num_violations, np.zeros_like(num_violations))

  def test_unresolved_violations(self):
    amber_relax = relax.AmberRelaxation(**self.test_config)
    with open(os.path.join(self.test_dir,
                                 'with_violations_casp14.pdb')) as f:
      test_prot = protein.from_pdb_string(f.read())
    _, _, num_violations = amber_relax.process(prot=test_prot)
    exp_num_violations = np.array(
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1,
         1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0,
         0, 0, 0, 0])
    # Check no violations were added. Can't check exactly due to stochasticity.
    self.assertTrue(np.all(np.array(num_violations) <= exp_num_violations))


if __name__ == '__main__':
  absltest.main()
