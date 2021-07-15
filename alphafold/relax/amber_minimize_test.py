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

"""Tests for amber_minimize."""
import os

from absl.testing import absltest
import numpy as np

from alphafold.common import protein
from alphafold.relax import amber_minimize
# Internal import (7716).


def _load_test_protein(data_path):
  pdb_path = os.path.join(absltest.get_default_test_srcdir(), data_path)
  with open(pdb_path, 'r') as f:
    return protein.from_pdb_string(f.read())


class AmberMinimizeTest(absltest.TestCase):

  def test_multiple_disulfides_target(self):
    prot = _load_test_protein(
        'alphafold/relax/testdata/multiple_disulfides_target.pdb'
        )
    ret = amber_minimize.run_pipeline(prot, max_iterations=10, max_attempts=1,
                                      stiffness=10.)
    self.assertIn('opt_time', ret)
    self.assertIn('min_attempts', ret)

  def test_raises_invalid_protein_assertion(self):
    prot = _load_test_protein(
        'alphafold/relax/testdata/multiple_disulfides_target.pdb'
        )
    prot.atom_mask[4, :] = 0
    with self.assertRaisesRegex(
        ValueError,
        'Amber minimization can only be performed on proteins with well-defined'
        ' residues. This protein contains at least one residue with no atoms.'):
      amber_minimize.run_pipeline(prot, max_iterations=10,
                                  stiffness=1.,
                                  max_attempts=1)

  def test_iterative_relax(self):
    # This test can occasionally fail because of nondeterminism in OpenMM.
    prot = _load_test_protein(
        'alphafold/relax/testdata/with_violations.pdb'
        )
    violations = amber_minimize.get_violation_metrics(prot)
    self.assertGreater(violations['num_residue_violations'], 0)
    out = amber_minimize.run_pipeline(
        prot=prot, max_outer_iterations=10, stiffness=10.)
    self.assertLess(out['efinal'], out['einit'])
    self.assertEqual(0, out['num_residue_violations'])

  def test_find_violations(self):
    prot = _load_test_protein(
        'alphafold/relax/testdata/multiple_disulfides_target.pdb'
        )
    viols, _ = amber_minimize.find_violations(prot)

    expected_between_residues_connection_mask = np.zeros((191,), np.float32)
    for residue in (42, 43, 59, 60, 135, 136):
      expected_between_residues_connection_mask[residue] = 1.0

    expected_clash_indices = np.array([
        [8, 4],
        [8, 5],
        [13, 3],
        [14, 1],
        [14, 4],
        [26, 4],
        [26, 5],
        [31, 8],
        [31, 10],
        [39, 0],
        [39, 1],
        [39, 2],
        [39, 3],
        [39, 4],
        [42, 5],
        [42, 6],
        [42, 7],
        [42, 8],
        [47, 7],
        [47, 8],
        [47, 9],
        [47, 10],
        [64, 4],
        [85, 5],
        [102, 4],
        [102, 5],
        [109, 13],
        [111, 5],
        [118, 6],
        [118, 7],
        [118, 8],
        [124, 4],
        [124, 5],
        [131, 5],
        [139, 7],
        [147, 4],
        [152, 7]], dtype=np.int32)
    expected_between_residues_clash_mask = np.zeros([191, 14])
    expected_between_residues_clash_mask[expected_clash_indices[:, 0],
                                         expected_clash_indices[:, 1]] += 1
    expected_per_atom_violations = np.zeros([191, 14])
    np.testing.assert_array_equal(
        viols['between_residues']['connections_per_residue_violation_mask'],
        expected_between_residues_connection_mask)
    np.testing.assert_array_equal(
        viols['between_residues']['clashes_per_atom_clash_mask'],
        expected_between_residues_clash_mask)
    np.testing.assert_array_equal(
        viols['within_residues']['per_atom_violations'],
        expected_per_atom_violations)


if __name__ == '__main__':
  absltest.main()
