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

"""Tests for lddt."""

from absl.testing import absltest
from absl.testing import parameterized
import numpy as np
from alphafold.model import lddt


class LddtTest(parameterized.TestCase, absltest.TestCase):

  @parameterized.named_parameters(
      ('same',
       [[0, 0, 0], [5, 0, 0], [10, 0, 0]],
       [[0, 0, 0], [5, 0, 0], [10, 0, 0]],
       [1, 1, 1]),
      ('all_shifted',
       [[0, 0, 0], [5, 0, 0], [10, 0, 0]],
       [[-1, 0, 0], [4, 0, 0], [9, 0, 0]],
       [1, 1, 1]),
      ('all_rotated',
       [[0, 0, 0], [5, 0, 0], [10, 0, 0]],
       [[0, 0, 0], [0, 5, 0], [0, 10, 0]],
       [1, 1, 1]),
      ('half_a_dist',
       [[0, 0, 0], [5, 0, 0]],
       [[0, 0, 0], [5.5-1e-5, 0, 0]],
       [1, 1]),
      ('one_a_dist',
       [[0, 0, 0], [5, 0, 0]],
       [[0, 0, 0], [6-1e-5, 0, 0]],
       [0.75, 0.75]),
      ('two_a_dist',
       [[0, 0, 0], [5, 0, 0]],
       [[0, 0, 0], [7-1e-5, 0, 0]],
       [0.5, 0.5]),
      ('four_a_dist',
       [[0, 0, 0], [5, 0, 0]],
       [[0, 0, 0], [9-1e-5, 0, 0]],
       [0.25, 0.25],),
      ('five_a_dist',
       [[0, 0, 0], [16-1e-5, 0, 0]],
       [[0, 0, 0], [11, 0, 0]],
       [0, 0]),
      ('no_pairs',
       [[0, 0, 0], [20, 0, 0]],
       [[0, 0, 0], [25-1e-5, 0, 0]],
       [1, 1]),
  )
  def test_lddt(
      self, predicted_pos, true_pos, exp_lddt):
    predicted_pos = np.array([predicted_pos], dtype=np.float32)
    true_points_mask = np.array([[[1]] * len(true_pos)], dtype=np.float32)
    true_pos = np.array([true_pos], dtype=np.float32)
    cutoff = 15.0
    per_residue = True

    result = lddt.lddt(
        predicted_pos, true_pos, true_points_mask, cutoff,
        per_residue)

    np.testing.assert_almost_equal(result, [exp_lddt], decimal=4)


if __name__ == '__main__':
  absltest.main()
