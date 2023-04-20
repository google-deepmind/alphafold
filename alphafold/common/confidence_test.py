# Copyright 2023 DeepMind Technologies Limited
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

"""Test confidence metrics."""


from absl.testing import absltest
from alphafold.common import confidence
import numpy as np


class ConfidenceTest(absltest.TestCase):

  def test_pae_json(self):
    pae = np.array([[0.01, 13.12345], [20.0987, 0.0]])
    pae_json = confidence.pae_json(pae=pae, max_pae=31.75)
    self.assertEqual(
        pae_json, '[{"predicted_aligned_error":[[0.0,13.1],[20.1,0.0]],'
        '"max_predicted_aligned_error":31.75}]')

  def test_confidence_json(self):
    plddt = np.array([42, 42.42])

    confidence_json = confidence.confidence_json(plddt=plddt)

    print(confidence_json)

    self.assertEqual(
        confidence_json,
        ('{"residueNumber":[1,2],'
         '"confidenceScore":[42.0,42.42],'
         '"confidenceCategory":["D","D"]}'),
    )


if __name__ == '__main__':
  absltest.main()
