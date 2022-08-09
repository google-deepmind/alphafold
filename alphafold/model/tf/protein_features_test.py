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

"""Tests for protein_features."""
import uuid

from absl.testing import absltest
from absl.testing import parameterized
from alphafold.model.tf import protein_features
import tensorflow.compat.v1 as tf


def _random_bytes():
  return str(uuid.uuid4()).encode('utf-8')


class FeaturesTest(parameterized.TestCase, tf.test.TestCase):

  def setUp(self):
    super().setUp()
    tf.disable_v2_behavior()

  def testFeatureNames(self):
    self.assertEqual(len(protein_features.FEATURE_SIZES),
                     len(protein_features.FEATURE_TYPES))
    sorted_size_names = sorted(protein_features.FEATURE_SIZES.keys())
    sorted_type_names = sorted(protein_features.FEATURE_TYPES.keys())
    for i, size_name in enumerate(sorted_size_names):
      self.assertEqual(size_name, sorted_type_names[i])

  def testReplacement(self):
    for name in protein_features.FEATURE_SIZES.keys():
      sizes = protein_features.shape(name,
                                     num_residues=12,
                                     msa_length=24,
                                     num_templates=3)
      for x in sizes:
        self.assertEqual(type(x), int)
        self.assertGreater(x, 0)


if __name__ == '__main__':
  absltest.main()
