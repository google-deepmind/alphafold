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

"""Tests for shape_helpers."""

import numpy as np
import tensorflow.compat.v1 as tf

from alphafold.model.tf import shape_helpers


class ShapeTest(tf.test.TestCase):

  def test_shape_list(self):
    """Test that shape_list can allow for reshaping to dynamic shapes."""
    a = tf.zeros([10, 4, 4, 2])
    p = tf.placeholder(tf.float32, shape=[None, None, 1, 4, 4])
    shape_dyn = shape_helpers.shape_list(p)[:2] + [4, 4]

    b = tf.reshape(a, shape_dyn)
    with self.session() as sess:
      out = sess.run(b, feed_dict={p: np.ones((20, 1, 1, 4, 4))})

    self.assertAllEqual(out.shape, (20, 1, 4, 4))


if __name__ == '__main__':
  tf.disable_v2_behavior()
  tf.test.main()
