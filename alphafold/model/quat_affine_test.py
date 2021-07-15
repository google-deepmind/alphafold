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

"""Tests for quat_affine."""

from absl import logging
from absl.testing import absltest
import jax
import jax.numpy as jnp
import numpy as np
from alphafold.model import quat_affine

VERBOSE = False
np.set_printoptions(precision=3, suppress=True)

r2t = quat_affine.rot_list_to_tensor
v2t = quat_affine.vec_list_to_tensor

q2r = lambda q: r2t(quat_affine.quat_to_rot(q))


class QuatAffineTest(absltest.TestCase):

  def _assert_check(self, to_check, tol=1e-5):
    for k, (correct, generated) in to_check.items():
      if VERBOSE:
        logging.info(k)
        logging.info('Correct %s', correct)
        logging.info('Predicted %s', generated)
      self.assertLess(np.max(np.abs(correct - generated)), tol)

  def test_conversion(self):
    quat = jnp.array([-2., 5., -1., 4.])

    rotation = jnp.array([
        [0.26087, 0.130435, 0.956522],
        [-0.565217, -0.782609, 0.26087],
        [0.782609, -0.608696, -0.130435]])

    translation = jnp.array([1., -3., 4.])
    point = jnp.array([0.7, 3.2, -2.9])

    a = quat_affine.QuatAffine(quat, translation, unstack_inputs=True)
    true_new_point = jnp.matmul(rotation, point[:, None])[:, 0] + translation

    self._assert_check({
        'rot': (rotation, r2t(a.rotation)),
        'trans': (translation, v2t(a.translation)),
        'point': (true_new_point,
                  v2t(a.apply_to_point(jnp.moveaxis(point, -1, 0)))),
        # Because of the double cover, we must be careful and compare rotations
        'quat': (q2r(a.quaternion),
                 q2r(quat_affine.rot_to_quat(a.rotation))),

    })

  def test_double_cover(self):
    """Test that -q is the same rotation as q."""
    rng = jax.random.PRNGKey(42)
    keys = jax.random.split(rng)
    q = jax.random.normal(keys[0], (2, 4))
    trans = jax.random.normal(keys[1], (2, 3))
    a1 = quat_affine.QuatAffine(q, trans, unstack_inputs=True)
    a2 = quat_affine.QuatAffine(-q, trans, unstack_inputs=True)

    self._assert_check({
        'rot': (r2t(a1.rotation),
                r2t(a2.rotation)),
        'trans': (v2t(a1.translation),
                  v2t(a2.translation)),
    })

  def test_homomorphism(self):
    rng = jax.random.PRNGKey(42)
    keys = jax.random.split(rng, 4)
    vec_q1 = jax.random.normal(keys[0], (2, 3))

    q1 = jnp.concatenate([
        jnp.ones_like(vec_q1)[:, :1],
        vec_q1], axis=-1)

    q2 = jax.random.normal(keys[1], (2, 4))
    t1 = jax.random.normal(keys[2], (2, 3))
    t2 = jax.random.normal(keys[3], (2, 3))

    a1 = quat_affine.QuatAffine(q1, t1, unstack_inputs=True)
    a2 = quat_affine.QuatAffine(q2, t2, unstack_inputs=True)
    a21 = a2.pre_compose(jnp.concatenate([vec_q1, t1], axis=-1))

    rng, key = jax.random.split(rng)
    x = jax.random.normal(key, (2, 3))
    new_x = a21.apply_to_point(jnp.moveaxis(x, -1, 0))
    new_x_apply2 = a2.apply_to_point(a1.apply_to_point(jnp.moveaxis(x, -1, 0)))

    self._assert_check({
        'quat': (q2r(quat_affine.quat_multiply(a2.quaternion, a1.quaternion)),
                 q2r(a21.quaternion)),
        'rot': (jnp.matmul(r2t(a2.rotation), r2t(a1.rotation)),
                r2t(a21.rotation)),
        'point': (v2t(new_x_apply2),
                  v2t(new_x)),
        'inverse': (x, v2t(a21.invert_point(new_x))),
    })

  def test_batching(self):
    """Test that affine applies batchwise."""
    rng = jax.random.PRNGKey(42)
    keys = jax.random.split(rng, 3)
    q = jax.random.uniform(keys[0], (5, 2, 4))
    t = jax.random.uniform(keys[1], (2, 3))
    x = jax.random.uniform(keys[2], (5, 1, 3))

    a = quat_affine.QuatAffine(q, t, unstack_inputs=True)
    y = v2t(a.apply_to_point(jnp.moveaxis(x, -1, 0)))

    y_list = []
    for i in range(5):
      for j in range(2):
        a_local = quat_affine.QuatAffine(q[i, j], t[j],
                                         unstack_inputs=True)
        y_local = v2t(a_local.apply_to_point(jnp.moveaxis(x[i, 0], -1, 0)))
        y_list.append(y_local)
    y_combine = jnp.reshape(jnp.stack(y_list, axis=0), (5, 2, 3))

    self._assert_check({
        'batch': (y_combine, y),
        'quat': (q2r(a.quaternion),
                 q2r(quat_affine.rot_to_quat(a.rotation))),
    })

  def assertAllClose(self, a, b, rtol=1e-06, atol=1e-06):
    self.assertTrue(np.allclose(a, b, rtol=rtol, atol=atol))

  def assertAllEqual(self, a, b):
    self.assertTrue(np.all(np.array(a) == np.array(b)))


if __name__ == '__main__':
  absltest.main()
