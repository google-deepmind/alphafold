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

from absl.testing import absltest
from absl.testing import parameterized
from alphafold.common import residue_constants
from alphafold.model import all_atom
from alphafold.model import r3
import jax
import jax.numpy as jnp
import numpy as np


L1_CLAMP_DISTANCE = 10

BL_C_N = residue_constants.between_res_bond_length_c_n
BL_STD_DEV_C_N = residue_constants.between_res_bond_length_stddev_c_n
COS_CA_C_N = residue_constants.between_res_cos_angles_ca_c_n
COS_C_N_CA = residue_constants.between_res_cos_angles_c_n_ca


def _relu(x):
  """Computes relu on a numpy array."""
  return np.maximum(0, x)


def _get_positions_for_ca_c_n_violation_mask():
  p = np.zeros((2, 37, 3), dtype=np.float32)
  p[1, 0, 0] = BL_C_N[0]
  return p


def _get_mask_for_ca_c_n_violation_mask():
  m = np.ones((2, 37), dtype=np.float32)
  m[1, 1] = 0.0
  return m


def get_identity_rigid(shape):
  """Returns identity rigid transform."""

  ones = np.ones(shape)
  zeros = np.zeros(shape)
  rot = r3.Rots(ones, zeros, zeros, zeros, ones, zeros, zeros, zeros, ones)
  trans = r3.Vecs(zeros, zeros, zeros)
  return r3.Rigids(rot, trans)


def get_global_rigid_transform(rot_angle, translation, bcast_dims):
  """Returns rigid transform that globally rotates/translates by same amount."""

  rot_angle = np.asarray(rot_angle)
  translation = np.asarray(translation)
  if bcast_dims:
    for _ in range(bcast_dims):
      rot_angle = np.expand_dims(rot_angle, 0)
      translation = np.expand_dims(translation, 0)
  sin_angle = np.sin(np.deg2rad(rot_angle))
  cos_angle = np.cos(np.deg2rad(rot_angle))
  ones = np.ones_like(sin_angle)
  zeros = np.zeros_like(sin_angle)
  rot = r3.Rots(
      ones,
      zeros,
      zeros,
      zeros,
      cos_angle,
      -sin_angle,
      zeros,
      sin_angle,
      cos_angle,
  )
  trans = r3.Vecs(translation[..., 0], translation[..., 1], translation[..., 2])
  return r3.Rigids(rot, trans)


class AllAtomTest(parameterized.TestCase):

  @parameterized.named_parameters(
      ('identity', 0, [0, 0, 0]),
      ('rot_90', 90, [0, 0, 0]),
      ('trans_10', 0, [0, 0, 10]),
      ('rot_174_trans_1', 174, [1, 1, 1]),
  )
  def test_frame_aligned_point_error_perfect_on_global_transform(
      self, rot_angle, translation
  ):
    """Tests global transform between target and preds gives perfect score."""

    # pylint: disable=bad-whitespace
    target_positions = np.array([
        [21.182, 23.095, 19.731],
        [22.055, 20.919, 17.294],
        [24.599, 20.005, 15.041],
        [25.567, 18.214, 12.166],
        [28.063, 17.082, 10.043],
        [28.779, 15.569, 6.985],
        [30.581, 13.815, 4.612],
        [29.258, 12.193, 2.296],
    ])
    # pylint: enable=bad-whitespace
    global_rigid_transform = get_global_rigid_transform(
        rot_angle, translation, 1
    )

    target_positions = r3.vecs_from_tensor(jax.numpy.array(target_positions))
    pred_positions = r3.rigids_mul_vecs(
        global_rigid_transform, target_positions
    )
    positions_mask = np.ones(target_positions.x.shape[0])

    target_frames = get_identity_rigid(10)
    pred_frames = r3.rigids_mul_rigids(global_rigid_transform, target_frames)
    frames_mask = np.ones(10)

    fape = all_atom.frame_aligned_point_error(
        pred_frames,
        target_frames,
        frames_mask,
        pred_positions,
        target_positions,
        positions_mask,
        L1_CLAMP_DISTANCE,
        L1_CLAMP_DISTANCE,
        epsilon=0,
    )
    self.assertAlmostEqual(fape, 0.0, places=6)

  @parameterized.named_parameters(
      (
          'identity',
          [[0, 0, 0], [5, 0, 0], [10, 0, 0]],
          [[0, 0, 0], [5, 0, 0], [10, 0, 0]],
          0.0,
      ),
      (
          'shift_2.5',
          [[0, 0, 0], [5, 0, 0], [10, 0, 0]],
          [[2.5, 0, 0], [7.5, 0, 0], [7.5, 0, 0]],
          0.25,
      ),
      (
          'shift_5',
          [[0, 0, 0], [5, 0, 0], [10, 0, 0]],
          [[5, 0, 0], [10, 0, 0], [15, 0, 0]],
          0.5,
      ),
      (
          'shift_10',
          [[0, 0, 0], [5, 0, 0], [10, 0, 0]],
          [[10, 0, 0], [15, 0, 0], [0, 0, 0]],
          1.0,
      ),
  )
  def test_frame_aligned_point_error_matches_expected(
      self, target_positions, pred_positions, expected_alddt
  ):
    """Tests score matches expected."""

    target_frames = get_identity_rigid(2)
    pred_frames = target_frames
    frames_mask = np.ones(2)

    target_positions = r3.vecs_from_tensor(jax.numpy.array(target_positions))
    pred_positions = r3.vecs_from_tensor(jax.numpy.array(pred_positions))
    positions_mask = np.ones(target_positions.x.shape[0])

    alddt = all_atom.frame_aligned_point_error(
        pred_frames,
        target_frames,
        frames_mask,
        pred_positions,
        target_positions,
        positions_mask,
        L1_CLAMP_DISTANCE,
        L1_CLAMP_DISTANCE,
        epsilon=0,
    )
    self.assertAlmostEqual(alddt, expected_alddt)

  @parameterized.named_parameters(
      dict(
          testcase_name='c_n_loss',
          key='c_n_loss_mean',
          pred_atom_positions=np.zeros((2, 37, 3), dtype=np.float32),
          pred_atom_mask=np.ones((2, 37), dtype=np.float32),
          residue_index=np.arange(2, dtype=np.int32),
          aatype=np.zeros(2, dtype=np.int32),
          expected_val=np.sum(
              _relu(
                  np.sqrt(1e-6 + np.square(0.001 - BL_C_N[0]))
                  - 12.0 * BL_STD_DEV_C_N[0]
              )
          ).astype(np.float32),
      ),
      dict(
          testcase_name='ca_c_n_loss',
          key='ca_c_n_loss_mean',
          pred_atom_positions=np.zeros((2, 37, 3), dtype=np.float32),
          pred_atom_mask=np.ones((2, 37), dtype=np.float32),
          residue_index=np.arange(2, dtype=np.int32),
          aatype=np.zeros(2, dtype=np.int32),
          expected_val=np.sum(
              _relu(
                  np.sqrt(1e-6 + np.square(-COS_CA_C_N[0]))
                  - 12.0 * COS_CA_C_N[1]
              )
          ).astype(np.float32),
      ),
      dict(
          testcase_name='c_n_ca_loss',
          key='c_n_ca_loss_mean',
          pred_atom_positions=np.zeros((2, 37, 3), dtype=np.float32),
          pred_atom_mask=np.ones((2, 37), dtype=np.float32),
          residue_index=np.arange(2, dtype=np.int32),
          aatype=np.zeros(2, dtype=np.int32),
          expected_val=np.sum(
              _relu(
                  np.sqrt(1e-6 + np.square(0.0 - COS_C_N_CA[0]))
                  - 12.0 * COS_C_N_CA[1]
              )
          ).astype(np.float32),
      ),
      dict(
          testcase_name='per_residue_loss_sum',
          key='per_residue_loss_sum',
          pred_atom_positions=np.zeros((2, 37, 3), dtype=np.float32),
          pred_atom_mask=np.ones((2, 37), dtype=np.float32),
          residue_index=np.arange(2, dtype=np.int32),
          aatype=np.zeros(2, dtype=np.int32),
          expected_val=np.array([0.665401, 0.665401], dtype=np.float32),
      ),
      dict(
          testcase_name='per_residue_violation_mask',
          key='per_residue_violation_mask',
          pred_atom_positions=np.zeros((2, 37, 3), dtype=np.float32),
          pred_atom_mask=np.ones((2, 37), dtype=np.float32),
          residue_index=np.arange(2, dtype=np.int32),
          aatype=np.zeros(2, dtype=np.int32),
          expected_val=np.array([1.0, 1.0], dtype=np.float32),
      ),
      dict(
          # This test verifies that the violation mask is correctly computed
          # for CA, C, N violations.
          testcase_name='ca_c_n_violation_mask',
          key='per_residue_violation_mask',
          pred_atom_positions=_get_positions_for_ca_c_n_violation_mask(),
          pred_atom_mask=_get_mask_for_ca_c_n_violation_mask(),
          residue_index=np.arange(2, dtype=np.int32),
          aatype=np.zeros(2, dtype=np.int32),
          expected_val=np.array([0.0, 0.0], dtype=np.float32),
          tolerance_factor_hard=15.0,
      ),
  )
  def test_between_residue_bond_loss(
      self,
      key,
      pred_atom_positions,
      pred_atom_mask,
      residue_index,
      aatype,
      expected_val,
      tolerance_factor_hard=12.0,
  ):
    got = all_atom.between_residue_bond_loss(
        pred_atom_positions=jnp.array(pred_atom_positions),
        pred_atom_mask=jnp.array(pred_atom_mask),
        residue_index=jnp.array(residue_index),
        aatype=jnp.array(aatype),
        tolerance_factor_hard=tolerance_factor_hard,
    )
    self.assertIn(key, got)
    self.assertEqual(
        got[key].shape,
        expected_val.shape,
        f'Shape mismatch for key "{key}"',
    )
    self.assertEqual(
        got[key].dtype,
        expected_val.dtype,
        f'Dtype mismatch for key "{key}"',
    )
    np.testing.assert_allclose(got[key], expected_val, rtol=2e-6)


if __name__ == '__main__':
  absltest.main()
