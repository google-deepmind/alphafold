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

"""Tests for all_atom."""

from absl.testing import absltest
from absl.testing import parameterized
import numpy as np
from alphafold.model import all_atom
from alphafold.model import r3

L1_CLAMP_DISTANCE = 10


def get_identity_rigid(shape):
  """Returns identity rigid transform."""

  ones = np.ones(shape)
  zeros = np.zeros(shape)
  rot = r3.Rots(ones, zeros, zeros,
                zeros, ones, zeros,
                zeros, zeros, ones)
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
  rot = r3.Rots(ones, zeros, zeros,
                zeros, cos_angle, -sin_angle,
                zeros, sin_angle, cos_angle)
  trans = r3.Vecs(translation[..., 0], translation[..., 1], translation[..., 2])
  return r3.Rigids(rot, trans)


class AllAtomTest(parameterized.TestCase, absltest.TestCase):

  @parameterized.named_parameters(
      ('identity', 0, [0, 0, 0]),
      ('rot_90', 90, [0, 0, 0]),
      ('trans_10', 0, [0, 0, 10]),
      ('rot_174_trans_1', 174, [1, 1, 1]))
  def test_frame_aligned_point_error_perfect_on_global_transform(
      self, rot_angle, translation):
    """Tests global transform between target and preds gives perfect score."""

    # pylint: disable=bad-whitespace
    target_positions = np.array(
        [[ 21.182,  23.095,  19.731],
         [ 22.055,  20.919,  17.294],
         [ 24.599,  20.005,  15.041],
         [ 25.567,  18.214,  12.166],
         [ 28.063,  17.082,  10.043],
         [ 28.779,  15.569,   6.985],
         [ 30.581,  13.815,   4.612],
         [ 29.258,  12.193,   2.296]])
    # pylint: enable=bad-whitespace
    global_rigid_transform = get_global_rigid_transform(
        rot_angle, translation, 1)

    target_positions = r3.vecs_from_tensor(target_positions)
    pred_positions = r3.rigids_mul_vecs(
        global_rigid_transform, target_positions)
    positions_mask = np.ones(target_positions.x.shape[0])

    target_frames = get_identity_rigid(10)
    pred_frames = r3.rigids_mul_rigids(global_rigid_transform, target_frames)
    frames_mask = np.ones(10)

    fape = all_atom.frame_aligned_point_error(
        pred_frames, target_frames, frames_mask, pred_positions,
        target_positions, positions_mask, L1_CLAMP_DISTANCE,
        L1_CLAMP_DISTANCE, epsilon=0)
    self.assertAlmostEqual(fape, 0.)

  @parameterized.named_parameters(
      ('identity',
       [[0, 0, 0], [5, 0, 0], [10, 0, 0]],
       [[0, 0, 0], [5, 0, 0], [10, 0, 0]],
       0.),
      ('shift_2.5',
       [[0, 0, 0], [5, 0, 0], [10, 0, 0]],
       [[2.5, 0, 0], [7.5, 0, 0], [7.5, 0, 0]],
       0.25),
      ('shift_5',
       [[0, 0, 0], [5, 0, 0], [10, 0, 0]],
       [[5, 0, 0], [10, 0, 0], [15, 0, 0]],
       0.5),
      ('shift_10',
       [[0, 0, 0], [5, 0, 0], [10, 0, 0]],
       [[10, 0, 0], [15, 0, 0], [0, 0, 0]],
       1.))
  def test_frame_aligned_point_error_matches_expected(
      self, target_positions, pred_positions, expected_alddt):
    """Tests score matches expected."""

    target_frames = get_identity_rigid(2)
    pred_frames = target_frames
    frames_mask = np.ones(2)

    target_positions = r3.vecs_from_tensor(np.array(target_positions))
    pred_positions = r3.vecs_from_tensor(np.array(pred_positions))
    positions_mask = np.ones(target_positions.x.shape[0])

    alddt = all_atom.frame_aligned_point_error(
        pred_frames, target_frames, frames_mask, pred_positions,
        target_positions, positions_mask, L1_CLAMP_DISTANCE,
        L1_CLAMP_DISTANCE, epsilon=0)
    self.assertAlmostEqual(alddt, expected_alddt)


if __name__ == '__main__':
  absltest.main()
