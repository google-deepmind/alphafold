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

"""Unit tests for all_atom.frame_aligned_point_error."""

import numpy as np
from absl.testing import absltest
from absl.testing import parameterized

from alphafold.model import all_atom
from alphafold.model import r3

L1_CLAMP_DISTANCE = 10.0


def vecs_from_xyz(coords):
    """Convert XYZ coords to Vecs object."""
    coords = np.asarray(coords)
    return r3.vecs_from_tensor(coords)


def identity_rigid_batch(size):
    """Generate identity rigid frames of given batch size."""
    ones = np.ones(size)
    zeros = np.zeros(size)
    rot = r3.Rots(
        ones, zeros, zeros,
        zeros, ones, zeros,
        zeros, zeros, ones,
    )
    trans = r3.Vecs(zeros, zeros, zeros)
    return r3.Rigids(rot, trans)


def global_transform_rigid(angle_deg, translation, batch_dims):
    """Create rigid transform with rotation (about x) and translation."""
    angle_deg = np.asarray(angle_deg)
    translation = np.asarray(translation)

    for _ in range(batch_dims):
        angle_deg = np.expand_dims(angle_deg, 0)
        translation = np.expand_dims(translation, 0)

    sin = np.sin(np.deg2rad(angle_deg))
    cos = np.cos(np.deg2rad(angle_deg))
    ones = np.ones_like(sin)
    zeros = np.zeros_like(sin)

    rot = r3.Rots(
        ones, zeros, zeros,
        zeros, cos, -sin,
        zeros, sin, cos
    )
    trans = r3.Vecs(translation[..., 0], translation[..., 1], translation[..., 2])
    return r3.Rigids(rot, trans)


class AllAtomTest(parameterized.TestCase):

    @parameterized.named_parameters(
        ('no_transform', 0, [0, 0, 0]),
        ('rotation_90', 90, [0, 0, 0]),
        ('translation_10z', 0, [0, 0, 10]),
        ('rot_174_trans_xyz', 174, [1, 1, 1]),
    )
    def test_fape_zero_on_perfect_alignment(self, rot_angle, translation):
        """FAPE should be zero if prediction = target under global transform."""
        coords = np.array([
            [21.182, 23.095, 19.731],
            [22.055, 20.919, 17.294],
            [24.599, 20.005, 15.041],
            [25.567, 18.214, 12.166],
            [28.063, 17.082, 10.043],
            [28.779, 15.569, 6.985],
            [30.581, 13.815, 4.612],
            [29.258, 12.193, 2.296],
        ])

        # Target and predicted positions
        target_positions = vecs_from_xyz(coords)
        transform = global_transform_rigid(rot_angle, translation, batch_dims=1)
        pred_positions = r3.rigids_mul_vecs(transform, target_positions)

        # Identity frames, with global transform applied to pred_frames
        target_frames = identity_rigid_batch(10)
        pred_frames = r3.rigids_mul_rigids(transform, target_frames)

        # All atoms and frames are valid
        positions_mask = np.ones(target_positions.x.shape[0])
        frames_mask = np.ones(10)

        error = all_atom.frame_aligned_point_error(
            pred_frames, target_frames, frames_mask,
            pred_positions, target_positions, positions_mask,
            L1_CLAMP_DISTANCE, L1_CLAMP_DISTANCE, epsilon=0
        )

        self.assertAlmostEqual(error, 0.0)

    @parameterized.named_parameters(
        ('exact_match',
         [[0, 0, 0], [5, 0, 0], [10, 0, 0]],
         [[0, 0, 0], [5, 0, 0], [10, 0, 0]],
         0.0),
        ('small_shift',
         [[0, 0, 0], [5, 0, 0], [10, 0, 0]],
         [[2.5, 0, 0], [7.5, 0, 0], [7.5, 0, 0]],
         0.25),
        ('half_shift',
         [[0, 0, 0], [5, 0, 0], [10, 0, 0]],
         [[5, 0, 0], [10, 0, 0], [15, 0, 0]],
         0.5),
        ('large_shift',
         [[0, 0, 0], [5, 0, 0], [10, 0, 0]],
         [[10, 0, 0], [15, 0, 0], [0, 0, 0]],
         1.0),
    )
    def test_fape_matches_expected_value(self, target_coords, pred_coords, expected_fape):
        """Verify FAPE outputs correct normalized error for simple inputs."""
        target_positions = vecs_from_xyz(target_coords)
        pred_positions = vecs_from_xyz(pred_coords)

        positions_mask = np.ones(len(target_coords))
        target_frames = identity_rigid_batch(2)
        pred_frames = target_frames
        frames_mask = np.ones(2)

        error = all_atom.frame_aligned_point_error(
            pred_frames, target_frames, frames_mask,
            pred_positions, target_positions, positions_mask,
            L1_CLAMP_DISTANCE, L1_CLAMP_DISTANCE, epsilon=0
        )

        self.assertAlmostEqual(error, expected_fape)


if __name__ == '__main__':
    absltest.main()
