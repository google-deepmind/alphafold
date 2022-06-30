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
"""Ops for all atom representations."""

from typing import Dict, Text

from alphafold.common import residue_constants
from alphafold.model import geometry
from alphafold.model import utils
import jax
import jax.numpy as jnp
import numpy as np


def squared_difference(x, y):
  return jnp.square(x - y)


def _make_chi_atom_indices():
  """Returns atom indices needed to compute chi angles for all residue types.

  Returns:
    A tensor of shape [residue_types=21, chis=4, atoms=4]. The residue types are
    in the order specified in residue_constants.restypes + unknown residue type
    at the end. For chi angles which are not defined on the residue, the
    positions indices are by default set to 0.
  """
  chi_atom_indices = []
  for residue_name in residue_constants.restypes:
    residue_name = residue_constants.restype_1to3[residue_name]
    residue_chi_angles = residue_constants.chi_angles_atoms[residue_name]
    atom_indices = []
    for chi_angle in residue_chi_angles:
      atom_indices.append(
          [residue_constants.atom_order[atom] for atom in chi_angle])
    for _ in range(4 - len(atom_indices)):
      atom_indices.append([0, 0, 0, 0])  # For chi angles not defined on the AA.
    chi_atom_indices.append(atom_indices)

  chi_atom_indices.append([[0, 0, 0, 0]] * 4)  # For UNKNOWN residue.

  return np.array(chi_atom_indices)


def _make_renaming_matrices():
  """Matrices to map atoms to symmetry partners in ambiguous case."""
  # As the atom naming is ambiguous for 7 of the 20 amino acids, provide
  # alternative groundtruth coordinates where the naming is swapped
  restype_3 = [
      residue_constants.restype_1to3[res] for res in residue_constants.restypes
  ]
  restype_3 += ['UNK']
  # Matrices for renaming ambiguous atoms.
  all_matrices = {res: np.eye(14, dtype=np.float32) for res in restype_3}
  for resname, swap in residue_constants.residue_atom_renaming_swaps.items():
    correspondences = np.arange(14)
    for source_atom_swap, target_atom_swap in swap.items():
      source_index = residue_constants.restype_name_to_atom14_names[
          resname].index(source_atom_swap)
      target_index = residue_constants.restype_name_to_atom14_names[
          resname].index(target_atom_swap)
      correspondences[source_index] = target_index
      correspondences[target_index] = source_index
      renaming_matrix = np.zeros((14, 14), dtype=np.float32)
      for index, correspondence in enumerate(correspondences):
        renaming_matrix[index, correspondence] = 1.
    all_matrices[resname] = renaming_matrix.astype(np.float32)
  renaming_matrices = np.stack([all_matrices[restype] for restype in restype_3])
  return renaming_matrices


def _make_restype_atom37_mask():
  """Mask of which atoms are present for which residue type in atom37."""
  # create the corresponding mask
  restype_atom37_mask = np.zeros([21, 37], dtype=np.float32)
  for restype, restype_letter in enumerate(residue_constants.restypes):
    restype_name = residue_constants.restype_1to3[restype_letter]
    atom_names = residue_constants.residue_atoms[restype_name]
    for atom_name in atom_names:
      atom_type = residue_constants.atom_order[atom_name]
      restype_atom37_mask[restype, atom_type] = 1
  return restype_atom37_mask


def _make_restype_atom14_mask():
  """Mask of which atoms are present for which residue type in atom14."""
  restype_atom14_mask = []

  for rt in residue_constants.restypes:
    atom_names = residue_constants.restype_name_to_atom14_names[
        residue_constants.restype_1to3[rt]]
    restype_atom14_mask.append([(1. if name else 0.) for name in atom_names])

  restype_atom14_mask.append([0.] * 14)
  restype_atom14_mask = np.array(restype_atom14_mask, dtype=np.float32)
  return restype_atom14_mask


def _make_restype_atom37_to_atom14():
  """Map from atom37 to atom14 per residue type."""
  restype_atom37_to_atom14 = []  # mapping (restype, atom37) --> atom14
  for rt in residue_constants.restypes:
    atom_names = residue_constants.restype_name_to_atom14_names[
        residue_constants.restype_1to3[rt]]
    atom_name_to_idx14 = {name: i for i, name in enumerate(atom_names)}
    restype_atom37_to_atom14.append([
        (atom_name_to_idx14[name] if name in atom_name_to_idx14 else 0)
        for name in residue_constants.atom_types
    ])

  restype_atom37_to_atom14.append([0] * 37)
  restype_atom37_to_atom14 = np.array(restype_atom37_to_atom14, dtype=np.int32)
  return restype_atom37_to_atom14


def _make_restype_atom14_to_atom37():
  """Map from atom14 to atom37 per residue type."""
  restype_atom14_to_atom37 = []  # mapping (restype, atom14) --> atom37
  for rt in residue_constants.restypes:
    atom_names = residue_constants.restype_name_to_atom14_names[
        residue_constants.restype_1to3[rt]]
    restype_atom14_to_atom37.append([
        (residue_constants.atom_order[name] if name else 0)
        for name in atom_names
    ])
  # Add dummy mapping for restype 'UNK'
  restype_atom14_to_atom37.append([0] * 14)
  restype_atom14_to_atom37 = np.array(restype_atom14_to_atom37, dtype=np.int32)
  return restype_atom14_to_atom37


def _make_restype_atom14_is_ambiguous():
  """Mask which atoms are ambiguous in atom14."""
  # create an ambiguous atoms mask.  shape: (21, 14)
  restype_atom14_is_ambiguous = np.zeros((21, 14), dtype=np.float32)
  for resname, swap in residue_constants.residue_atom_renaming_swaps.items():
    for atom_name1, atom_name2 in swap.items():
      restype = residue_constants.restype_order[
          residue_constants.restype_3to1[resname]]
      atom_idx1 = residue_constants.restype_name_to_atom14_names[resname].index(
          atom_name1)
      atom_idx2 = residue_constants.restype_name_to_atom14_names[resname].index(
          atom_name2)
      restype_atom14_is_ambiguous[restype, atom_idx1] = 1
      restype_atom14_is_ambiguous[restype, atom_idx2] = 1

  return restype_atom14_is_ambiguous


def _make_restype_rigidgroup_base_atom37_idx():
  """Create Map from rigidgroups to atom37 indices."""
  # Create an array with the atom names.
  # shape (num_restypes, num_rigidgroups, 3_atoms): (21, 8, 3)
  base_atom_names = np.full([21, 8, 3], '', dtype=object)

  # 0: backbone frame
  base_atom_names[:, 0, :] = ['C', 'CA', 'N']

  # 3: 'psi-group'
  base_atom_names[:, 3, :] = ['CA', 'C', 'O']

  # 4,5,6,7: 'chi1,2,3,4-group'
  for restype, restype_letter in enumerate(residue_constants.restypes):
    resname = residue_constants.restype_1to3[restype_letter]
    for chi_idx in range(4):
      if residue_constants.chi_angles_mask[restype][chi_idx]:
        atom_names = residue_constants.chi_angles_atoms[resname][chi_idx]
        base_atom_names[restype, chi_idx + 4, :] = atom_names[1:]

  # Translate atom names into atom37 indices.
  lookuptable = residue_constants.atom_order.copy()
  lookuptable[''] = 0
  restype_rigidgroup_base_atom37_idx = np.vectorize(lambda x: lookuptable[x])(
      base_atom_names)
  return restype_rigidgroup_base_atom37_idx


CHI_ATOM_INDICES = _make_chi_atom_indices()
RENAMING_MATRICES = _make_renaming_matrices()
RESTYPE_ATOM14_TO_ATOM37 = _make_restype_atom14_to_atom37()
RESTYPE_ATOM37_TO_ATOM14 = _make_restype_atom37_to_atom14()
RESTYPE_ATOM37_MASK = _make_restype_atom37_mask()
RESTYPE_ATOM14_MASK = _make_restype_atom14_mask()
RESTYPE_ATOM14_IS_AMBIGUOUS = _make_restype_atom14_is_ambiguous()
RESTYPE_RIGIDGROUP_BASE_ATOM37_IDX = _make_restype_rigidgroup_base_atom37_idx()

# Create mask for existing rigid groups.
RESTYPE_RIGIDGROUP_MASK = np.zeros([21, 8], dtype=np.float32)
RESTYPE_RIGIDGROUP_MASK[:, 0] = 1
RESTYPE_RIGIDGROUP_MASK[:, 3] = 1
RESTYPE_RIGIDGROUP_MASK[:20, 4:] = residue_constants.chi_angles_mask


def get_atom37_mask(aatype):
  return utils.batched_gather(jnp.asarray(RESTYPE_ATOM37_MASK), aatype)


def get_atom14_mask(aatype):
  return utils.batched_gather(jnp.asarray(RESTYPE_ATOM14_MASK), aatype)


def get_atom14_is_ambiguous(aatype):
  return utils.batched_gather(jnp.asarray(RESTYPE_ATOM14_IS_AMBIGUOUS), aatype)


def get_atom14_to_atom37_map(aatype):
  return utils.batched_gather(jnp.asarray(RESTYPE_ATOM14_TO_ATOM37), aatype)


def get_atom37_to_atom14_map(aatype):
  return utils.batched_gather(jnp.asarray(RESTYPE_ATOM37_TO_ATOM14), aatype)


def atom14_to_atom37(atom14_data: jnp.ndarray,  # (N, 14, ...)
                     aatype: jnp.ndarray
                    ) -> jnp.ndarray:  # (N, 37, ...)
  """Convert atom14 to atom37 representation."""
  assert len(atom14_data.shape) in [2, 3]
  idx_atom37_to_atom14 = get_atom37_to_atom14_map(aatype)
  atom37_data = utils.batched_gather(
      atom14_data, idx_atom37_to_atom14, batch_dims=1)
  atom37_mask = get_atom37_mask(aatype)
  if len(atom14_data.shape) == 2:
    atom37_data *= atom37_mask
  elif len(atom14_data.shape) == 3:
    atom37_data *= atom37_mask[:, :, None].astype(atom37_data.dtype)
  return atom37_data


def atom37_to_atom14(aatype, all_atom_pos, all_atom_mask):
  """Convert Atom37 positions to Atom14 positions."""
  residx_atom14_to_atom37 = utils.batched_gather(
      jnp.asarray(RESTYPE_ATOM14_TO_ATOM37), aatype)
  atom14_mask = utils.batched_gather(
      all_atom_mask, residx_atom14_to_atom37, batch_dims=1).astype(jnp.float32)
  # create a mask for known groundtruth positions
  atom14_mask *= utils.batched_gather(jnp.asarray(RESTYPE_ATOM14_MASK), aatype)
  # gather the groundtruth positions
  atom14_positions = jax.tree_map(
      lambda x: utils.batched_gather(x, residx_atom14_to_atom37, batch_dims=1),
      all_atom_pos)
  atom14_positions = atom14_mask * atom14_positions
  return atom14_positions, atom14_mask


def get_alt_atom14(aatype, positions: geometry.Vec3Array, mask):
  """Get alternative atom14 positions."""
  # pick the transformation matrices for the given residue sequence
  # shape (num_res, 14, 14)
  renaming_transform = utils.batched_gather(
      jnp.asarray(RENAMING_MATRICES), aatype)

  alternative_positions = jax.tree_map(
      lambda x: jnp.sum(x, axis=1), positions[:, :, None] * renaming_transform)

  # Create the mask for the alternative ground truth (differs from the
  # ground truth mask, if only one of the atoms in an ambiguous pair has a
  # ground truth position)
  alternative_mask = jnp.sum(mask[..., None] * renaming_transform, axis=1)

  return alternative_positions, alternative_mask


def atom37_to_frames(
    aatype: jnp.ndarray,  # (...)
    all_atom_positions: geometry.Vec3Array,  # (..., 37)
    all_atom_mask: jnp.ndarray,  # (..., 37)
) -> Dict[Text, jnp.ndarray]:
  """Computes the frames for the up to 8 rigid groups for each residue."""
  # 0: 'backbone group',
  # 1: 'pre-omega-group', (empty)
  # 2: 'phi-group', (currently empty, because it defines only hydrogens)
  # 3: 'psi-group',
  # 4,5,6,7: 'chi1,2,3,4-group'
  aatype_in_shape = aatype.shape

  # If there is a batch axis, just flatten it away, and reshape everything
  # back at the end of the function.
  aatype = jnp.reshape(aatype, [-1])
  all_atom_positions = jax.tree_map(lambda x: jnp.reshape(x, [-1, 37]),
                                    all_atom_positions)
  all_atom_mask = jnp.reshape(all_atom_mask, [-1, 37])

  # Compute the gather indices for all residues in the chain.
  # shape (N, 8, 3)
  residx_rigidgroup_base_atom37_idx = utils.batched_gather(
      RESTYPE_RIGIDGROUP_BASE_ATOM37_IDX, aatype)

  # Gather the base atom positions for each rigid group.
  base_atom_pos = jax.tree_map(
      lambda x: utils.batched_gather(  # pylint: disable=g-long-lambda
          x, residx_rigidgroup_base_atom37_idx, batch_dims=1),
      all_atom_positions)

  # Compute the Rigids.
  point_on_neg_x_axis = base_atom_pos[:, :, 0]
  origin = base_atom_pos[:, :, 1]
  point_on_xy_plane = base_atom_pos[:, :, 2]
  gt_rotation = geometry.Rot3Array.from_two_vectors(
      origin - point_on_neg_x_axis, point_on_xy_plane - origin)

  gt_frames = geometry.Rigid3Array(gt_rotation, origin)

  # Compute a mask whether the group exists.
  # (N, 8)
  group_exists = utils.batched_gather(RESTYPE_RIGIDGROUP_MASK, aatype)

  # Compute a mask whether ground truth exists for the group
  gt_atoms_exist = utils.batched_gather(  # shape (N, 8, 3)
      all_atom_mask.astype(jnp.float32),
      residx_rigidgroup_base_atom37_idx,
      batch_dims=1)
  gt_exists = jnp.min(gt_atoms_exist, axis=-1) * group_exists  # (N, 8)

  # Adapt backbone frame to old convention (mirror x-axis and z-axis).
  rots = np.tile(np.eye(3, dtype=np.float32), [8, 1, 1])
  rots[0, 0, 0] = -1
  rots[0, 2, 2] = -1
  gt_frames = gt_frames.compose_rotation(
      geometry.Rot3Array.from_array(rots))

  # The frames for ambiguous rigid groups are just rotated by 180 degree around
  # the x-axis. The ambiguous group is always the last chi-group.
  restype_rigidgroup_is_ambiguous = np.zeros([21, 8], dtype=np.float32)
  restype_rigidgroup_rots = np.tile(np.eye(3, dtype=np.float32), [21, 8, 1, 1])

  for resname, _ in residue_constants.residue_atom_renaming_swaps.items():
    restype = residue_constants.restype_order[
        residue_constants.restype_3to1[resname]]
    chi_idx = int(sum(residue_constants.chi_angles_mask[restype]) - 1)
    restype_rigidgroup_is_ambiguous[restype, chi_idx + 4] = 1
    restype_rigidgroup_rots[restype, chi_idx + 4, 1, 1] = -1
    restype_rigidgroup_rots[restype, chi_idx + 4, 2, 2] = -1

  # Gather the ambiguity information for each residue.
  residx_rigidgroup_is_ambiguous = utils.batched_gather(
      restype_rigidgroup_is_ambiguous, aatype)
  ambiguity_rot = utils.batched_gather(restype_rigidgroup_rots, aatype)
  ambiguity_rot = geometry.Rot3Array.from_array(ambiguity_rot)

  # Create the alternative ground truth frames.
  alt_gt_frames = gt_frames.compose_rotation(ambiguity_rot)

  fix_shape = lambda x: jnp.reshape(x, aatype_in_shape + (8,))

  # reshape back to original residue layout
  gt_frames = jax.tree_map(fix_shape, gt_frames)
  gt_exists = fix_shape(gt_exists)
  group_exists = fix_shape(group_exists)
  residx_rigidgroup_is_ambiguous = fix_shape(residx_rigidgroup_is_ambiguous)
  alt_gt_frames = jax.tree_map(fix_shape, alt_gt_frames)

  return {
      'rigidgroups_gt_frames': gt_frames,  # Rigid (..., 8)
      'rigidgroups_gt_exists': gt_exists,  # (..., 8)
      'rigidgroups_group_exists': group_exists,  # (..., 8)
      'rigidgroups_group_is_ambiguous':
          residx_rigidgroup_is_ambiguous,  # (..., 8)
      'rigidgroups_alt_gt_frames': alt_gt_frames,  # Rigid (..., 8)
  }


def torsion_angles_to_frames(
    aatype: jnp.ndarray,  # (N)
    backb_to_global: geometry.Rigid3Array,  # (N)
    torsion_angles_sin_cos: jnp.ndarray  # (N, 7, 2)
) -> geometry.Rigid3Array:  # (N, 8)
  """Compute rigid group frames from torsion angles."""
  assert len(aatype.shape) == 1, (
      f'Expected array of rank 1, got array with shape: {aatype.shape}.')
  assert len(backb_to_global.rotation.shape) == 1, (
      f'Expected array of rank 1, got array with shape: '
      f'{backb_to_global.rotation.shape}')
  assert len(torsion_angles_sin_cos.shape) == 3, (
      f'Expected array of rank 3, got array with shape: '
      f'{torsion_angles_sin_cos.shape}')
  assert torsion_angles_sin_cos.shape[1] == 7, (
      f'wrong shape {torsion_angles_sin_cos.shape}')
  assert torsion_angles_sin_cos.shape[2] == 2, (
      f'wrong shape {torsion_angles_sin_cos.shape}')

  # Gather the default frames for all rigid groups.
  # geometry.Rigid3Array with shape (N, 8)
  m = utils.batched_gather(residue_constants.restype_rigid_group_default_frame,
                           aatype)
  default_frames = geometry.Rigid3Array.from_array4x4(m)

  # Create the rotation matrices according to the given angles (each frame is
  # defined such that its rotation is around the x-axis).
  sin_angles = torsion_angles_sin_cos[..., 0]
  cos_angles = torsion_angles_sin_cos[..., 1]

  # insert zero rotation for backbone group.
  num_residues, = aatype.shape
  sin_angles = jnp.concatenate([jnp.zeros([num_residues, 1]), sin_angles],
                               axis=-1)
  cos_angles = jnp.concatenate([jnp.ones([num_residues, 1]), cos_angles],
                               axis=-1)
  zeros = jnp.zeros_like(sin_angles)
  ones = jnp.ones_like(sin_angles)

  # all_rots are geometry.Rot3Array with shape (N, 8)
  all_rots = geometry.Rot3Array(ones, zeros, zeros,
                                zeros, cos_angles, -sin_angles,
                                zeros, sin_angles, cos_angles)

  # Apply rotations to the frames.
  all_frames = default_frames.compose_rotation(all_rots)

  # chi2, chi3, and chi4 frames do not transform to the backbone frame but to
  # the previous frame. So chain them up accordingly.

  chi1_frame_to_backb = all_frames[:, 4]
  chi2_frame_to_backb = chi1_frame_to_backb @ all_frames[:, 5]
  chi3_frame_to_backb = chi2_frame_to_backb @ all_frames[:, 6]
  chi4_frame_to_backb = chi3_frame_to_backb @ all_frames[:, 7]

  all_frames_to_backb = jax.tree_map(
      lambda *x: jnp.concatenate(x, axis=-1), all_frames[:, 0:5],
      chi2_frame_to_backb[:, None], chi3_frame_to_backb[:, None],
      chi4_frame_to_backb[:, None])

  # Create the global frames.
  # shape (N, 8)
  all_frames_to_global = backb_to_global[:, None] @ all_frames_to_backb

  return all_frames_to_global


def frames_and_literature_positions_to_atom14_pos(
    aatype: jnp.ndarray,  # (N)
    all_frames_to_global: geometry.Rigid3Array  # (N, 8)
) -> geometry.Vec3Array:  # (N, 14)
  """Put atom literature positions (atom14 encoding) in each rigid group."""

  # Pick the appropriate transform for every atom.
  residx_to_group_idx = utils.batched_gather(
      residue_constants.restype_atom14_to_rigid_group, aatype)
  group_mask = jax.nn.one_hot(
      residx_to_group_idx, num_classes=8)  # shape (N, 14, 8)

  # geometry.Rigid3Array with shape (N, 14)
  map_atoms_to_global = jax.tree_map(
      lambda x: jnp.sum(x[:, None, :] * group_mask, axis=-1),
      all_frames_to_global)

  # Gather the literature atom positions for each residue.
  # geometry.Vec3Array with shape (N, 14)
  lit_positions = geometry.Vec3Array.from_array(
      utils.batched_gather(
          residue_constants.restype_atom14_rigid_group_positions, aatype))

  # Transform each atom from its local frame to the global frame.
  # geometry.Vec3Array with shape (N, 14)
  pred_positions = map_atoms_to_global.apply_to_point(lit_positions)

  # Mask out non-existing atoms.
  mask = utils.batched_gather(residue_constants.restype_atom14_mask, aatype)
  pred_positions = pred_positions * mask

  return pred_positions


def extreme_ca_ca_distance_violations(
    positions: geometry.Vec3Array,  # (N, 37(14))
    mask: jnp.ndarray,  # (N, 37(14))
    residue_index: jnp.ndarray,  # (N)
    max_angstrom_tolerance=1.5
    ) -> jnp.ndarray:
  """Counts residues whose Ca is a large distance from its neighbor."""
  this_ca_pos = positions[:-1, 1]  # (N - 1,)
  this_ca_mask = mask[:-1, 1]         # (N - 1)
  next_ca_pos = positions[1:, 1]  # (N - 1,)
  next_ca_mask = mask[1:, 1]  # (N - 1)
  has_no_gap_mask = ((residue_index[1:] - residue_index[:-1]) == 1.0).astype(
      jnp.float32)
  ca_ca_distance = geometry.euclidean_distance(this_ca_pos, next_ca_pos, 1e-6)
  violations = (ca_ca_distance -
                residue_constants.ca_ca) > max_angstrom_tolerance
  mask = this_ca_mask * next_ca_mask * has_no_gap_mask
  return utils.mask_mean(mask=mask, value=violations)


def between_residue_bond_loss(
    pred_atom_positions: geometry.Vec3Array,  # (N, 37(14))
    pred_atom_mask: jnp.ndarray,  # (N, 37(14))
    residue_index: jnp.ndarray,  # (N)
    aatype: jnp.ndarray,  # (N)
    tolerance_factor_soft=12.0,
    tolerance_factor_hard=12.0) -> Dict[Text, jnp.ndarray]:
  """Flat-bottom loss to penalize structural violations between residues."""
  assert len(pred_atom_positions.shape) == 2
  assert len(pred_atom_mask.shape) == 2
  assert len(residue_index.shape) == 1
  assert len(aatype.shape) == 1

  # Get the positions of the relevant backbone atoms.
  this_ca_pos = pred_atom_positions[:-1, 1]  # (N - 1)
  this_ca_mask = pred_atom_mask[:-1, 1]         # (N - 1)
  this_c_pos = pred_atom_positions[:-1, 2]  # (N - 1)
  this_c_mask = pred_atom_mask[:-1, 2]          # (N - 1)
  next_n_pos = pred_atom_positions[1:, 0]  # (N - 1)
  next_n_mask = pred_atom_mask[1:, 0]           # (N - 1)
  next_ca_pos = pred_atom_positions[1:, 1]  # (N - 1)
  next_ca_mask = pred_atom_mask[1:, 1]          # (N - 1)
  has_no_gap_mask = ((residue_index[1:] - residue_index[:-1]) == 1.0).astype(
      jnp.float32)

  # Compute loss for the C--N bond.
  c_n_bond_length = geometry.euclidean_distance(this_c_pos, next_n_pos, 1e-6)

  # The C-N bond to proline has slightly different length because of the ring.
  next_is_proline = (
      aatype[1:] == residue_constants.restype_order['P']).astype(jnp.float32)
  gt_length = (
      (1. - next_is_proline) * residue_constants.between_res_bond_length_c_n[0]
      + next_is_proline * residue_constants.between_res_bond_length_c_n[1])
  gt_stddev = (
      (1. - next_is_proline) *
      residue_constants.between_res_bond_length_stddev_c_n[0] +
      next_is_proline * residue_constants.between_res_bond_length_stddev_c_n[1])
  c_n_bond_length_error = jnp.sqrt(1e-6 +
                                   jnp.square(c_n_bond_length - gt_length))
  c_n_loss_per_residue = jax.nn.relu(
      c_n_bond_length_error - tolerance_factor_soft * gt_stddev)
  mask = this_c_mask * next_n_mask * has_no_gap_mask
  c_n_loss = jnp.sum(mask * c_n_loss_per_residue) / (jnp.sum(mask) + 1e-6)
  c_n_violation_mask = mask * (
      c_n_bond_length_error > (tolerance_factor_hard * gt_stddev))

  # Compute loss for the angles.
  c_ca_unit_vec = (this_ca_pos - this_c_pos).normalized(1e-6)
  c_n_unit_vec = (next_n_pos - this_c_pos) / c_n_bond_length
  n_ca_unit_vec = (next_ca_pos - next_n_pos).normalized(1e-6)

  ca_c_n_cos_angle = c_ca_unit_vec.dot(c_n_unit_vec)
  gt_angle = residue_constants.between_res_cos_angles_ca_c_n[0]
  gt_stddev = residue_constants.between_res_bond_length_stddev_c_n[0]
  ca_c_n_cos_angle_error = jnp.sqrt(
      1e-6 + jnp.square(ca_c_n_cos_angle - gt_angle))
  ca_c_n_loss_per_residue = jax.nn.relu(
      ca_c_n_cos_angle_error - tolerance_factor_soft * gt_stddev)
  mask = this_ca_mask * this_c_mask * next_n_mask * has_no_gap_mask
  ca_c_n_loss = jnp.sum(mask * ca_c_n_loss_per_residue) / (jnp.sum(mask) + 1e-6)
  ca_c_n_violation_mask = mask * (ca_c_n_cos_angle_error >
                                  (tolerance_factor_hard * gt_stddev))

  c_n_ca_cos_angle = (-c_n_unit_vec).dot(n_ca_unit_vec)
  gt_angle = residue_constants.between_res_cos_angles_c_n_ca[0]
  gt_stddev = residue_constants.between_res_cos_angles_c_n_ca[1]
  c_n_ca_cos_angle_error = jnp.sqrt(
      1e-6 + jnp.square(c_n_ca_cos_angle - gt_angle))
  c_n_ca_loss_per_residue = jax.nn.relu(
      c_n_ca_cos_angle_error - tolerance_factor_soft * gt_stddev)
  mask = this_c_mask * next_n_mask * next_ca_mask * has_no_gap_mask
  c_n_ca_loss = jnp.sum(mask * c_n_ca_loss_per_residue) / (jnp.sum(mask) + 1e-6)
  c_n_ca_violation_mask = mask * (
      c_n_ca_cos_angle_error > (tolerance_factor_hard * gt_stddev))

  # Compute a per residue loss (equally distribute the loss to both
  # neighbouring residues).
  per_residue_loss_sum = (c_n_loss_per_residue +
                          ca_c_n_loss_per_residue +
                          c_n_ca_loss_per_residue)
  per_residue_loss_sum = 0.5 * (jnp.pad(per_residue_loss_sum, [[0, 1]]) +
                                jnp.pad(per_residue_loss_sum, [[1, 0]]))

  # Compute hard violations.
  violation_mask = jnp.max(
      jnp.stack([c_n_violation_mask,
                 ca_c_n_violation_mask,
                 c_n_ca_violation_mask]), axis=0)
  violation_mask = jnp.maximum(
      jnp.pad(violation_mask, [[0, 1]]),
      jnp.pad(violation_mask, [[1, 0]]))

  return {'c_n_loss_mean': c_n_loss,  # shape ()
          'ca_c_n_loss_mean': ca_c_n_loss,  # shape ()
          'c_n_ca_loss_mean': c_n_ca_loss,  # shape ()
          'per_residue_loss_sum': per_residue_loss_sum,  # shape (N)
          'per_residue_violation_mask': violation_mask  # shape (N)
         }


def between_residue_clash_loss(
    pred_positions: geometry.Vec3Array,  # (N, 14)
    atom_exists: jnp.ndarray,  # (N, 14)
    atom_radius: jnp.ndarray,  # (N, 14)
    residue_index: jnp.ndarray,  # (N)
    asym_id: jnp.ndarray,  # (N)
    overlap_tolerance_soft=1.5,
    overlap_tolerance_hard=1.5) -> Dict[Text, jnp.ndarray]:
  """Loss to penalize steric clashes between residues."""
  assert len(pred_positions.shape) == 2
  assert len(atom_exists.shape) == 2
  assert len(atom_radius.shape) == 2
  assert len(residue_index.shape) == 1

  # Create the distance matrix.
  # (N, N, 14, 14)
  dists = geometry.euclidean_distance(pred_positions[:, None, :, None],
                                      pred_positions[None, :, None, :], 1e-10)

  # Create the mask for valid distances.
  # shape (N, N, 14, 14)
  dists_mask = (atom_exists[:, None, :, None] * atom_exists[None, :, None, :])

  # Mask out all the duplicate entries in the lower triangular matrix.
  # Also mask out the diagonal (atom-pairs from the same residue) -- these atoms
  # are handled separately.
  dists_mask *= (
      residue_index[:, None, None, None] < residue_index[None, :, None, None])

  # Backbone C--N bond between subsequent residues is no clash.
  c_one_hot = jax.nn.one_hot(2, num_classes=14)
  n_one_hot = jax.nn.one_hot(0, num_classes=14)
  neighbour_mask = ((residue_index[:, None] + 1) == residue_index[None, :])
  neighbour_mask &= (asym_id[:, None] == asym_id[None, :])
  neighbour_mask = neighbour_mask[..., None, None]
  c_n_bonds = neighbour_mask * c_one_hot[None, None, :,
                                         None] * n_one_hot[None, None, None, :]
  dists_mask *= (1. - c_n_bonds)

  # Disulfide bridge between two cysteines is no clash.
  cys_sg_idx = residue_constants.restype_name_to_atom14_names['CYS'].index('SG')
  cys_sg_one_hot = jax.nn.one_hot(cys_sg_idx, num_classes=14)
  disulfide_bonds = (cys_sg_one_hot[None, None, :, None] *
                     cys_sg_one_hot[None, None, None, :])
  dists_mask *= (1. - disulfide_bonds)

  # Compute the lower bound for the allowed distances.
  # shape (N, N, 14, 14)
  dists_lower_bound = dists_mask * (
      atom_radius[:, None, :, None] + atom_radius[None, :, None, :])

  # Compute the error.
  # shape (N, N, 14, 14)
  dists_to_low_error = dists_mask * jax.nn.relu(
      dists_lower_bound - overlap_tolerance_soft - dists)

  # Compute the mean loss.
  # shape ()
  mean_loss = (jnp.sum(dists_to_low_error)
               / (1e-6 + jnp.sum(dists_mask)))

  # Compute the per atom loss sum.
  # shape (N, 14)
  per_atom_loss_sum = (jnp.sum(dists_to_low_error, axis=[0, 2]) +
                       jnp.sum(dists_to_low_error, axis=[1, 3]))

  # Compute the hard clash mask.
  # shape (N, N, 14, 14)
  clash_mask = dists_mask * (
      dists < (dists_lower_bound - overlap_tolerance_hard))

  # Compute the per atom clash.
  # shape (N, 14)
  per_atom_clash_mask = jnp.maximum(
      jnp.max(clash_mask, axis=[0, 2]),
      jnp.max(clash_mask, axis=[1, 3]))

  return {'mean_loss': mean_loss,  # shape ()
          'per_atom_loss_sum': per_atom_loss_sum,  # shape (N, 14)
          'per_atom_clash_mask': per_atom_clash_mask  # shape (N, 14)
         }


def within_residue_violations(
    pred_positions: geometry.Vec3Array,  # (N, 14)
    atom_exists: jnp.ndarray,  # (N, 14)
    dists_lower_bound: jnp.ndarray,  # (N, 14, 14)
    dists_upper_bound: jnp.ndarray,  # (N, 14, 14)
    tighten_bounds_for_loss=0.0,
) -> Dict[Text, jnp.ndarray]:
  """Find within-residue violations."""
  assert len(pred_positions.shape) == 2
  assert len(atom_exists.shape) == 2
  assert len(dists_lower_bound.shape) == 3
  assert len(dists_upper_bound.shape) == 3

  # Compute the mask for each residue.
  # shape (N, 14, 14)
  dists_masks = (1. - jnp.eye(14, 14)[None])
  dists_masks *= (atom_exists[:, :, None] * atom_exists[:, None, :])

  # Distance matrix
  # shape (N, 14, 14)
  dists = geometry.euclidean_distance(pred_positions[:, :, None],
                                      pred_positions[:, None, :], 1e-10)

  # Compute the loss.
  # shape (N, 14, 14)
  dists_to_low_error = jax.nn.relu(
      dists_lower_bound + tighten_bounds_for_loss - dists)
  dists_to_high_error = jax.nn.relu(
      dists + tighten_bounds_for_loss - dists_upper_bound)
  loss = dists_masks * (dists_to_low_error + dists_to_high_error)

  # Compute the per atom loss sum.
  # shape (N, 14)
  per_atom_loss_sum = (jnp.sum(loss, axis=1) +
                       jnp.sum(loss, axis=2))

  # Compute the violations mask.
  # shape (N, 14, 14)
  violations = dists_masks * ((dists < dists_lower_bound) |
                              (dists > dists_upper_bound))

  # Compute the per atom violations.
  # shape (N, 14)
  per_atom_violations = jnp.maximum(
      jnp.max(violations, axis=1), jnp.max(violations, axis=2))

  return {'per_atom_loss_sum': per_atom_loss_sum,  # shape (N, 14)
          'per_atom_violations': per_atom_violations  # shape (N, 14)
         }


def find_optimal_renaming(
    gt_positions: geometry.Vec3Array,  # (N, 14)
    alt_gt_positions: geometry.Vec3Array,  # (N, 14)
    atom_is_ambiguous: jnp.ndarray,  # (N, 14)
    gt_exists: jnp.ndarray,  # (N, 14)
    pred_positions: geometry.Vec3Array,  # (N, 14)
) -> jnp.ndarray:  # (N):
  """Find optimal renaming for ground truth that maximizes LDDT."""
  assert len(gt_positions.shape) == 2
  assert len(alt_gt_positions.shape) == 2
  assert len(atom_is_ambiguous.shape) == 2
  assert len(gt_exists.shape) == 2
  assert len(pred_positions.shape) == 2

  # Create the pred distance matrix.
  # shape (N, N, 14, 14)
  pred_dists = geometry.euclidean_distance(pred_positions[:, None, :, None],
                                           pred_positions[None, :, None, :],
                                           1e-10)

  # Compute distances for ground truth with original and alternative names.
  # shape (N, N, 14, 14)
  gt_dists = geometry.euclidean_distance(gt_positions[:, None, :, None],
                                         gt_positions[None, :, None, :], 1e-10)

  alt_gt_dists = geometry.euclidean_distance(alt_gt_positions[:, None, :, None],
                                             alt_gt_positions[None, :, None, :],
                                             1e-10)

  # Compute LDDT's.
  # shape (N, N, 14, 14)
  lddt = jnp.sqrt(1e-10 + squared_difference(pred_dists, gt_dists))
  alt_lddt = jnp.sqrt(1e-10 + squared_difference(pred_dists, alt_gt_dists))

  # Create a mask for ambiguous atoms in rows vs. non-ambiguous atoms
  # in cols.
  # shape (N ,N, 14, 14)
  mask = (
      gt_exists[:, None, :, None] *  # rows
      atom_is_ambiguous[:, None, :, None] *  # rows
      gt_exists[None, :, None, :] *  # cols
      (1. - atom_is_ambiguous[None, :, None, :]))  # cols

  # Aggregate distances for each residue to the non-amibuguous atoms.
  # shape (N)
  per_res_lddt = jnp.sum(mask * lddt, axis=[1, 2, 3])
  alt_per_res_lddt = jnp.sum(mask * alt_lddt, axis=[1, 2, 3])

  # Decide for each residue, whether alternative naming is better.
  # shape (N)
  alt_naming_is_better = (alt_per_res_lddt < per_res_lddt).astype(jnp.float32)

  return alt_naming_is_better  # shape (N)


def frame_aligned_point_error(
    pred_frames: geometry.Rigid3Array,  # shape (num_frames)
    target_frames: geometry.Rigid3Array,  # shape (num_frames)
    frames_mask: jnp.ndarray,  # shape (num_frames)
    pred_positions: geometry.Vec3Array,  # shape (num_positions)
    target_positions: geometry.Vec3Array,  # shape (num_positions)
    positions_mask: jnp.ndarray,  # shape (num_positions)
    pair_mask: jnp.ndarray,  # shape (num_frames, num_posiitons)
    l1_clamp_distance: float,
    length_scale=20.,
    epsilon=1e-4) -> jnp.ndarray:  # shape ()
  """Measure point error under different alignements.

  Computes error between two structures with B points
  under A alignments derived form the given pairs of frames.
  Args:
    pred_frames: num_frames reference frames for 'pred_positions'.
    target_frames: num_frames reference frames for 'target_positions'.
    frames_mask: Mask for frame pairs to use.
    pred_positions: num_positions predicted positions of the structure.
    target_positions: num_positions target positions of the structure.
    positions_mask: Mask on which positions to score.
    pair_mask: A (num_frames, num_positions) mask to use in the loss, useful
      for separating intra from inter chain losses.
    l1_clamp_distance: Distance cutoff on error beyond which gradients will
      be zero.
    length_scale: length scale to divide loss by.
    epsilon: small value used to regularize denominator for masked average.
  Returns:
    Masked Frame aligned point error.
  """
  # For now we do not allow any batch dimensions.
  assert len(pred_frames.rotation.shape) == 1
  assert len(target_frames.rotation.shape) == 1
  assert frames_mask.ndim == 1
  assert pred_positions.x.ndim == 1
  assert target_positions.x.ndim == 1
  assert positions_mask.ndim == 1

  # Compute array of predicted positions in the predicted frames.
  # geometry.Vec3Array (num_frames, num_positions)
  local_pred_pos = pred_frames[:, None].inverse().apply_to_point(
      pred_positions[None, :])

  # Compute array of target positions in the target frames.
  # geometry.Vec3Array (num_frames, num_positions)
  local_target_pos = target_frames[:, None].inverse().apply_to_point(
      target_positions[None, :])

  # Compute errors between the structures.
  # jnp.ndarray (num_frames, num_positions)
  error_dist = geometry.euclidean_distance(local_pred_pos, local_target_pos,
                                           epsilon)

  clipped_error_dist = jnp.clip(error_dist, 0, l1_clamp_distance)

  normed_error = clipped_error_dist / length_scale
  normed_error *= jnp.expand_dims(frames_mask, axis=-1)
  normed_error *= jnp.expand_dims(positions_mask, axis=-2)
  if pair_mask is not None:
    normed_error *= pair_mask

  mask = (jnp.expand_dims(frames_mask, axis=-1) *
          jnp.expand_dims(positions_mask, axis=-2))
  if pair_mask is not None:
    mask *= pair_mask
  normalization_factor = jnp.sum(mask, axis=(-1, -2))
  return (jnp.sum(normed_error, axis=(-2, -1)) /
          (epsilon + normalization_factor))


def get_chi_atom_indices():
  """Returns atom indices needed to compute chi angles for all residue types.

  Returns:
    A tensor of shape [residue_types=21, chis=4, atoms=4]. The residue types are
    in the order specified in residue_constants.restypes + unknown residue type
    at the end. For chi angles which are not defined on the residue, the
    positions indices are by default set to 0.
  """
  chi_atom_indices = []
  for residue_name in residue_constants.restypes:
    residue_name = residue_constants.restype_1to3[residue_name]
    residue_chi_angles = residue_constants.chi_angles_atoms[residue_name]
    atom_indices = []
    for chi_angle in residue_chi_angles:
      atom_indices.append(
          [residue_constants.atom_order[atom] for atom in chi_angle])
    for _ in range(4 - len(atom_indices)):
      atom_indices.append([0, 0, 0, 0])  # For chi angles not defined on the AA.
    chi_atom_indices.append(atom_indices)

  chi_atom_indices.append([[0, 0, 0, 0]] * 4)  # For UNKNOWN residue.

  return jnp.asarray(chi_atom_indices)


def compute_chi_angles(positions: geometry.Vec3Array,
                       mask: geometry.Vec3Array,
                       aatype: geometry.Vec3Array):
  """Computes the chi angles given all atom positions and the amino acid type.

  Args:
    positions: A Vec3Array of shape
      [num_res, residue_constants.atom_type_num], with positions of
      atoms needed to calculate chi angles. Supports up to 1 batch dimension.
    mask: An optional tensor of shape
      [num_res, residue_constants.atom_type_num] that masks which atom
      positions are set for each residue. If given, then the chi mask will be
      set to 1 for a chi angle only if the amino acid has that chi angle and all
      the chi atoms needed to calculate that chi angle are set. If not given
      (set to None), the chi mask will be set to 1 for a chi angle if the amino
      acid has that chi angle and whether the actual atoms needed to calculate
      it were set will be ignored.
    aatype: A tensor of shape [num_res] with amino acid type integer
      code (0 to 21). Supports up to 1 batch dimension.

  Returns:
    A tuple of tensors (chi_angles, mask), where both have shape
    [num_res, 4]. The mask masks out unused chi angles for amino acid
    types that have less than 4 chi angles. If atom_positions_mask is set, the
    chi mask will also mask out uncomputable chi angles.
  """

  # Don't assert on the num_res and batch dimensions as they might be unknown.
  assert positions.shape[-1] == residue_constants.atom_type_num
  assert mask.shape[-1] == residue_constants.atom_type_num

  # Compute the table of chi angle indices. Shape: [restypes, chis=4, atoms=4].
  chi_atom_indices = get_chi_atom_indices()
  # Select atoms to compute chis. Shape: [num_res, chis=4, atoms=4].
  atom_indices = utils.batched_gather(
      params=chi_atom_indices, indices=aatype, axis=0)
  # Gather atom positions. Shape: [num_res, chis=4, atoms=4, xyz=3].
  chi_angle_atoms = jax.tree_map(
      lambda x: utils.batched_gather(  # pylint: disable=g-long-lambda
          params=x, indices=atom_indices, axis=-1, batch_dims=1), positions)
  a, b, c, d = [chi_angle_atoms[..., i] for i in range(4)]

  chi_angles = geometry.dihedral_angle(a, b, c, d)

  # Copy the chi angle mask, add the UNKNOWN residue. Shape: [restypes, 4].
  chi_angles_mask = list(residue_constants.chi_angles_mask)
  chi_angles_mask.append([0.0, 0.0, 0.0, 0.0])
  chi_angles_mask = jnp.asarray(chi_angles_mask)
  # Compute the chi angle mask. Shape [num_res, chis=4].
  chi_mask = utils.batched_gather(params=chi_angles_mask, indices=aatype,
                                  axis=0)

  # The chi_mask is set to 1 only when all necessary chi angle atoms were set.
  # Gather the chi angle atoms mask. Shape: [num_res, chis=4, atoms=4].
  chi_angle_atoms_mask = utils.batched_gather(
      params=mask, indices=atom_indices, axis=-1, batch_dims=1)
  # Check if all 4 chi angle atoms were set. Shape: [num_res, chis=4].
  chi_angle_atoms_mask = jnp.prod(chi_angle_atoms_mask, axis=[-1])
  chi_mask = chi_mask * chi_angle_atoms_mask.astype(jnp.float32)

  return chi_angles, chi_mask


def make_transform_from_reference(
    a_xyz: geometry.Vec3Array,
    b_xyz: geometry.Vec3Array,
    c_xyz: geometry.Vec3Array) -> geometry.Rigid3Array:
  """Returns rotation and translation matrices to convert from reference.

  Note that this method does not take care of symmetries. If you provide the
  coordinates in the non-standard way, the A atom will end up in the negative
  y-axis rather than in the positive y-axis. You need to take care of such
  cases in your code.

  Args:
    a_xyz: A Vec3Array.
    b_xyz: A Vec3Array.
    c_xyz: A Vec3Array.

  Returns:
    A Rigid3Array which, when applied to coordinates in a canonicalized
    reference frame, will give coordinates approximately equal
    the original coordinates (in the global frame).
  """
  rotation = geometry.Rot3Array.from_two_vectors(c_xyz - b_xyz,
                                                 a_xyz - b_xyz)
  return geometry.Rigid3Array(rotation, b_xyz)
