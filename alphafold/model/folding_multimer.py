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

"""Modules and utilities for the structure module in the multimer system."""

import functools
import numbers
from typing import Any, Dict, Iterable, Mapping, Optional, Tuple, Union

from alphafold.common import residue_constants
from alphafold.model import all_atom_multimer
from alphafold.model import common_modules
from alphafold.model import geometry
from alphafold.model import modules
from alphafold.model import prng
from alphafold.model import utils
from alphafold.model.geometry import utils as geometry_utils
import haiku as hk
import jax
import jax.numpy as jnp
import ml_collections
import numpy as np


EPSILON = 1e-8
Float = Union[float, jnp.ndarray]


def squared_difference(x: jnp.ndarray, y: jnp.ndarray) -> jnp.ndarray:
  """Computes Squared difference between two arrays."""
  return jnp.square(x - y)


def make_backbone_affine(
    positions: geometry.Vec3Array,
    mask: jnp.ndarray,
    aatype: jnp.ndarray,
    ) -> Tuple[geometry.Rigid3Array, jnp.ndarray]:
  """Make backbone Rigid3Array and mask."""
  del aatype
  a = residue_constants.atom_order['N']
  b = residue_constants.atom_order['CA']
  c = residue_constants.atom_order['C']

  rigid_mask = (mask[:, a] * mask[:, b] * mask[:, c]).astype(
      jnp.float32)

  rigid = all_atom_multimer.make_transform_from_reference(
      a_xyz=positions[:, a], b_xyz=positions[:, b], c_xyz=positions[:, c])

  return rigid, rigid_mask


class QuatRigid(hk.Module):
  """Module for projecting Rigids via a quaternion."""

  def __init__(self,
               global_config: ml_collections.ConfigDict,
               rigid_shape: Union[int, Iterable[int]] = tuple(),
               full_quat: bool = False,
               init: str = 'zeros',
               name: str = 'quat_rigid'):
    """Module projecting a Rigid Object.

    For this Module the Rotation is parametrized as a quaternion,
    If 'full_quat' is True a 4 vector is produced for the rotation which is
    normalized and treated as a quaternion.
    When 'full_quat' is False a 3 vector is produced and the 1st component of
    the quaternion is set to 1.

    Args:
      global_config: Global Config, used to set certain properties of underlying
        Linear module, see common_modules.Linear for details.
      rigid_shape: Shape of Rigids relative to shape of activations, e.g. when
        activations have shape (n,) and this is (m,) output will be (n, m)
      full_quat: Whether to parametrize rotation using full quaternion.
      init: initializer to use, see common_modules.Linear for details
      name: Name to use for module.
    """
    self.init = init
    self.global_config = global_config
    if isinstance(rigid_shape, int):
      self.rigid_shape = (rigid_shape,)
    else:
      self.rigid_shape = tuple(rigid_shape)
    self.full_quat = full_quat
    super(QuatRigid, self).__init__(name=name)

  def __call__(self, activations: jnp.ndarray) -> geometry.Rigid3Array:
    """Executes Module.

    This returns a set of rigid with the same shape as activations, projecting
    the channel dimension, rigid_shape controls the trailing dimensions.
    For example when activations is shape (12, 5) and rigid_shape is (3, 2)
    then the shape of the output rigids will be (12, 3, 2).
    This also supports passing in an empty tuple for rigid shape, in that case
    the example would produce a rigid of shape (12,).

    Args:
      activations: Activations to use for projection, shape [..., num_channel]
    Returns:
      Rigid transformations with shape [...] + rigid_shape
    """
    if self.full_quat:
      rigid_dim = 7
    else:
      rigid_dim = 6
    linear_dims = self.rigid_shape + (rigid_dim,)
    rigid_flat = common_modules.Linear(
        linear_dims,
        initializer=self.init,
        precision=jax.lax.Precision.HIGHEST,
        name='rigid')(
            activations)
    rigid_flat = geometry_utils.unstack(rigid_flat)
    if self.full_quat:
      qw, qx, qy, qz = rigid_flat[:4]
      translation = rigid_flat[4:]
    else:
      qx, qy, qz = rigid_flat[:3]
      qw = jnp.ones_like(qx)
      translation = rigid_flat[3:]
    rotation = geometry.Rot3Array.from_quaternion(
        qw, qx, qy, qz, normalize=True)
    translation = geometry.Vec3Array(*translation)
    return geometry.Rigid3Array(rotation, translation)


class PointProjection(hk.Module):
  """Given input reprensentation and frame produces points in global frame."""

  def __init__(self,
               num_points: Union[Iterable[int], int],
               global_config: ml_collections.ConfigDict,
               return_local_points: bool = False,
               name: str = 'point_projection'):
    """Constructs Linear Module.

    Args:
      num_points: number of points to project. Can be tuple when outputting
          multiple dimensions
      global_config: Global Config, passed through to underlying Linear
      return_local_points: Whether to return points in local frame as well.
      name: name of module, used for name scopes.
    """
    if isinstance(num_points, numbers.Integral):
      self.num_points = (num_points,)
    else:
      self.num_points = tuple(num_points)

    self.return_local_points = return_local_points

    self.global_config = global_config

    super().__init__(name=name)

  def __call__(
      self, activations: jnp.ndarray, rigids: geometry.Rigid3Array
  ) -> Union[geometry.Vec3Array, Tuple[geometry.Vec3Array, geometry.Vec3Array]]:
    output_shape = self.num_points
    output_shape = output_shape[:-1] + (3 * output_shape[-1],)
    points_local = common_modules.Linear(
        output_shape,
        precision=jax.lax.Precision.HIGHEST,
        name='point_projection')(
            activations)
    points_local = jnp.split(points_local, 3, axis=-1)
    points_local = geometry.Vec3Array(*points_local)
    rigids = rigids[(...,) + (None,) * len(output_shape)]
    points_global = rigids.apply_to_point(points_local)
    if self.return_local_points:
      return points_global, points_local
    else:
      return points_global


class InvariantPointAttention(hk.Module):
  """Invariant point attention module.

  The high-level idea is that this attention module works over a set of points
  and associated orientations in 3D space (e.g. protein residues).

  Each residue outputs a set of queries and keys as points in their local
  reference frame.  The attention is then defined as the euclidean distance
  between the queries and keys in the global frame.
  """

  def __init__(self,
               config: ml_collections.ConfigDict,
               global_config: ml_collections.ConfigDict,
               dist_epsilon: float = 1e-8,
               name: str = 'invariant_point_attention'):
    """Initialize.

    Args:
      config: iterative Fold Head Config
      global_config: Global Config of Model.
      dist_epsilon: Small value to avoid NaN in distance calculation.
      name: Sonnet name.
    """
    super().__init__(name=name)

    self._dist_epsilon = dist_epsilon
    self._zero_initialize_last = global_config.zero_init

    self.config = config

    self.global_config = global_config

  def __call__(
      self,
      inputs_1d: jnp.ndarray,
      inputs_2d: jnp.ndarray,
      mask: jnp.ndarray,
      rigid: geometry.Rigid3Array,
      ) -> jnp.ndarray:
    """Compute geometric aware attention.

    Given a set of query residues (defined by affines and associated scalar
    features), this function computes geometric aware attention between the
    query residues and target residues.

    The residues produce points in their local reference frame, which
    are converted into the global frame to get attention via euclidean distance.

    Equivalently the target residues produce points in their local frame to be
    used as attention values, which are converted into the query residues local
    frames.

    Args:
      inputs_1d: (N, C) 1D input embedding that is the basis for the
        scalar queries.
      inputs_2d: (N, M, C') 2D input embedding, used for biases values in the
        attention between query_inputs_1d and target_inputs_1d.
      mask: (N, 1) mask to indicate query_inputs_1d that participate in
        the attention.
      rigid: Rigid object describing the position and orientation of
        every element in query_inputs_1d.

    Returns:
      Transformation of the input embedding.
    """

    num_head = self.config.num_head

    attn_logits = 0.

    num_point_qk = self.config.num_point_qk
    # Each point pair (q, k) contributes Var [0.5 ||q||^2 - <q, k>] = 9 / 2
    point_variance = max(num_point_qk, 1) * 9. / 2
    point_weights = np.sqrt(1.0 / point_variance)

    # This is equivalent to jax.nn.softplus, but avoids a bug in the test...
    softplus = lambda x: jnp.logaddexp(x, jnp.zeros_like(x))
    raw_point_weights = hk.get_parameter(
        'trainable_point_weights',
        shape=[num_head],
        # softplus^{-1} (1)
        init=hk.initializers.Constant(np.log(np.exp(1.) - 1.)))

    # Trainable per-head weights for points.
    trainable_point_weights = softplus(raw_point_weights)
    point_weights *= trainable_point_weights
    q_point = PointProjection([num_head, num_point_qk],
                              self.global_config,
                              name='q_point_projection')(inputs_1d,
                                                         rigid)

    k_point = PointProjection([num_head, num_point_qk],
                              self.global_config,
                              name='k_point_projection')(inputs_1d,
                                                         rigid)

    dist2 = geometry.square_euclidean_distance(
        q_point[:, None, :, :], k_point[None, :, :, :], epsilon=0.)
    attn_qk_point = -0.5 * jnp.sum(point_weights[:, None] * dist2, axis=-1)
    attn_logits += attn_qk_point

    num_scalar_qk = self.config.num_scalar_qk
    # We assume that all queries and keys come iid from N(0, 1) distribution
    # and compute the variances of the attention logits.
    # Each scalar pair (q, k) contributes Var q*k = 1
    scalar_variance = max(num_scalar_qk, 1) * 1.
    scalar_weights = np.sqrt(1.0 / scalar_variance)
    q_scalar = common_modules.Linear([num_head, num_scalar_qk],
                                     use_bias=False,
                                     name='q_scalar_projection')(
                                         inputs_1d)

    k_scalar = common_modules.Linear([num_head, num_scalar_qk],
                                     use_bias=False,
                                     name='k_scalar_projection')(
                                         inputs_1d)
    q_scalar *= scalar_weights
    attn_logits += jnp.einsum('qhc,khc->qkh', q_scalar, k_scalar)

    attention_2d = common_modules.Linear(
        num_head, name='attention_2d')(inputs_2d)
    attn_logits += attention_2d

    mask_2d = mask * jnp.swapaxes(mask, -1, -2)
    attn_logits -= 1e5 * (1. - mask_2d[..., None])

    attn_logits *= np.sqrt(1. / 3)     # Normalize by number of logit terms (3)
    attn = jax.nn.softmax(attn_logits, axis=-2)

    num_scalar_v = self.config.num_scalar_v

    v_scalar = common_modules.Linear([num_head, num_scalar_v],
                                     use_bias=False,
                                     name='v_scalar_projection')(
                                         inputs_1d)

    # [num_query_residues, num_head, num_scalar_v]
    result_scalar = jnp.einsum('qkh, khc->qhc', attn, v_scalar)

    num_point_v = self.config.num_point_v
    v_point = PointProjection([num_head, num_point_v],
                              self.global_config,
                              name='v_point_projection')(inputs_1d,
                                                         rigid)

    result_point_global = jax.tree.map(
        lambda x: jnp.sum(attn[..., None] * x, axis=-3), v_point[None])

    # Features used in the linear output projection. Should have the size
    # [num_query_residues, ?]
    output_features = []
    num_query_residues, _ = inputs_1d.shape

    flat_shape = [num_query_residues, -1]

    result_scalar = jnp.reshape(result_scalar, flat_shape)
    output_features.append(result_scalar)

    result_point_global = jax.tree.map(lambda r: jnp.reshape(r, flat_shape),
                                       result_point_global)
    result_point_local = rigid[..., None].apply_inverse_to_point(
        result_point_global)
    output_features.extend(
        [result_point_local.x, result_point_local.y, result_point_local.z])

    point_norms = result_point_local.norm(self._dist_epsilon)
    output_features.append(point_norms)

    # Dimensions: h = heads, i and j = residues,
    # c = inputs_2d channels
    # Contraction happens over the second residue dimension, similarly to how
    # the usual attention is performed.
    result_attention_over_2d = jnp.einsum('ijh, ijc->ihc', attn, inputs_2d)
    output_features.append(jnp.reshape(result_attention_over_2d, flat_shape))

    final_init = 'zeros' if self._zero_initialize_last else 'linear'

    final_act = jnp.concatenate(output_features, axis=-1)

    return common_modules.Linear(
        self.config.num_channel,
        initializer=final_init,
        name='output_projection')(final_act)


class FoldIteration(hk.Module):
  """A single iteration of iterative folding.

  First, each residue attends to all residues using InvariantPointAttention.
  Then, we apply transition layers to update the hidden representations.
  Finally, we use the hidden representations to produce an update to the
  affine of each residue.
  """

  def __init__(self,
               config: ml_collections.ConfigDict,
               global_config: ml_collections.ConfigDict,
               name: str = 'fold_iteration'):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(
      self,
      activations: Mapping[str, Any],
      aatype: jnp.ndarray,
      sequence_mask: jnp.ndarray,
      update_rigid: bool,
      is_training: bool,
      initial_act: jnp.ndarray,
      safe_key: Optional[prng.SafeKey] = None,
      static_feat_2d: Optional[jnp.ndarray] = None,
  ) -> Tuple[Dict[str, Any], Dict[str, Any]]:

    c = self.config

    if safe_key is None:
      safe_key = prng.SafeKey(hk.next_rng_key())

    def safe_dropout_fn(tensor, safe_key):
      return modules.apply_dropout(
          tensor=tensor,
          safe_key=safe_key,
          rate=0.0 if self.global_config.deterministic else c.dropout,
          is_training=is_training)

    rigid = activations['rigid']

    act = activations['act']
    attention_module = InvariantPointAttention(
        self.config, self.global_config)
    # Attention
    act += attention_module(
        inputs_1d=act,
        inputs_2d=static_feat_2d,
        mask=sequence_mask,
        rigid=rigid)

    safe_key, *sub_keys = safe_key.split(3)
    sub_keys = iter(sub_keys)
    act = safe_dropout_fn(act, next(sub_keys))
    act = common_modules.LayerNorm(
        axis=-1,
        create_scale=True,
        create_offset=True,
        name='attention_layer_norm')(
            act)
    final_init = 'zeros' if self.global_config.zero_init else 'linear'

    # Transition
    input_act = act
    for i in range(c.num_layer_in_transition):
      init = 'relu' if i < c.num_layer_in_transition - 1 else final_init
      act = common_modules.Linear(
          c.num_channel,
          initializer=init,
          name='transition')(
              act)
      if i < c.num_layer_in_transition - 1:
        act = jax.nn.relu(act)
    act += input_act
    act = safe_dropout_fn(act, next(sub_keys))
    act = common_modules.LayerNorm(
        axis=-1,
        create_scale=True,
        create_offset=True,
        name='transition_layer_norm')(act)
    if update_rigid:
      # Rigid update
      rigid_update = QuatRigid(
          self.global_config, init=final_init)(
              act)
      rigid = rigid @ rigid_update

    sc = MultiRigidSidechain(c.sidechain, self.global_config)(
        rigid.scale_translation(c.position_scale), [act, initial_act], aatype)

    outputs = {'rigid': rigid, 'sc': sc}

    rotation = jax.tree.map(jax.lax.stop_gradient, rigid.rotation)
    rigid = geometry.Rigid3Array(rotation, rigid.translation)

    new_activations = {
        'act': act,
        'rigid': rigid
    }
    return new_activations, outputs


def generate_monomer_rigids(representations: Mapping[str, jnp.ndarray],
                            batch: Mapping[str, jnp.ndarray],
                            config: ml_collections.ConfigDict,
                            global_config: ml_collections.ConfigDict,
                            is_training: bool,
                            safe_key: prng.SafeKey
                            ) -> Dict[str, Any]:
  """Generate predicted Rigid's for a single chain.

  This is the main part of the iterative fold head - it iteratively applies
  folding to produce a set of predicted residue positions.

  Args:
    representations: Embeddings dictionary.
    batch: Batch dictionary.
    config: config for the iterative fold head.
    global_config: global config.
    is_training: is training.
    safe_key: A prng.SafeKey object that wraps a PRNG key.

  Returns:
    A dictionary containing residue Rigid's and sidechain positions.
  """
  c = config
  sequence_mask = batch['seq_mask'][:, None]
  act = common_modules.LayerNorm(
      axis=-1, create_scale=True, create_offset=True, name='single_layer_norm')(
          representations['single'])

  initial_act = act
  act = common_modules.Linear(
      c.num_channel, name='initial_projection')(act)

  # Sequence Mask has extra 1 at the end.
  rigid = geometry.Rigid3Array.identity(sequence_mask.shape[:-1])

  fold_iteration = FoldIteration(
      c, global_config, name='fold_iteration')

  assert len(batch['seq_mask'].shape) == 1

  activations = {
      'act':
          act,
      'rigid':
          rigid
  }

  act_2d = common_modules.LayerNorm(
      axis=-1,
      create_scale=True,
      create_offset=True,
      name='pair_layer_norm')(
          representations['pair'])

  safe_keys = safe_key.split(c.num_layer)
  outputs = []
  for key in safe_keys:

    activations, output = fold_iteration(
        activations,
        initial_act=initial_act,
        static_feat_2d=act_2d,
        aatype=batch['aatype'],
        safe_key=key,
        sequence_mask=sequence_mask,
        update_rigid=True,
        is_training=is_training,
        )
    outputs.append(output)

  output = jax.tree.map(lambda *x: jnp.stack(x), *outputs)
  # Pass along for LDDT-Head.
  output['act'] = activations['act']

  return output


class StructureModule(hk.Module):
  """StructureModule as a network head.

  Jumper et al. (2021) Suppl. Alg. 20 "StructureModule"
  """

  def __init__(self,
               config: ml_collections.ConfigDict,
               global_config: ml_collections.ConfigDict,
               name: str = 'structure_module'):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(self,
               representations: Mapping[str, jnp.ndarray],
               batch: Mapping[str, Any],
               is_training: bool,
               safe_key: Optional[prng.SafeKey] = None,
               compute_loss: bool = False
               ) -> Dict[str, Any]:
    c = self.config
    ret = {}

    if safe_key is None:
      safe_key = prng.SafeKey(hk.next_rng_key())

    output = generate_monomer_rigids(
        representations=representations,
        batch=batch,
        config=self.config,
        global_config=self.global_config,
        is_training=is_training,
        safe_key=safe_key)

    ret['traj'] = output['rigid'].scale_translation(c.position_scale).to_array()
    ret['sidechains'] = output['sc']
    ret['sidechains']['atom_pos'] = ret['sidechains']['atom_pos'].to_array()
    ret['sidechains']['frames'] = ret['sidechains']['frames'].to_array()
    if 'local_atom_pos' in ret['sidechains']:
      ret['sidechains']['local_atom_pos'] = ret['sidechains'][
          'local_atom_pos'].to_array()
      ret['sidechains']['local_frames'] = ret['sidechains'][
          'local_frames'].to_array()

    aatype = batch['aatype']
    seq_mask = batch['seq_mask']

    atom14_pred_mask = all_atom_multimer.get_atom14_mask(
        aatype) * seq_mask[:, None]
    atom14_pred_positions = output['sc']['atom_pos'][-1]
    ret['final_atom14_positions'] = atom14_pred_positions  # (N, 14, 3)
    ret['final_atom14_mask'] = atom14_pred_mask  # (N, 14)

    atom37_mask = all_atom_multimer.get_atom37_mask(aatype) * seq_mask[:, None]
    atom37_pred_positions = all_atom_multimer.atom14_to_atom37(
        atom14_pred_positions, aatype)
    atom37_pred_positions *= atom37_mask[:, :, None]
    ret['final_atom_positions'] = atom37_pred_positions  # (N, 37, 3)
    ret['final_atom_mask'] = atom37_mask  # (N, 37)
    ret['final_rigids'] = ret['traj'][-1]

    ret['act'] = output['act']

    if compute_loss:
      return ret
    else:
      no_loss_features = ['final_atom_positions', 'final_atom_mask', 'act']
      no_loss_ret = {k: ret[k] for k in no_loss_features}
      return no_loss_ret

  def loss(self,
           value: Mapping[str, Any],
           batch: Mapping[str, Any]
           ) -> Dict[str, Any]:

    raise NotImplementedError(
        'This function should be called on a batch with reordered chains (see '
        'Evans et al (2021) Section 7.3. Multi-Chain Permutation Alignment.')

    ret = {'loss': 0.}

    ret['metrics'] = {}

    aatype = batch['aatype']
    all_atom_positions = batch['all_atom_positions']
    all_atom_positions = geometry.Vec3Array.from_array(all_atom_positions)
    all_atom_mask = batch['all_atom_mask']
    seq_mask = batch['seq_mask']
    residue_index = batch['residue_index']

    gt_rigid, gt_affine_mask = make_backbone_affine(all_atom_positions,
                                                    all_atom_mask,
                                                    aatype)

    chi_angles, chi_mask = all_atom_multimer.compute_chi_angles(
        all_atom_positions, all_atom_mask, aatype)

    pred_mask = all_atom_multimer.get_atom14_mask(aatype)
    pred_mask *= seq_mask[:, None]
    pred_positions = value['final_atom14_positions']
    pred_positions = geometry.Vec3Array.from_array(pred_positions)

    gt_positions, gt_mask, alt_naming_is_better = compute_atom14_gt(
        aatype, all_atom_positions, all_atom_mask, pred_positions)

    violations = find_structural_violations(
        aatype=aatype,
        residue_index=residue_index,
        mask=pred_mask,
        pred_positions=pred_positions,
        config=self.config,
        asym_id=batch['asym_id'])

    sidechains = value['sidechains']

    gt_chi_angles = get_renamed_chi_angles(aatype, chi_angles,
                                           alt_naming_is_better)

    # Several violation metrics:
    violation_metrics = compute_violation_metrics(
        residue_index=residue_index,
        mask=pred_mask,
        seq_mask=seq_mask,
        pred_positions=pred_positions,
        violations=violations)
    ret['metrics'].update(violation_metrics)

    target_rigid = geometry.Rigid3Array.from_array(value['traj'])
    gt_frames_mask = gt_affine_mask

    # Split the loss into within-chain and between-chain components.
    intra_chain_mask = batch['asym_id'][:, None] == batch['asym_id'][None, :]
    intra_chain_bb_loss, intra_chain_fape = backbone_loss(
        gt_rigid=gt_rigid,
        gt_frames_mask=gt_frames_mask,
        gt_positions_mask=gt_affine_mask,
        target_rigid=target_rigid,
        config=self.config.intra_chain_fape,
        pair_mask=intra_chain_mask)
    interface_bb_loss, interface_fape = backbone_loss(
        gt_rigid=gt_rigid,
        gt_frames_mask=gt_frames_mask,
        gt_positions_mask=gt_affine_mask,
        target_rigid=target_rigid,
        config=self.config.interface_fape,
        pair_mask=1. - intra_chain_mask)

    bb_loss = intra_chain_bb_loss + interface_bb_loss
    ret['fape'] = intra_chain_fape + interface_fape
    ret['bb_loss'] = bb_loss
    ret['loss'] += bb_loss

    pred_frames = geometry.Rigid3Array.from_array(sidechains['frames'])
    pred_positions = geometry.Vec3Array.from_array(sidechains['atom_pos'])
    gt_sc_frames, gt_sc_frames_mask = compute_frames(
        aatype=aatype,
        all_atom_positions=all_atom_positions,
        all_atom_mask=all_atom_mask,
        use_alt=alt_naming_is_better)

    sc_loss = sidechain_loss(
        gt_frames=gt_sc_frames,
        gt_frames_mask=gt_sc_frames_mask,
        gt_positions=gt_positions,
        gt_mask=gt_mask,
        pred_frames=pred_frames,
        pred_positions=pred_positions,
        config=self.config)

    ret['loss'] = ((1 - self.config.sidechain.weight_frac) * ret['loss'] +
                   self.config.sidechain.weight_frac * sc_loss['loss'])
    ret['sidechain_fape'] = sc_loss['fape']

    unnormed_angles = sidechains['unnormalized_angles_sin_cos']
    pred_angles = sidechains['angles_sin_cos']

    sup_chi_loss, ret['chi_loss'], ret[
        'angle_norm_loss'] = supervised_chi_loss(
            sequence_mask=seq_mask,
            target_chi_mask=chi_mask,
            target_chi_angles=gt_chi_angles,
            aatype=aatype,
            pred_angles=pred_angles,
            unnormed_angles=unnormed_angles,
            config=self.config)
    ret['loss'] += sup_chi_loss

    if self.config.structural_violation_loss_weight:

      ret['loss'] += structural_violation_loss(
          mask=pred_mask, violations=violations, config=self.config)

    return ret


def compute_atom14_gt(
    aatype: jnp.ndarray,
    all_atom_positions: geometry.Vec3Array,
    all_atom_mask: jnp.ndarray,
    pred_pos: geometry.Vec3Array
) -> Tuple[geometry.Vec3Array, jnp.ndarray, jnp.ndarray]:
  """Find atom14 positions, this includes finding the correct renaming."""
  gt_positions, gt_mask = all_atom_multimer.atom37_to_atom14(
      aatype, all_atom_positions,
      all_atom_mask)
  alt_gt_positions, alt_gt_mask = all_atom_multimer.get_alt_atom14(
      aatype, gt_positions, gt_mask)
  atom_is_ambiguous = all_atom_multimer.get_atom14_is_ambiguous(aatype)

  alt_naming_is_better = all_atom_multimer.find_optimal_renaming(
      gt_positions=gt_positions,
      alt_gt_positions=alt_gt_positions,
      atom_is_ambiguous=atom_is_ambiguous,
      gt_exists=gt_mask,
      pred_positions=pred_pos)

  use_alt = alt_naming_is_better[:, None]

  gt_mask = (1. - use_alt) * gt_mask + use_alt * alt_gt_mask
  gt_positions = (1. - use_alt) * gt_positions + use_alt * alt_gt_positions

  return gt_positions, gt_mask, alt_naming_is_better


def backbone_loss(gt_rigid: geometry.Rigid3Array,
                  gt_frames_mask: jnp.ndarray,
                  gt_positions_mask: jnp.ndarray,
                  target_rigid: geometry.Rigid3Array,
                  config: ml_collections.ConfigDict,
                  pair_mask: jnp.ndarray
                  ) -> Tuple[Float, jnp.ndarray]:
  """Backbone FAPE Loss."""
  loss_fn = functools.partial(
      all_atom_multimer.frame_aligned_point_error,
      l1_clamp_distance=config.atom_clamp_distance,
      length_scale=config.loss_unit_distance)

  loss_fn = jax.vmap(loss_fn, (0, None, None, 0, None, None, None))
  fape = loss_fn(target_rigid, gt_rigid, gt_frames_mask,
                 target_rigid.translation, gt_rigid.translation,
                 gt_positions_mask, pair_mask)

  return jnp.mean(fape), fape[-1]


def compute_frames(
    aatype: jnp.ndarray,
    all_atom_positions: geometry.Vec3Array,
    all_atom_mask: jnp.ndarray,
    use_alt: jnp.ndarray
    ) -> Tuple[geometry.Rigid3Array, jnp.ndarray]:
  """Compute Frames from all atom positions.

  Args:
    aatype: array of aatypes, int of [N]
    all_atom_positions: Vector of all atom positions, shape [N, 37]
    all_atom_mask: mask, shape [N]
    use_alt: whether to use alternative orientation for ambiguous aatypes
             shape [N]
  Returns:
    Rigid corresponding to Frames w shape [N, 8],
    mask which Rigids are present w shape [N, 8]
  """
  frames_batch = all_atom_multimer.atom37_to_frames(aatype, all_atom_positions,
                                                    all_atom_mask)
  gt_frames = frames_batch['rigidgroups_gt_frames']
  alt_gt_frames = frames_batch['rigidgroups_alt_gt_frames']
  use_alt = use_alt[:, None]

  renamed_gt_frames = jax.tree.map(
      lambda x, y: (1. - use_alt) * x + use_alt * y, gt_frames, alt_gt_frames)

  return renamed_gt_frames, frames_batch['rigidgroups_gt_exists']


def sidechain_loss(gt_frames: geometry.Rigid3Array,
                   gt_frames_mask: jnp.ndarray,
                   gt_positions: geometry.Vec3Array,
                   gt_mask: jnp.ndarray,
                   pred_frames: geometry.Rigid3Array,
                   pred_positions: geometry.Vec3Array,
                   config: ml_collections.ConfigDict
                   ) -> Dict[str, jnp.ndarray]:
  """Sidechain Loss using cleaned up rigids."""

  flat_gt_frames = jax.tree.map(jnp.ravel, gt_frames)
  flat_frames_mask = jnp.ravel(gt_frames_mask)

  flat_gt_positions = jax.tree.map(jnp.ravel, gt_positions)
  flat_positions_mask = jnp.ravel(gt_mask)

  # Compute frame_aligned_point_error score for the final layer.
  def _slice_last_layer_and_flatten(x):
    return jnp.ravel(x[-1])

  flat_pred_frames = jax.tree.map(_slice_last_layer_and_flatten, pred_frames)
  flat_pred_positions = jax.tree.map(_slice_last_layer_and_flatten,
                                     pred_positions)
  fape = all_atom_multimer.frame_aligned_point_error(
      pred_frames=flat_pred_frames,
      target_frames=flat_gt_frames,
      frames_mask=flat_frames_mask,
      pred_positions=flat_pred_positions,
      target_positions=flat_gt_positions,
      positions_mask=flat_positions_mask,
      pair_mask=None,
      length_scale=config.sidechain.loss_unit_distance,
      l1_clamp_distance=config.sidechain.atom_clamp_distance)

  return {
      'fape': fape,
      'loss': fape}


def structural_violation_loss(mask: jnp.ndarray,
                              violations: Mapping[str, Float],
                              config: ml_collections.ConfigDict
                              ) -> Float:
  """Computes Loss for structural Violations."""
  # Put all violation losses together to one large loss.
  num_atoms = jnp.sum(mask).astype(jnp.float32) + 1e-6
  between_residues = violations['between_residues']
  within_residues = violations['within_residues']
  return (config.structural_violation_loss_weight *
          (between_residues['bonds_c_n_loss_mean'] +
           between_residues['angles_ca_c_n_loss_mean']  +
           between_residues['angles_c_n_ca_loss_mean'] +
           jnp.sum(between_residues['clashes_per_atom_loss_sum'] +
                   within_residues['per_atom_loss_sum']) / num_atoms
           ))


def find_structural_violations(
    aatype: jnp.ndarray,
    residue_index: jnp.ndarray,
    mask: jnp.ndarray,
    pred_positions: geometry.Vec3Array,  # (N, 14)
    config: ml_collections.ConfigDict,
    asym_id: jnp.ndarray,
    ) -> Dict[str, Any]:
  """Computes several checks for structural Violations."""

  # Compute between residue backbone violations of bonds and angles.
  connection_violations = all_atom_multimer.between_residue_bond_loss(
      pred_atom_positions=pred_positions,
      pred_atom_mask=mask.astype(jnp.float32),
      residue_index=residue_index.astype(jnp.float32),
      aatype=aatype,
      tolerance_factor_soft=config.violation_tolerance_factor,
      tolerance_factor_hard=config.violation_tolerance_factor)

  # Compute the van der Waals radius for every atom
  # (the first letter of the atom name is the element type).
  # shape (N, 14)
  atomtype_radius = jnp.array([
      residue_constants.van_der_waals_radius[name[0]]
      for name in residue_constants.atom_types
  ])
  residx_atom14_to_atom37 = all_atom_multimer.get_atom14_to_atom37_map(aatype)
  atom_radius = mask * utils.batched_gather(atomtype_radius,
                                            residx_atom14_to_atom37)

  # Compute the between residue clash loss.
  between_residue_clashes = all_atom_multimer.between_residue_clash_loss(
      pred_positions=pred_positions,
      atom_exists=mask,
      atom_radius=atom_radius,
      residue_index=residue_index,
      overlap_tolerance_soft=config.clash_overlap_tolerance,
      overlap_tolerance_hard=config.clash_overlap_tolerance,
      asym_id=asym_id)

  # Compute all within-residue violations (clashes,
  # bond length and angle violations).
  restype_atom14_bounds = residue_constants.make_atom14_dists_bounds(
      overlap_tolerance=config.clash_overlap_tolerance,
      bond_length_tolerance_factor=config.violation_tolerance_factor)
  dists_lower_bound = utils.batched_gather(restype_atom14_bounds['lower_bound'],
                                           aatype)
  dists_upper_bound = utils.batched_gather(restype_atom14_bounds['upper_bound'],
                                           aatype)
  within_residue_violations = all_atom_multimer.within_residue_violations(
      pred_positions=pred_positions,
      atom_exists=mask,
      dists_lower_bound=dists_lower_bound,
      dists_upper_bound=dists_upper_bound,
      tighten_bounds_for_loss=0.0)

  # Combine them to a single per-residue violation mask (used later for LDDT).
  per_residue_violations_mask = jnp.max(jnp.stack([
      connection_violations['per_residue_violation_mask'],
      jnp.max(between_residue_clashes['per_atom_clash_mask'], axis=-1),
      jnp.max(within_residue_violations['per_atom_violations'],
              axis=-1)]), axis=0)

  return {
      'between_residues': {
          'bonds_c_n_loss_mean':
              connection_violations['c_n_loss_mean'],  # ()
          'angles_ca_c_n_loss_mean':
              connection_violations['ca_c_n_loss_mean'],  # ()
          'angles_c_n_ca_loss_mean':
              connection_violations['c_n_ca_loss_mean'],  # ()
          'connections_per_residue_loss_sum':
              connection_violations['per_residue_loss_sum'],  # (N)
          'connections_per_residue_violation_mask':
              connection_violations['per_residue_violation_mask'],  # (N)
          'clashes_mean_loss':
              between_residue_clashes['mean_loss'],  # ()
          'clashes_per_atom_loss_sum':
              between_residue_clashes['per_atom_loss_sum'],  # (N, 14)
          'clashes_per_atom_clash_mask':
              between_residue_clashes['per_atom_clash_mask'],  # (N, 14)
      },
      'within_residues': {
          'per_atom_loss_sum':
              within_residue_violations['per_atom_loss_sum'],  # (N, 14)
          'per_atom_violations':
              within_residue_violations['per_atom_violations'],  # (N, 14),
      },
      'total_per_residue_violations_mask':
          per_residue_violations_mask,  # (N)
  }


def compute_violation_metrics(
    residue_index: jnp.ndarray,
    mask: jnp.ndarray,
    seq_mask: jnp.ndarray,
    pred_positions: geometry.Vec3Array,  # (N, 14)
    violations: Mapping[str, jnp.ndarray],
) -> Dict[str, jnp.ndarray]:
  """Compute several metrics to assess the structural violations."""
  ret = {}
  between_residues = violations['between_residues']
  within_residues = violations['within_residues']
  extreme_ca_ca_violations = all_atom_multimer.extreme_ca_ca_distance_violations(
      positions=pred_positions,
      mask=mask.astype(jnp.float32),
      residue_index=residue_index.astype(jnp.float32))
  ret['violations_extreme_ca_ca_distance'] = extreme_ca_ca_violations
  ret['violations_between_residue_bond'] = utils.mask_mean(
      mask=seq_mask,
      value=between_residues['connections_per_residue_violation_mask'])
  ret['violations_between_residue_clash'] = utils.mask_mean(
      mask=seq_mask,
      value=jnp.max(between_residues['clashes_per_atom_clash_mask'], axis=-1))
  ret['violations_within_residue'] = utils.mask_mean(
      mask=seq_mask,
      value=jnp.max(within_residues['per_atom_violations'], axis=-1))
  ret['violations_per_residue'] = utils.mask_mean(
      mask=seq_mask, value=violations['total_per_residue_violations_mask'])
  return ret


def supervised_chi_loss(
    sequence_mask: jnp.ndarray,
    target_chi_mask: jnp.ndarray,
    aatype: jnp.ndarray,
    target_chi_angles: jnp.ndarray,
    pred_angles: jnp.ndarray,
    unnormed_angles: jnp.ndarray,
    config: ml_collections.ConfigDict) -> Tuple[Float, Float, Float]:
  """Computes loss for direct chi angle supervision."""
  eps = 1e-6
  chi_mask = target_chi_mask.astype(jnp.float32)

  pred_angles = pred_angles[:, :, 3:]

  residue_type_one_hot = jax.nn.one_hot(
      aatype, residue_constants.restype_num + 1, dtype=jnp.float32)[None]
  chi_pi_periodic = jnp.einsum('ijk, kl->ijl', residue_type_one_hot,
                               jnp.asarray(residue_constants.chi_pi_periodic))

  true_chi = target_chi_angles[None]
  sin_true_chi = jnp.sin(true_chi)
  cos_true_chi = jnp.cos(true_chi)
  sin_cos_true_chi = jnp.stack([sin_true_chi, cos_true_chi], axis=-1)

  # This is -1 if chi is pi periodic and +1 if it's 2 pi periodic
  shifted_mask = (1 - 2 * chi_pi_periodic)[..., None]
  sin_cos_true_chi_shifted = shifted_mask * sin_cos_true_chi

  sq_chi_error = jnp.sum(
      squared_difference(sin_cos_true_chi, pred_angles), -1)
  sq_chi_error_shifted = jnp.sum(
      squared_difference(sin_cos_true_chi_shifted, pred_angles), -1)
  sq_chi_error = jnp.minimum(sq_chi_error, sq_chi_error_shifted)

  sq_chi_loss = utils.mask_mean(mask=chi_mask[None], value=sq_chi_error)
  angle_norm = jnp.sqrt(jnp.sum(jnp.square(unnormed_angles), axis=-1) + eps)
  norm_error = jnp.abs(angle_norm - 1.)
  angle_norm_loss = utils.mask_mean(mask=sequence_mask[None, :, None],
                                    value=norm_error)
  loss = (config.chi_weight * sq_chi_loss
          + config.angle_norm_weight * angle_norm_loss)
  return loss, sq_chi_loss, angle_norm_loss


def l2_normalize(x: jnp.ndarray,
                 axis: int = -1,
                 epsilon: float = 1e-12
                 ) -> jnp.ndarray:
  return x / jnp.sqrt(
      jnp.maximum(jnp.sum(x**2, axis=axis, keepdims=True), epsilon))


def get_renamed_chi_angles(aatype: jnp.ndarray,
                           chi_angles: jnp.ndarray,
                           alt_is_better: jnp.ndarray
                           ) -> jnp.ndarray:
  """Return renamed chi angles."""
  chi_angle_is_ambiguous = utils.batched_gather(
      jnp.array(residue_constants.chi_pi_periodic, dtype=jnp.float32), aatype)
  alt_chi_angles = chi_angles + np.pi * chi_angle_is_ambiguous
  # Map back to [-pi, pi].
  alt_chi_angles = alt_chi_angles - 2 * np.pi * (alt_chi_angles > np.pi).astype(
      jnp.float32)
  alt_is_better = alt_is_better[:, None]
  return (1. - alt_is_better) * chi_angles + alt_is_better * alt_chi_angles


class MultiRigidSidechain(hk.Module):
  """Class to make side chain atoms."""

  def __init__(self,
               config: ml_collections.ConfigDict,
               global_config: ml_collections.ConfigDict,
               name: str = 'rigid_sidechain'):
    super().__init__(name=name)
    self.config = config
    self.global_config = global_config

  def __call__(self,
               rigid: geometry.Rigid3Array,
               representations_list: Iterable[jnp.ndarray],
               aatype: jnp.ndarray
               ) -> Dict[str, Any]:
    """Predict sidechains using multi-rigid representations.

    Args:
      rigid: The Rigid's for each residue (translations in angstoms)
      representations_list: A list of activations to predict sidechains from.
      aatype: amino acid types.

    Returns:
      dict containing atom positions and frames (in angstrom)
    """
    act = [
        common_modules.Linear(  # pylint: disable=g-complex-comprehension
            self.config.num_channel,
            name='input_projection')(jax.nn.relu(x))
        for x in representations_list]
    # Sum the activation list (equivalent to concat then Conv1D)
    act = sum(act)

    final_init = 'zeros' if self.global_config.zero_init else 'linear'

    # Mapping with some residual blocks.
    for _ in range(self.config.num_residual_block):
      old_act = act
      act = common_modules.Linear(
          self.config.num_channel,
          initializer='relu',
          name='resblock1')(
              jax.nn.relu(act))
      act = common_modules.Linear(
          self.config.num_channel,
          initializer=final_init,
          name='resblock2')(
              jax.nn.relu(act))
      act += old_act

    # Map activations to torsion angles.
    # [batch_size, num_res, 14]
    num_res = act.shape[0]
    unnormalized_angles = common_modules.Linear(
        14, name='unnormalized_angles')(
            jax.nn.relu(act))
    unnormalized_angles = jnp.reshape(
        unnormalized_angles, [num_res, 7, 2])
    angles = l2_normalize(unnormalized_angles, axis=-1)

    outputs = {
        'angles_sin_cos': angles,  # jnp.ndarray (N, 7, 2)
        'unnormalized_angles_sin_cos':
            unnormalized_angles,  # jnp.ndarray (N, 7, 2)
    }

    # Map torsion angles to frames.
    # geometry.Rigid3Array with shape (N, 8)
    all_frames_to_global = all_atom_multimer.torsion_angles_to_frames(
        aatype,
        rigid,
        angles)

    # Use frames and literature positions to create the final atom coordinates.
    # geometry.Vec3Array with shape (N, 14)
    pred_positions = all_atom_multimer.frames_and_literature_positions_to_atom14_pos(
        aatype, all_frames_to_global)

    outputs.update({
        'atom_pos': pred_positions,  # geometry.Vec3Array (N, 14)
        'frames': all_frames_to_global,  # geometry.Rigid3Array (N, 8)
    })
    return outputs
