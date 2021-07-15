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

"""Data for AlphaFold."""
import numpy as np
import tensorflow.compat.v1 as tf

from alphafold.common import residue_constants
from alphafold.model.tf import shape_helpers
from alphafold.model.tf import shape_placeholders
from alphafold.model.tf import utils

# Pylint gets confused by the curry1 decorator because it changes the number
#   of arguments to the function.
# pylint:disable=no-value-for-parameter


NUM_RES = shape_placeholders.NUM_RES
NUM_MSA_SEQ = shape_placeholders.NUM_MSA_SEQ
NUM_EXTRA_SEQ = shape_placeholders.NUM_EXTRA_SEQ
NUM_TEMPLATES = shape_placeholders.NUM_TEMPLATES


def cast_64bit_ints(protein):

  for k, v in protein.items():
    if v.dtype == tf.int64:
      protein[k] = tf.cast(v, tf.int32)
  return protein


_MSA_FEATURE_NAMES = [
    'msa', 'deletion_matrix', 'msa_mask', 'msa_row_mask', 'bert_mask',
    'true_msa'
]


def make_seq_mask(protein):
  protein['seq_mask'] = tf.ones(
      shape_helpers.shape_list(protein['aatype']), dtype=tf.float32)
  return protein


def make_template_mask(protein):
  protein['template_mask'] = tf.ones(
      shape_helpers.shape_list(protein['template_domain_names']),
      dtype=tf.float32)
  return protein


def curry1(f):
  """Supply all arguments but the first."""

  def fc(*args, **kwargs):
    return lambda x: f(x, *args, **kwargs)

  return fc


@curry1
def add_distillation_flag(protein, distillation):
  protein['is_distillation'] = tf.constant(float(distillation),
                                           shape=[],
                                           dtype=tf.float32)
  return protein


def make_all_atom_aatype(protein):
  protein['all_atom_aatype'] = protein['aatype']
  return protein


def fix_templates_aatype(protein):
  """Fixes aatype encoding of templates."""
  # Map one-hot to indices.
  protein['template_aatype'] = tf.argmax(
      protein['template_aatype'], output_type=tf.int32, axis=-1)
  # Map hhsearch-aatype to our aatype.
  new_order_list = residue_constants.MAP_HHBLITS_AATYPE_TO_OUR_AATYPE
  new_order = tf.constant(new_order_list, dtype=tf.int32)
  protein['template_aatype'] = tf.gather(params=new_order,
                                         indices=protein['template_aatype'])
  return protein


def correct_msa_restypes(protein):
  """Correct MSA restype to have the same order as residue_constants."""
  new_order_list = residue_constants.MAP_HHBLITS_AATYPE_TO_OUR_AATYPE
  new_order = tf.constant(new_order_list, dtype=protein['msa'].dtype)
  protein['msa'] = tf.gather(new_order, protein['msa'], axis=0)

  perm_matrix = np.zeros((22, 22), dtype=np.float32)
  perm_matrix[range(len(new_order_list)), new_order_list] = 1.

  for k in protein:
    if 'profile' in k:  # Include both hhblits and psiblast profiles
      num_dim = protein[k].shape.as_list()[-1]
      assert num_dim in [20, 21, 22], (
          'num_dim for %s out of expected range: %s' % (k, num_dim))
      protein[k] = tf.tensordot(protein[k], perm_matrix[:num_dim, :num_dim], 1)
  return protein


def squeeze_features(protein):
  """Remove singleton and repeated dimensions in protein features."""
  protein['aatype'] = tf.argmax(
      protein['aatype'], axis=-1, output_type=tf.int32)
  for k in [
      'domain_name', 'msa', 'num_alignments', 'seq_length', 'sequence',
      'superfamily', 'deletion_matrix', 'resolution',
      'between_segment_residues', 'residue_index', 'template_all_atom_masks']:
    if k in protein:
      final_dim = shape_helpers.shape_list(protein[k])[-1]
      if isinstance(final_dim, int) and final_dim == 1:
        protein[k] = tf.squeeze(protein[k], axis=-1)

  for k in ['seq_length', 'num_alignments']:
    if k in protein:
      protein[k] = protein[k][0]  # Remove fake sequence dimension
  return protein


def make_random_crop_to_size_seed(protein):
  """Random seed for cropping residues and templates."""
  protein['random_crop_to_size_seed'] = utils.make_random_seed()
  return protein


@curry1
def randomly_replace_msa_with_unknown(protein, replace_proportion):
  """Replace a proportion of the MSA with 'X'."""
  msa_mask = (tf.random.uniform(shape_helpers.shape_list(protein['msa'])) <
              replace_proportion)
  x_idx = 20
  gap_idx = 21
  msa_mask = tf.logical_and(msa_mask, protein['msa'] != gap_idx)
  protein['msa'] = tf.where(msa_mask,
                            tf.ones_like(protein['msa']) * x_idx,
                            protein['msa'])
  aatype_mask = (
      tf.random.uniform(shape_helpers.shape_list(protein['aatype'])) <
      replace_proportion)

  protein['aatype'] = tf.where(aatype_mask,
                               tf.ones_like(protein['aatype']) * x_idx,
                               protein['aatype'])
  return protein


@curry1
def sample_msa(protein, max_seq, keep_extra):
  """Sample MSA randomly, remaining sequences are stored as `extra_*`.

  Args:
    protein: batch to sample msa from.
    max_seq: number of sequences to sample.
    keep_extra: When True sequences not sampled are put into fields starting
      with 'extra_*'.

  Returns:
    Protein with sampled msa.
  """
  num_seq = tf.shape(protein['msa'])[0]
  shuffled = tf.random_shuffle(tf.range(1, num_seq))
  index_order = tf.concat([[0], shuffled], axis=0)
  num_sel = tf.minimum(max_seq, num_seq)

  sel_seq, not_sel_seq = tf.split(index_order, [num_sel, num_seq - num_sel])

  for k in _MSA_FEATURE_NAMES:
    if k in protein:
      if keep_extra:
        protein['extra_' + k] = tf.gather(protein[k], not_sel_seq)
      protein[k] = tf.gather(protein[k], sel_seq)

  return protein


@curry1
def crop_extra_msa(protein, max_extra_msa):
  """MSA features are cropped so only `max_extra_msa` sequences are kept."""
  num_seq = tf.shape(protein['extra_msa'])[0]
  num_sel = tf.minimum(max_extra_msa, num_seq)
  select_indices = tf.random_shuffle(tf.range(0, num_seq))[:num_sel]
  for k in _MSA_FEATURE_NAMES:
    if 'extra_' + k in protein:
      protein['extra_' + k] = tf.gather(protein['extra_' + k], select_indices)

  return protein


def delete_extra_msa(protein):
  for k in _MSA_FEATURE_NAMES:
    if 'extra_' + k in protein:
      del protein['extra_' + k]
  return protein


@curry1
def block_delete_msa(protein, config):
  """Sample MSA by deleting contiguous blocks.

  Jumper et al. (2021) Suppl. Alg. 1 "MSABlockDeletion"

  Arguments:
    protein: batch dict containing the msa
    config: ConfigDict with parameters

  Returns:
    updated protein
  """
  num_seq = shape_helpers.shape_list(protein['msa'])[0]
  block_num_seq = tf.cast(
      tf.floor(tf.cast(num_seq, tf.float32) * config.msa_fraction_per_block),
      tf.int32)

  if config.randomize_num_blocks:
    nb = tf.random.uniform([], 0, config.num_blocks + 1, dtype=tf.int32)
  else:
    nb = config.num_blocks

  del_block_starts = tf.random.uniform([nb], 0, num_seq, dtype=tf.int32)
  del_blocks = del_block_starts[:, None] + tf.range(block_num_seq)
  del_blocks = tf.clip_by_value(del_blocks, 0, num_seq - 1)
  del_indices = tf.unique(tf.sort(tf.reshape(del_blocks, [-1])))[0]

  # Make sure we keep the original sequence
  sparse_diff = tf.sets.difference(tf.range(1, num_seq)[None],
                                   del_indices[None])
  keep_indices = tf.squeeze(tf.sparse.to_dense(sparse_diff), 0)
  keep_indices = tf.concat([[0], keep_indices], axis=0)

  for k in _MSA_FEATURE_NAMES:
    if k in protein:
      protein[k] = tf.gather(protein[k], keep_indices)

  return protein


@curry1
def nearest_neighbor_clusters(protein, gap_agreement_weight=0.):
  """Assign each extra MSA sequence to its nearest neighbor in sampled MSA."""

  # Determine how much weight we assign to each agreement.  In theory, we could
  # use a full blosum matrix here, but right now let's just down-weight gap
  # agreement because it could be spurious.
  # Never put weight on agreeing on BERT mask
  weights = tf.concat([
      tf.ones(21),
      gap_agreement_weight * tf.ones(1),
      np.zeros(1)], 0)

  # Make agreement score as weighted Hamming distance
  sample_one_hot = (protein['msa_mask'][:, :, None] *
                    tf.one_hot(protein['msa'], 23))
  extra_one_hot = (protein['extra_msa_mask'][:, :, None] *
                   tf.one_hot(protein['extra_msa'], 23))

  num_seq, num_res, _ = shape_helpers.shape_list(sample_one_hot)
  extra_num_seq, _, _ = shape_helpers.shape_list(extra_one_hot)

  # Compute tf.einsum('mrc,nrc,c->mn', sample_one_hot, extra_one_hot, weights)
  # in an optimized fashion to avoid possible memory or computation blowup.
  agreement = tf.matmul(
      tf.reshape(extra_one_hot, [extra_num_seq, num_res * 23]),
      tf.reshape(sample_one_hot * weights, [num_seq, num_res * 23]),
      transpose_b=True)

  # Assign each sequence in the extra sequences to the closest MSA sample
  protein['extra_cluster_assignment'] = tf.argmax(
      agreement, axis=1, output_type=tf.int32)

  return protein


@curry1
def summarize_clusters(protein):
  """Produce profile and deletion_matrix_mean within each cluster."""
  num_seq = shape_helpers.shape_list(protein['msa'])[0]
  def csum(x):
    return tf.math.unsorted_segment_sum(
        x, protein['extra_cluster_assignment'], num_seq)

  mask = protein['extra_msa_mask']
  mask_counts = 1e-6 + protein['msa_mask'] + csum(mask)  # Include center

  msa_sum = csum(mask[:, :, None] * tf.one_hot(protein['extra_msa'], 23))
  msa_sum += tf.one_hot(protein['msa'], 23)  # Original sequence
  protein['cluster_profile'] = msa_sum / mask_counts[:, :, None]

  del msa_sum

  del_sum = csum(mask * protein['extra_deletion_matrix'])
  del_sum += protein['deletion_matrix']  # Original sequence
  protein['cluster_deletion_mean'] = del_sum / mask_counts
  del del_sum

  return protein


def make_msa_mask(protein):
  """Mask features are all ones, but will later be zero-padded."""
  protein['msa_mask'] = tf.ones(
      shape_helpers.shape_list(protein['msa']), dtype=tf.float32)
  protein['msa_row_mask'] = tf.ones(
      shape_helpers.shape_list(protein['msa'])[0], dtype=tf.float32)
  return protein


def pseudo_beta_fn(aatype, all_atom_positions, all_atom_masks):
  """Create pseudo beta features."""
  is_gly = tf.equal(aatype, residue_constants.restype_order['G'])
  ca_idx = residue_constants.atom_order['CA']
  cb_idx = residue_constants.atom_order['CB']
  pseudo_beta = tf.where(
      tf.tile(is_gly[..., None], [1] * len(is_gly.shape) + [3]),
      all_atom_positions[..., ca_idx, :],
      all_atom_positions[..., cb_idx, :])

  if all_atom_masks is not None:
    pseudo_beta_mask = tf.where(
        is_gly, all_atom_masks[..., ca_idx], all_atom_masks[..., cb_idx])
    pseudo_beta_mask = tf.cast(pseudo_beta_mask, tf.float32)
    return pseudo_beta, pseudo_beta_mask
  else:
    return pseudo_beta


@curry1
def make_pseudo_beta(protein, prefix=''):
  """Create pseudo-beta (alpha for glycine) position and mask."""
  assert prefix in ['', 'template_']
  protein[prefix + 'pseudo_beta'], protein[prefix + 'pseudo_beta_mask'] = (
      pseudo_beta_fn(
          protein['template_aatype' if prefix else 'all_atom_aatype'],
          protein[prefix + 'all_atom_positions'],
          protein['template_all_atom_masks' if prefix else 'all_atom_mask']))
  return protein


@curry1
def add_constant_field(protein, key, value):
  protein[key] = tf.convert_to_tensor(value)
  return protein


def shaped_categorical(probs, epsilon=1e-10):
  ds = shape_helpers.shape_list(probs)
  num_classes = ds[-1]
  counts = tf.random.categorical(
      tf.reshape(tf.log(probs + epsilon), [-1, num_classes]),
      1,
      dtype=tf.int32)
  return tf.reshape(counts, ds[:-1])


def make_hhblits_profile(protein):
  """Compute the HHblits MSA profile if not already present."""
  if 'hhblits_profile' in protein:
    return protein

  # Compute the profile for every residue (over all MSA sequences).
  protein['hhblits_profile'] = tf.reduce_mean(
      tf.one_hot(protein['msa'], 22), axis=0)
  return protein


@curry1
def make_masked_msa(protein, config, replace_fraction):
  """Create data for BERT on raw MSA."""
  # Add a random amino acid uniformly
  random_aa = tf.constant([0.05] * 20 + [0., 0.], dtype=tf.float32)

  categorical_probs = (
      config.uniform_prob * random_aa +
      config.profile_prob * protein['hhblits_profile'] +
      config.same_prob * tf.one_hot(protein['msa'], 22))

  # Put all remaining probability on [MASK] which is a new column
  pad_shapes = [[0, 0] for _ in range(len(categorical_probs.shape))]
  pad_shapes[-1][1] = 1
  mask_prob = 1. - config.profile_prob - config.same_prob - config.uniform_prob
  assert mask_prob >= 0.
  categorical_probs = tf.pad(
      categorical_probs, pad_shapes, constant_values=mask_prob)

  sh = shape_helpers.shape_list(protein['msa'])
  mask_position = tf.random.uniform(sh) < replace_fraction

  bert_msa = shaped_categorical(categorical_probs)
  bert_msa = tf.where(mask_position, bert_msa, protein['msa'])

  # Mix real and masked MSA
  protein['bert_mask'] = tf.cast(mask_position, tf.float32)
  protein['true_msa'] = protein['msa']
  protein['msa'] = bert_msa

  return protein


@curry1
def make_fixed_size(protein, shape_schema, msa_cluster_size, extra_msa_size,
                    num_res, num_templates=0):
  """Guess at the MSA and sequence dimensions to make fixed size."""

  pad_size_map = {
      NUM_RES: num_res,
      NUM_MSA_SEQ: msa_cluster_size,
      NUM_EXTRA_SEQ: extra_msa_size,
      NUM_TEMPLATES: num_templates,
  }

  for k, v in protein.items():
    # Don't transfer this to the accelerator.
    if k == 'extra_cluster_assignment':
      continue
    shape = v.shape.as_list()
    schema = shape_schema[k]
    assert len(shape) == len(schema), (
        f'Rank mismatch between shape and shape schema for {k}: '
        f'{shape} vs {schema}')
    pad_size = [
        pad_size_map.get(s2, None) or s1 for (s1, s2) in zip(shape, schema)
    ]
    padding = [(0, p - tf.shape(v)[i]) for i, p in enumerate(pad_size)]
    if padding:
      protein[k] = tf.pad(
          v, padding, name=f'pad_to_fixed_{k}')
      protein[k].set_shape(pad_size)

  return protein


@curry1
def make_msa_feat(protein):
  """Create and concatenate MSA features."""
  # Whether there is a domain break. Always zero for chains, but keeping
  # for compatibility with domain datasets.
  has_break = tf.clip_by_value(
      tf.cast(protein['between_segment_residues'], tf.float32),
      0, 1)
  aatype_1hot = tf.one_hot(protein['aatype'], 21, axis=-1)

  target_feat = [
      tf.expand_dims(has_break, axis=-1),
      aatype_1hot,  # Everyone gets the original sequence.
  ]

  msa_1hot = tf.one_hot(protein['msa'], 23, axis=-1)
  has_deletion = tf.clip_by_value(protein['deletion_matrix'], 0., 1.)
  deletion_value = tf.atan(protein['deletion_matrix'] / 3.) * (2. / np.pi)

  msa_feat = [
      msa_1hot,
      tf.expand_dims(has_deletion, axis=-1),
      tf.expand_dims(deletion_value, axis=-1),
  ]

  if 'cluster_profile' in protein:
    deletion_mean_value = (
        tf.atan(protein['cluster_deletion_mean'] / 3.) * (2. / np.pi))
    msa_feat.extend([
        protein['cluster_profile'],
        tf.expand_dims(deletion_mean_value, axis=-1),
    ])

  if 'extra_deletion_matrix' in protein:
    protein['extra_has_deletion'] = tf.clip_by_value(
        protein['extra_deletion_matrix'], 0., 1.)
    protein['extra_deletion_value'] = tf.atan(
        protein['extra_deletion_matrix'] / 3.) * (2. / np.pi)

  protein['msa_feat'] = tf.concat(msa_feat, axis=-1)
  protein['target_feat'] = tf.concat(target_feat, axis=-1)
  return protein


@curry1
def select_feat(protein, feature_list):
  return {k: v for k, v in protein.items() if k in feature_list}


@curry1
def crop_templates(protein, max_templates):
  for k, v in protein.items():
    if k.startswith('template_'):
      protein[k] = v[:max_templates]
  return protein


@curry1
def random_crop_to_size(protein, crop_size, max_templates, shape_schema,
                        subsample_templates=False):
  """Crop randomly to `crop_size`, or keep as is if shorter than that."""
  seq_length = protein['seq_length']
  if 'template_mask' in protein:
    num_templates = tf.cast(
        shape_helpers.shape_list(protein['template_mask'])[0], tf.int32)
  else:
    num_templates = tf.constant(0, dtype=tf.int32)
  num_res_crop_size = tf.math.minimum(seq_length, crop_size)

  # Ensures that the cropping of residues and templates happens in the same way
  # across ensembling iterations.
  # Do not use for randomness that should vary in ensembling.
  seed_maker = utils.SeedMaker(initial_seed=protein['random_crop_to_size_seed'])

  if subsample_templates:
    templates_crop_start = tf.random.stateless_uniform(
        shape=(), minval=0, maxval=num_templates + 1, dtype=tf.int32,
        seed=seed_maker())
  else:
    templates_crop_start = 0

  num_templates_crop_size = tf.math.minimum(
      num_templates - templates_crop_start, max_templates)

  num_res_crop_start = tf.random.stateless_uniform(
      shape=(), minval=0, maxval=seq_length - num_res_crop_size + 1,
      dtype=tf.int32, seed=seed_maker())

  templates_select_indices = tf.argsort(tf.random.stateless_uniform(
      [num_templates], seed=seed_maker()))

  for k, v in protein.items():
    if k not in shape_schema or (
        'template' not in k and NUM_RES not in shape_schema[k]):
      continue

    # randomly permute the templates before cropping them.
    if k.startswith('template') and subsample_templates:
      v = tf.gather(v, templates_select_indices)

    crop_sizes = []
    crop_starts = []
    for i, (dim_size, dim) in enumerate(zip(shape_schema[k],
                                            shape_helpers.shape_list(v))):
      is_num_res = (dim_size == NUM_RES)
      if i == 0 and k.startswith('template'):
        crop_size = num_templates_crop_size
        crop_start = templates_crop_start
      else:
        crop_start = num_res_crop_start if is_num_res else 0
        crop_size = (num_res_crop_size if is_num_res else
                     (-1 if dim is None else dim))
      crop_sizes.append(crop_size)
      crop_starts.append(crop_start)
    protein[k] = tf.slice(v, crop_starts, crop_sizes)

  protein['seq_length'] = num_res_crop_size
  return protein


def make_atom14_masks(protein):
  """Construct denser atom positions (14 dimensions instead of 37)."""
  restype_atom14_to_atom37 = []  # mapping (restype, atom14) --> atom37
  restype_atom37_to_atom14 = []  # mapping (restype, atom37) --> atom14
  restype_atom14_mask = []

  for rt in residue_constants.restypes:
    atom_names = residue_constants.restype_name_to_atom14_names[
        residue_constants.restype_1to3[rt]]

    restype_atom14_to_atom37.append([
        (residue_constants.atom_order[name] if name else 0)
        for name in atom_names
    ])

    atom_name_to_idx14 = {name: i for i, name in enumerate(atom_names)}
    restype_atom37_to_atom14.append([
        (atom_name_to_idx14[name] if name in atom_name_to_idx14 else 0)
        for name in residue_constants.atom_types
    ])

    restype_atom14_mask.append([(1. if name else 0.) for name in atom_names])

  # Add dummy mapping for restype 'UNK'
  restype_atom14_to_atom37.append([0] * 14)
  restype_atom37_to_atom14.append([0] * 37)
  restype_atom14_mask.append([0.] * 14)

  restype_atom14_to_atom37 = np.array(restype_atom14_to_atom37, dtype=np.int32)
  restype_atom37_to_atom14 = np.array(restype_atom37_to_atom14, dtype=np.int32)
  restype_atom14_mask = np.array(restype_atom14_mask, dtype=np.float32)

  # create the mapping for (residx, atom14) --> atom37, i.e. an array
  # with shape (num_res, 14) containing the atom37 indices for this protein
  residx_atom14_to_atom37 = tf.gather(restype_atom14_to_atom37,
                                      protein['aatype'])
  residx_atom14_mask = tf.gather(restype_atom14_mask,
                                 protein['aatype'])

  protein['atom14_atom_exists'] = residx_atom14_mask
  protein['residx_atom14_to_atom37'] = residx_atom14_to_atom37

  # create the gather indices for mapping back
  residx_atom37_to_atom14 = tf.gather(restype_atom37_to_atom14,
                                      protein['aatype'])
  protein['residx_atom37_to_atom14'] = residx_atom37_to_atom14

  # create the corresponding mask
  restype_atom37_mask = np.zeros([21, 37], dtype=np.float32)
  for restype, restype_letter in enumerate(residue_constants.restypes):
    restype_name = residue_constants.restype_1to3[restype_letter]
    atom_names = residue_constants.residue_atoms[restype_name]
    for atom_name in atom_names:
      atom_type = residue_constants.atom_order[atom_name]
      restype_atom37_mask[restype, atom_type] = 1

  residx_atom37_mask = tf.gather(restype_atom37_mask,
                                 protein['aatype'])
  protein['atom37_atom_exists'] = residx_atom37_mask

  return protein

