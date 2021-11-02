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

"""A collection of JAX utility functions for use in protein folding."""

import collections
import functools
import numbers
from typing import Mapping

import haiku as hk
import jax
import jax.numpy as jnp
import numpy as np


def final_init(config):
  if config.zero_init:
    return 'zeros'
  else:
    return 'linear'


def batched_gather(params, indices, axis=0, batch_dims=0):
  """Implements a JAX equivalent of `tf.gather` with `axis` and `batch_dims`."""
  take_fn = lambda p, i: jnp.take(p, i, axis=axis)
  for _ in range(batch_dims):
    take_fn = jax.vmap(take_fn)
  return take_fn(params, indices)


def mask_mean(mask, value, axis=None, drop_mask_channel=False, eps=1e-10):
  """Masked mean."""
  if drop_mask_channel:
    mask = mask[..., 0]

  mask_shape = mask.shape
  value_shape = value.shape

  assert len(mask_shape) == len(value_shape)

  if isinstance(axis, numbers.Integral):
    axis = [axis]
  elif axis is None:
    axis = list(range(len(mask_shape)))
  assert isinstance(axis, collections.Iterable), (
      'axis needs to be either an iterable, integer or "None"')

  broadcast_factor = 1.
  for axis_ in axis:
    value_size = value_shape[axis_]
    mask_size = mask_shape[axis_]
    if mask_size == 1:
      broadcast_factor *= value_size
    else:
      assert mask_size == value_size

  return (jnp.sum(mask * value, axis=axis) /
          (jnp.sum(mask, axis=axis) * broadcast_factor + eps))


def flat_params_to_haiku(params: Mapping[str, np.ndarray]) -> hk.Params:
  """Convert a dictionary of NumPy arrays to Haiku parameters."""
  hk_params = {}
  for path, array in params.items():
    scope, name = path.split('//')
    if scope not in hk_params:
      hk_params[scope] = {}
    hk_params[scope][name] = jnp.array(array)

  return hk_params


def padding_consistent_rng(f):
  """Modify any element-wise random function to be consistent with padding.

  Normally if you take a function like jax.random.normal and generate an array,
  say of size (10,10), you will get a different set of random numbers to if you
  add padding and take the first (10,10) sub-array.

  This function makes a random function that is consistent regardless of the
  amount of padding added.

  Note: The padding-consistent function is likely to be slower to compile and
  run than the function it is wrapping, but these slowdowns are likely to be
  negligible in a large network.

  Args:
    f: Any element-wise function that takes (PRNG key, shape) as the first 2
      arguments.

  Returns:
    An equivalent function to f, that is now consistent for different amounts of
    padding.
  """
  def grid_keys(key, shape):
    """Generate a grid of rng keys that is consistent with different padding.

    Generate random keys such that the keys will be identical, regardless of
    how much padding is added to any dimension.

    Args:
      key: A PRNG key.
      shape: The shape of the output array of keys that will be generated.

    Returns:
      An array of shape `shape` consisting of random keys.
    """
    if not shape:
      return key
    new_keys = jax.vmap(functools.partial(jax.random.fold_in, key))(
        jnp.arange(shape[0]))
    return jax.vmap(functools.partial(grid_keys, shape=shape[1:]))(new_keys)

  def inner(key, shape, **kwargs):
    return jnp.vectorize(
        lambda key: f(key, shape=(), **kwargs),
        signature='(2)->()')(
            grid_keys(key, shape))
  return inner
