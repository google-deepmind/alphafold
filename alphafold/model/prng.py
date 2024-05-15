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

"""A collection of utilities surrounding PRNG usage in protein folding."""

import haiku as hk
import jax


def safe_dropout(*, tensor, safe_key, rate, is_deterministic, is_training):
  if is_training and rate != 0.0 and not is_deterministic:
    return hk.dropout(safe_key.get(), rate, tensor)
  else:
    return tensor


class SafeKey:
  """Safety wrapper for PRNG keys."""

  def __init__(self, key):
    self._key = key
    self._used = False

  def _assert_not_used(self):
    if self._used:
      raise RuntimeError('Random key has been used previously.')

  def get(self):
    self._assert_not_used()
    self._used = True
    return self._key

  def split(self, num_keys=2):
    self._assert_not_used()
    self._used = True
    new_keys = jax.random.split(self._key, num_keys)
    return jax.tree.map(SafeKey, tuple(new_keys))

  def duplicate(self, num_keys=2):
    self._assert_not_used()
    self._used = True
    return tuple(SafeKey(self._key) for _ in range(num_keys))


def _safe_key_flatten(safe_key):
  # Flatten transfers "ownership" to the tree
  return (safe_key._key,), safe_key._used  # pylint: disable=protected-access


def _safe_key_unflatten(aux_data, children):
  ret = SafeKey(children[0])
  ret._used = aux_data  # pylint: disable=protected-access
  return ret


jax.tree_util.register_pytree_node(
    SafeKey, _safe_key_flatten, _safe_key_unflatten)

