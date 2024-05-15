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
"""Rigid3Array Transformations represented by a Matrix and a Vector."""

from __future__ import annotations
from typing import Union

from alphafold.model.geometry import rotation_matrix
from alphafold.model.geometry import struct_of_array
from alphafold.model.geometry import vector
import jax
import jax.numpy as jnp

Float = Union[float, jnp.ndarray]

VERSION = '0.1'


@struct_of_array.StructOfArray(same_dtype=True)
class Rigid3Array:
  """Rigid Transformation, i.e. element of special euclidean group."""

  rotation: rotation_matrix.Rot3Array
  translation: vector.Vec3Array

  def __matmul__(self, other: Rigid3Array) -> Rigid3Array:
    new_rotation = self.rotation @ other.rotation
    new_translation = self.apply_to_point(other.translation)
    return Rigid3Array(new_rotation, new_translation)

  def inverse(self) -> Rigid3Array:
    """Return Rigid3Array corresponding to inverse transform."""
    inv_rotation = self.rotation.inverse()
    inv_translation = inv_rotation.apply_to_point(-self.translation)
    return Rigid3Array(inv_rotation, inv_translation)

  def apply_to_point(self, point: vector.Vec3Array) -> vector.Vec3Array:
    """Apply Rigid3Array transform to point."""
    return self.rotation.apply_to_point(point) + self.translation

  def apply_inverse_to_point(self, point: vector.Vec3Array) -> vector.Vec3Array:
    """Apply inverse Rigid3Array transform to point."""
    new_point = point - self.translation
    return self.rotation.apply_inverse_to_point(new_point)

  def compose_rotation(self, other_rotation):
    rot = self.rotation @ other_rotation
    trans = jax.tree.map(lambda x: jnp.broadcast_to(x, rot.shape),
                         self.translation)
    return Rigid3Array(rot, trans)

  @classmethod
  def identity(cls, shape, dtype=jnp.float32) -> Rigid3Array:
    """Return identity Rigid3Array of given shape."""
    return cls(
        rotation_matrix.Rot3Array.identity(shape, dtype=dtype),
        vector.Vec3Array.zeros(shape, dtype=dtype))  # pytype: disable=wrong-arg-count  # trace-all-classes

  def scale_translation(self, factor: Float) -> Rigid3Array:
    """Scale translation in Rigid3Array by 'factor'."""
    return Rigid3Array(self.rotation, self.translation * factor)

  def to_array(self):
    rot_array = self.rotation.to_array()
    vec_array = self.translation.to_array()
    return jnp.concatenate([rot_array, vec_array[..., None]], axis=-1)

  @classmethod
  def from_array(cls, array):
    rot = rotation_matrix.Rot3Array.from_array(array[..., :3])
    vec = vector.Vec3Array.from_array(array[..., -1])
    return cls(rot, vec)  # pytype: disable=wrong-arg-count  # trace-all-classes

  @classmethod
  def from_array4x4(cls, array: jnp.ndarray) -> Rigid3Array:
    """Construct Rigid3Array from homogeneous 4x4 array."""
    assert array.shape[-1] == 4
    assert array.shape[-2] == 4
    rotation = rotation_matrix.Rot3Array(
        array[..., 0, 0], array[..., 0, 1], array[..., 0, 2],
        array[..., 1, 0], array[..., 1, 1], array[..., 1, 2],
        array[..., 2, 0], array[..., 2, 1], array[..., 2, 2]
        )
    translation = vector.Vec3Array(
        array[..., 0, 3], array[..., 1, 3], array[..., 2, 3])
    return cls(rotation, translation)  # pytype: disable=wrong-arg-count  # trace-all-classes

  def __getstate__(self):
    return (VERSION, (self.rotation, self.translation))

  def __setstate__(self, state):
    version, (rot, trans) = state
    del version
    object.__setattr__(self, 'rotation', rot)
    object.__setattr__(self, 'translation', trans)
