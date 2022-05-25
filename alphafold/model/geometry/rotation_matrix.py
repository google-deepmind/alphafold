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
"""Rot3Array Matrix Class."""

from __future__ import annotations
import dataclasses

from alphafold.model.geometry import struct_of_array
from alphafold.model.geometry import utils
from alphafold.model.geometry import vector
import jax
import jax.numpy as jnp
import numpy as np

COMPONENTS = ['xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz']

VERSION = '0.1'


@struct_of_array.StructOfArray(same_dtype=True)
class Rot3Array:
  """Rot3Array Matrix in 3 dimensional Space implemented as struct of arrays."""

  xx: jnp.ndarray = dataclasses.field(metadata={'dtype': jnp.float32})
  xy: jnp.ndarray
  xz: jnp.ndarray
  yx: jnp.ndarray
  yy: jnp.ndarray
  yz: jnp.ndarray
  zx: jnp.ndarray
  zy: jnp.ndarray
  zz: jnp.ndarray

  __array_ufunc__ = None

  def inverse(self) -> Rot3Array:
    """Returns inverse of Rot3Array."""
    return Rot3Array(self.xx, self.yx, self.zx,
                     self.xy, self.yy, self.zy,
                     self.xz, self.yz, self.zz)

  def apply_to_point(self, point: vector.Vec3Array) -> vector.Vec3Array:
    """Applies Rot3Array to point."""
    return vector.Vec3Array(
        self.xx * point.x + self.xy * point.y + self.xz * point.z,
        self.yx * point.x + self.yy * point.y + self.yz * point.z,
        self.zx * point.x + self.zy * point.y + self.zz * point.z)

  def apply_inverse_to_point(self, point: vector.Vec3Array) -> vector.Vec3Array:
    """Applies inverse Rot3Array to point."""
    return self.inverse().apply_to_point(point)

  def __matmul__(self, other: Rot3Array) -> Rot3Array:
    """Composes two Rot3Arrays."""
    c0 = self.apply_to_point(vector.Vec3Array(other.xx, other.yx, other.zx))
    c1 = self.apply_to_point(vector.Vec3Array(other.xy, other.yy, other.zy))
    c2 = self.apply_to_point(vector.Vec3Array(other.xz, other.yz, other.zz))
    return Rot3Array(c0.x, c1.x, c2.x, c0.y, c1.y, c2.y, c0.z, c1.z, c2.z)

  @classmethod
  def identity(cls, shape, dtype=jnp.float32) -> Rot3Array:
    """Returns identity of given shape."""
    ones = jnp.ones(shape, dtype=dtype)
    zeros = jnp.zeros(shape, dtype=dtype)
    return cls(ones, zeros, zeros, zeros, ones, zeros, zeros, zeros, ones)  # pytype: disable=wrong-arg-count  # trace-all-classes

  @classmethod
  def from_two_vectors(cls, e0: vector.Vec3Array,
                       e1: vector.Vec3Array) -> Rot3Array:
    """Construct Rot3Array from two Vectors.

    Rot3Array is constructed such that in the corresponding frame 'e0' lies on
    the positive x-Axis and 'e1' lies in the xy plane with positive sign of y.

    Args:
      e0: Vector
      e1: Vector
    Returns:
      Rot3Array
    """
    # Normalize the unit vector for the x-axis, e0.
    e0 = e0.normalized()
    # make e1 perpendicular to e0.
    c = e1.dot(e0)
    e1 = (e1 - c * e0).normalized()
    # Compute e2 as cross product of e0 and e1.
    e2 = e0.cross(e1)
    return cls(e0.x, e1.x, e2.x, e0.y, e1.y, e2.y, e0.z, e1.z, e2.z)  # pytype: disable=wrong-arg-count  # trace-all-classes

  @classmethod
  def from_array(cls, array: jnp.ndarray) -> Rot3Array:
    """Construct Rot3Array Matrix from array of shape. [..., 3, 3]."""
    unstacked = utils.unstack(array, axis=-2)
    unstacked = sum([utils.unstack(x, axis=-1) for x in unstacked], [])
    return cls(*unstacked)

  def to_array(self) -> jnp.ndarray:
    """Convert Rot3Array to array of shape [..., 3, 3]."""
    return jnp.stack(
        [jnp.stack([self.xx, self.xy, self.xz], axis=-1),
         jnp.stack([self.yx, self.yy, self.yz], axis=-1),
         jnp.stack([self.zx, self.zy, self.zz], axis=-1)],
        axis=-2)

  @classmethod
  def from_quaternion(cls,
                      w: jnp.ndarray,
                      x: jnp.ndarray,
                      y: jnp.ndarray,
                      z: jnp.ndarray,
                      normalize: bool = True,
                      epsilon: float = 1e-6) -> Rot3Array:
    """Construct Rot3Array from components of quaternion."""
    if normalize:
      inv_norm = jax.lax.rsqrt(jnp.maximum(epsilon, w**2 + x**2 + y**2 + z**2))
      w *= inv_norm
      x *= inv_norm
      y *= inv_norm
      z *= inv_norm
    xx = 1 - 2 * (jnp.square(y) + jnp.square(z))
    xy = 2 * (x * y - w * z)
    xz = 2 * (x * z + w * y)
    yx = 2 * (x * y + w * z)
    yy = 1 - 2 * (jnp.square(x) + jnp.square(z))
    yz = 2 * (y * z - w * x)
    zx = 2 * (x * z - w * y)
    zy = 2 * (y * z + w * x)
    zz = 1 - 2 * (jnp.square(x) + jnp.square(y))
    return cls(xx, xy, xz, yx, yy, yz, zx, zy, zz)  # pytype: disable=wrong-arg-count  # trace-all-classes

  @classmethod
  def random_uniform(cls, key, shape, dtype=jnp.float32) -> Rot3Array:
    """Samples uniform random Rot3Array according to Haar Measure."""
    quat_array = jax.random.normal(key, tuple(shape) + (4,), dtype=dtype)
    quats = utils.unstack(quat_array)
    return cls.from_quaternion(*quats)

  def __getstate__(self):
    return (VERSION,
            [np.asarray(getattr(self, field)) for field in COMPONENTS])

  def __setstate__(self, state):
    version, state = state
    del version
    for i, field in enumerate(COMPONENTS):
      object.__setattr__(self, field, state[i])
