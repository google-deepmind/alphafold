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

"""Utilities for dealing with shapes of TensorFlow tensors."""
import tensorflow.compat.v1 as tf


def shape_list(x):
  """Return list of dimensions of a tensor, statically where possible.

  Like `x.shape.as_list()` but with tensors instead of `None`s.

  Args:
    x: A tensor.
  Returns:
    A list with length equal to the rank of the tensor. The n-th element of the
    list is an integer when that dimension is statically known otherwise it is
    the n-th element of `tf.shape(x)`.
  """
  x = tf.convert_to_tensor(x)

  # If unknown rank, return dynamic shape
  if x.get_shape().dims is None:
    return tf.shape(x)

  static = x.get_shape().as_list()
  shape = tf.shape(x)

  ret = []
  for i in range(len(static)):
    dim = static[i]
    if dim is None:
      dim = shape[i]
    ret.append(dim)
  return ret

