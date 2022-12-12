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

"""A collection of common Haiku modules for use in protein folding."""
import numbers
from typing import Union, Sequence

import haiku as hk
import jax.numpy as jnp
import numpy as np


# Constant from scipy.stats.truncnorm.std(a=-2, b=2, loc=0., scale=1.)
TRUNCATED_NORMAL_STDDEV_FACTOR = np.asarray(.87962566103423978,
                                            dtype=np.float32)


def get_initializer_scale(initializer_name, input_shape):
  """Get Initializer for weights and scale to multiply activations by."""

  if initializer_name == 'zeros':
    w_init = hk.initializers.Constant(0.0)
  else:
    # fan-in scaling
    scale = 1.
    for channel_dim in input_shape:
      scale /= channel_dim
    if initializer_name == 'relu':
      scale *= 2

    noise_scale = scale

    stddev = np.sqrt(noise_scale)
    # Adjust stddev for truncation.
    stddev = stddev / TRUNCATED_NORMAL_STDDEV_FACTOR
    w_init = hk.initializers.TruncatedNormal(mean=0.0, stddev=stddev)

  return w_init


class Linear(hk.Module):
  """Protein folding specific Linear module.

  This differs from the standard Haiku Linear in a few ways:
    * It supports inputs and outputs of arbitrary rank
    * Initializers are specified by strings
  """

  def __init__(self,
               num_output: Union[int, Sequence[int]],
               initializer: str = 'linear',
               num_input_dims: int = 1,
               use_bias: bool = True,
               bias_init: float = 0.,
               precision = None,
               name: str = 'linear'):
    """Constructs Linear Module.

    Args:
      num_output: Number of output channels. Can be tuple when outputting
          multiple dimensions.
      initializer: What initializer to use, should be one of {'linear', 'relu',
        'zeros'}
      num_input_dims: Number of dimensions from the end to project.
      use_bias: Whether to include trainable bias
      bias_init: Value used to initialize bias.
      precision: What precision to use for matrix multiplication, defaults
        to None.
      name: Name of module, used for name scopes.
    """
    super().__init__(name=name)
    if isinstance(num_output, numbers.Integral):
      self.output_shape = (num_output,)
    else:
      self.output_shape = tuple(num_output)
    self.initializer = initializer
    self.use_bias = use_bias
    self.bias_init = bias_init
    self.num_input_dims = num_input_dims
    self.num_output_dims = len(self.output_shape)
    self.precision = precision

  def __call__(self, inputs):
    """Connects Module.

    Args:
      inputs: Tensor with at least num_input_dims dimensions.

    Returns:
      output of shape [...] + num_output.
    """

    num_input_dims = self.num_input_dims

    if self.num_input_dims > 0:
      in_shape = inputs.shape[-self.num_input_dims:]
    else:
      in_shape = ()

    weight_init = get_initializer_scale(self.initializer, in_shape)

    in_letters = 'abcde'[:self.num_input_dims]
    out_letters = 'hijkl'[:self.num_output_dims]

    weight_shape = in_shape + self.output_shape
    weights = hk.get_parameter('weights', weight_shape, inputs.dtype,
                               weight_init)

    equation = f'...{in_letters}, {in_letters}{out_letters}->...{out_letters}'

    output = jnp.einsum(equation, inputs, weights, precision=self.precision)

    if self.use_bias:
      bias = hk.get_parameter('bias', self.output_shape, inputs.dtype,
                              hk.initializers.Constant(self.bias_init))
      output += bias

    return output


class LayerNorm(hk.LayerNorm):
  """LayerNorm module.

  Equivalent to hk.LayerNorm but with different parameter shapes: they are
  always vectors rather than possibly higher-rank tensors. This makes it easier
  to change the layout whilst keep the model weight-compatible.
  """

  def __init__(self,
               axis,
               create_scale: bool,
               create_offset: bool,
               eps: float = 1e-5,
               scale_init=None,
               offset_init=None,
               use_fast_variance: bool = False,
               name=None,
               param_axis=None):
    super().__init__(
        axis=axis,
        create_scale=False,
        create_offset=False,
        eps=eps,
        scale_init=None,
        offset_init=None,
        use_fast_variance=use_fast_variance,
        name=name,
        param_axis=param_axis)
    self._temp_create_scale = create_scale
    self._temp_create_offset = create_offset

  def __call__(self, x: jnp.ndarray) -> jnp.ndarray:
    is_bf16 = (x.dtype == jnp.bfloat16)
    if is_bf16:
      x = x.astype(jnp.float32)

    param_axis = self.param_axis[0] if self.param_axis else -1
    param_shape = (x.shape[param_axis],)

    param_broadcast_shape = [1] * x.ndim
    param_broadcast_shape[param_axis] = x.shape[param_axis]
    scale = None
    offset = None
    if self._temp_create_scale:
      scale = hk.get_parameter(
          'scale', param_shape, x.dtype, init=self.scale_init)
      scale = scale.reshape(param_broadcast_shape)

    if self._temp_create_offset:
      offset = hk.get_parameter(
          'offset', param_shape, x.dtype, init=self.offset_init)
      offset = offset.reshape(param_broadcast_shape)

    out = super().__call__(x, scale=scale, offset=offset)

    if is_bf16:
      out = out.astype(jnp.bfloat16)

    return out
  