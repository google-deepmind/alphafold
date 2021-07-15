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
import haiku as hk
import jax.numpy as jnp


class Linear(hk.Module):
  """Protein folding specific Linear Module.

  This differs from the standard Haiku Linear in a few ways:
    * It supports inputs of arbitrary rank
    * Initializers are specified by strings
  """

  def __init__(self,
               num_output: int,
               initializer: str = 'linear',
               use_bias: bool = True,
               bias_init: float = 0.,
               name: str = 'linear'):
    """Constructs Linear Module.

    Args:
      num_output: number of output channels.
      initializer: What initializer to use, should be one of {'linear', 'relu',
        'zeros'}
      use_bias: Whether to include trainable bias
      bias_init: Value used to initialize bias.
      name: name of module, used for name scopes.
    """

    super().__init__(name=name)
    self.num_output = num_output
    self.initializer = initializer
    self.use_bias = use_bias
    self.bias_init = bias_init

  def __call__(self, inputs: jnp.ndarray) -> jnp.ndarray:
    """Connects Module.

    Args:
      inputs: Tensor of shape [..., num_channel]

    Returns:
      output of shape [..., num_output]
    """
    n_channels = int(inputs.shape[-1])

    weight_shape = [n_channels, self.num_output]
    if self.initializer == 'linear':
      weight_init = hk.initializers.VarianceScaling(mode='fan_in', scale=1.)
    elif self.initializer == 'relu':
      weight_init = hk.initializers.VarianceScaling(mode='fan_in', scale=2.)
    elif self.initializer == 'zeros':
      weight_init = hk.initializers.Constant(0.0)

    weights = hk.get_parameter('weights', weight_shape, inputs.dtype,
                               weight_init)

    # this is equivalent to einsum('...c,cd->...d', inputs, weights)
    # but turns out to be slightly faster
    inputs = jnp.swapaxes(inputs, -1, -2)
    output = jnp.einsum('...cb,cd->...db', inputs, weights)
    output = jnp.swapaxes(output, -1, -2)

    if self.use_bias:
      bias = hk.get_parameter('bias', [self.num_output], inputs.dtype,
                              hk.initializers.Constant(self.bias_init))
      output += bias

    return output
