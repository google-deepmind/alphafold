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

# Advanced Haiku Protein Folding Modules
import haiku as hk
import jax.numpy as jnp
import numpy as np

class ProteinLinear(hk.Module):
    def __init__(self, num_output, initializer='linear', num_input_dims=1, use_bias=True, bias_init=0., name='protein_linear'):
        super().__init__(name=name)
        self.num_output = num_output
        self.initializer = initializer
        self.num_input_dims = num_input_dims
        self.use_bias = use_bias
        self.bias_init = bias_init

    def __call__(self, inputs):
        input_shape = inputs.shape[-self.num_input_dims:]
        weight_init = self._get_initializer_scale(input_shape)
        weights = hk.get_parameter('weights', input_shape + self.num_output, inputs.dtype, weight_init)
        output = jnp.matmul(inputs, weights)

        if self.use_bias:
            bias = hk.get_parameter('bias', self.num_output, inputs.dtype, hk.initializers.Constant(self.bias_init))
            output += bias

        return output

    def _get_initializer_scale(self, input_shape):
        if self.initializer == 'zeros':
            return hk.initializers.Constant(0.0)
        else:
            scale = 1.0 / jnp.prod(input_shape)
            if self.initializer == 'relu':
                scale *= 2
            stddev = jnp.sqrt(scale)
            stddev /= np.sqrt(0.87962566103423978)  # Adjusted for truncation
            return hk.initializers.TruncatedNormal(mean=0.0, stddev=stddev)

class ProteinLayerNorm(hk.LayerNorm):
    def __init__(self, axis, create_scale=True, create_offset=True, eps=1e-5, name=None):
        super().__init__(axis=axis, create_scale=False, create_offset=False, eps=eps, name=name)
        self._temp_create_scale = create_scale
        self._temp_create_offset = create_offset

    def __call__(self, x):
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
            scale = hk.get_parameter('scale', param_shape, x.dtype, init=self.scale_init)
            scale = scale.reshape(param_broadcast_shape)

        if self._temp_create_offset:
            offset = hk.get_parameter('offset', param_shape, x.dtype, init=self.offset_init)
            offset = offset.reshape(param_broadcast_shape)

        out = super().__call__(x, scale=scale, offset=offset)

        if is_bf16:
            out = out.astype(jnp.bfloat16)

        return out
