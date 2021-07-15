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

"""Tests for layer_stack."""

import functools
from absl.testing import absltest
from absl.testing import parameterized
import haiku as hk
import jax
import jax.numpy as jnp
import numpy as np
import scipy

from alphafold.model import layer_stack


# Suffixes applied by Haiku for repeated module names.
suffixes = [''] + [f'_{i}' for i in range(1, 100)]


def _slice_layers_params(layers_params):
  sliced_layers_params = {}
  for k, v in layers_params.items():
    for inner_k in v:
      for var_slice, suffix in zip(v[inner_k], suffixes):
        k_new = k.split('/')[-1] + suffix
        if k_new not in sliced_layers_params:
          sliced_layers_params[k_new] = {}
        sliced_layers_params[k_new][inner_k] = var_slice
  return sliced_layers_params


class LayerStackTest(parameterized.TestCase):

  @parameterized.parameters([1, 2, 4])
  def test_layer_stack(self, unroll):
    """Compare layer_stack to the equivalent unrolled stack.

    Tests that the layer_stack application of a Haiku layer function is
    equivalent to repeatedly applying the layer function in an unrolled loop.

    Args:
      unroll: Number of unrolled layers.
    """
    num_layers = 20

    def inner_fn(x):
      x += hk.Linear(100, name='linear1')(x)
      x += hk.Linear(100, name='linear2')(x)
      return x

    def outer_fn_unrolled(x):
      for _ in range(num_layers):
        x = inner_fn(x)
      return x

    def outer_fn_layer_stack(x):
      stack = layer_stack.layer_stack(num_layers, unroll=unroll)(inner_fn)
      return stack(x)

    unrolled_fn = hk.transform(outer_fn_unrolled)
    layer_stack_fn = hk.transform(outer_fn_layer_stack)

    x = jax.random.uniform(jax.random.PRNGKey(0), [10, 256, 100])

    rng_init = jax.random.PRNGKey(42)

    params = layer_stack_fn.init(rng_init, x)

    sliced_params = _slice_layers_params(params)

    unrolled_pred = unrolled_fn.apply(sliced_params, None, x)
    layer_stack_pred = layer_stack_fn.apply(params, None, x)

    np.testing.assert_allclose(unrolled_pred, layer_stack_pred)

  def test_layer_stack_multi_args(self):
    """Compare layer_stack to the equivalent unrolled stack.

    Similar to `test_layer_stack`, but use a function that takes more than one
    argument.
    """
    num_layers = 20

    def inner_fn(x, y):
      x_out = x + hk.Linear(100, name='linear1')(y)
      y_out = y + hk.Linear(100, name='linear2')(x)
      return x_out, y_out

    def outer_fn_unrolled(x, y):
      for _ in range(num_layers):
        x, y = inner_fn(x, y)
      return x, y

    def outer_fn_layer_stack(x, y):
      stack = layer_stack.layer_stack(num_layers)(inner_fn)
      return stack(x, y)

    unrolled_fn = hk.transform(outer_fn_unrolled)
    layer_stack_fn = hk.transform(outer_fn_layer_stack)

    x = jax.random.uniform(jax.random.PRNGKey(0), [10, 256, 100])
    y = jax.random.uniform(jax.random.PRNGKey(1), [10, 256, 100])

    rng_init = jax.random.PRNGKey(42)

    params = layer_stack_fn.init(rng_init, x, y)

    sliced_params = _slice_layers_params(params)

    unrolled_x, unrolled_y = unrolled_fn.apply(sliced_params, None, x, y)
    layer_stack_x, layer_stack_y = layer_stack_fn.apply(params, None, x, y)

    np.testing.assert_allclose(unrolled_x, layer_stack_x)
    np.testing.assert_allclose(unrolled_y, layer_stack_y)

  def test_layer_stack_no_varargs(self):
    """Test an error is raised when using a function with varargs."""

    class VarArgsModule(hk.Module):
      """When used, this module should cause layer_stack to raise an Error."""

      def __call__(self, *args):
        return args

    class NoVarArgsModule(hk.Module):
      """This module should be fine to use with layer_stack."""

      def __call__(self, x):
        return x

    def build_and_init_stack(module_class):
      def stack_fn(x):
        module = module_class()
        return layer_stack.layer_stack(1)(module)(x)

      stack = hk.without_apply_rng(hk.transform(stack_fn))
      stack.init(jax.random.PRNGKey(1729), jnp.ones([5]))

    build_and_init_stack(NoVarArgsModule)
    with self.assertRaisesRegex(
        ValueError, 'The function `f` should not have any `varargs`'):
      build_and_init_stack(VarArgsModule)

  @parameterized.parameters([1, 2, 4])
  def test_layer_stack_grads(self, unroll):
    """Compare layer_stack gradients to the equivalent unrolled stack.

    Tests that the layer_stack application of a Haiku layer function is
    equivalent to repeatedly applying the layer function in an unrolled loop.

    Args:
      unroll: Number of unrolled layers.
    """
    num_layers = 20

    def inner_fn(x):
      x += hk.Linear(100, name='linear1')(x)
      x += hk.Linear(100, name='linear2')(x)
      return x

    def outer_fn_unrolled(x):
      for _ in range(num_layers):
        x = inner_fn(x)
      return x

    def outer_fn_layer_stack(x):
      stack = layer_stack.layer_stack(num_layers, unroll=unroll)(inner_fn)
      return stack(x)

    unrolled_fn = hk.transform(outer_fn_unrolled)
    layer_stack_fn = hk.transform(outer_fn_layer_stack)

    x = jax.random.uniform(jax.random.PRNGKey(0), [10, 256, 100])

    rng_init = jax.random.PRNGKey(42)

    params = layer_stack_fn.init(rng_init, x)

    sliced_params = _slice_layers_params(params)

    unrolled_grad = jax.grad(
        lambda p, x: jnp.mean(unrolled_fn.apply(p, None, x)))(sliced_params, x)
    layer_stack_grad = jax.grad(
        lambda p, x: jnp.mean(layer_stack_fn.apply(p, None, x)))(params, x)

    assert_fn = functools.partial(
        np.testing.assert_allclose, atol=1e-4, rtol=1e-4)

    jax.tree_multimap(assert_fn, unrolled_grad,
                      _slice_layers_params(layer_stack_grad))

  def test_random(self):
    """Random numbers should be handled correctly."""
    n = 100

    @hk.transform
    @layer_stack.layer_stack(n)
    def add_random(x):
      x = x + jax.random.normal(hk.next_rng_key())
      return x

    # Evaluate a bunch of times
    key, *keys = jax.random.split(jax.random.PRNGKey(7), 1024 + 1)
    params = add_random.init(key, 0.)
    apply_fn = jax.jit(add_random.apply)
    values = [apply_fn(params, key, 0.) for key in keys]

    # Should be roughly N(0, sqrt(n))
    cdf = scipy.stats.norm(scale=np.sqrt(n)).cdf
    _, p = scipy.stats.kstest(values, cdf)
    self.assertLess(0.3, p)

  def test_threading(self):
    """Test @layer_stack when the function gets per-layer state."""
    n = 5

    @layer_stack.layer_stack(n, with_state=True)
    def f(x, y):
      x = x + y * jax.nn.one_hot(y, len(x)) / 10
      return x, 2 * y

    @hk.without_apply_rng
    @hk.transform
    def g(x, ys):
      x, zs = f(x, ys)
      # Check here to catch issues at init time
      self.assertEqual(zs.shape, (n,))
      return x, zs

    rng = jax.random.PRNGKey(7)
    x = np.zeros(n)
    ys = np.arange(n).astype(np.float32)
    params = g.init(rng, x, ys)
    x, zs = g.apply(params, x, ys)
    self.assertTrue(np.allclose(x, [0, .1, .2, .3, .4]))
    self.assertTrue(np.all(zs == 2 * ys))

  def test_nested_stacks(self):
    def stack_fn(x):
      def layer_fn(x):
        return hk.Linear(100)(x)

      outer_fn = layer_stack.layer_stack(10)(layer_fn)

      layer_outer = layer_stack.layer_stack(20)(outer_fn)
      return layer_outer(x)

    hk_mod = hk.transform(stack_fn)
    apply_rng, init_rng = jax.random.split(jax.random.PRNGKey(0))

    params = hk_mod.init(init_rng, jnp.zeros([10, 100]))

    hk_mod.apply(params, apply_rng, jnp.zeros([10, 100]))

    p, = params.values()

    assert p['w'].shape == (10, 20, 100, 100)
    assert p['b'].shape == (10, 20, 100)

  def test_with_state_multi_args(self):
    """Test layer_stack with state with multiple arguments."""
    width = 4
    batch_size = 5
    stack_height = 3

    def f_with_multi_args(x, a, b):
      return hk.Linear(
          width, w_init=hk.initializers.Constant(
              jnp.eye(width)))(x) * a + b, None

    @hk.without_apply_rng
    @hk.transform
    def hk_fn(x):
      return layer_stack.layer_stack(
          stack_height,
          with_state=True)(f_with_multi_args)(x, jnp.full([stack_height], 2.),
                                              jnp.ones([stack_height]))

    x = jnp.zeros([batch_size, width])
    key_seq = hk.PRNGSequence(19)
    params = hk_fn.init(next(key_seq), x)
    output, z = hk_fn.apply(params, x)
    self.assertIsNone(z)
    self.assertEqual(output.shape, (batch_size, width))
    np.testing.assert_equal(output, np.full([batch_size, width], 7.))

  def test_with_container_state(self):
    width = 2
    batch_size = 2
    stack_height = 3

    def f_with_container_state(x):
      hk_layer = hk.Linear(
          width, w_init=hk.initializers.Constant(jnp.eye(width)))
      layer_output = hk_layer(x)
      layer_state = {
          'raw_output': layer_output,
          'output_projection': jnp.sum(layer_output)
      }
      return layer_output + jnp.ones_like(layer_output), layer_state

    @hk.without_apply_rng
    @hk.transform
    def hk_fn(x):
      return layer_stack.layer_stack(
          stack_height,
          with_state=True)(f_with_container_state)(x)

    x = jnp.zeros([batch_size, width])
    key_seq = hk.PRNGSequence(19)
    params = hk_fn.init(next(key_seq), x)
    output, z = hk_fn.apply(params, x)
    self.assertEqual(z['raw_output'].shape, (stack_height, batch_size, width))
    self.assertEqual(output.shape, (batch_size, width))
    self.assertEqual(z['output_projection'].shape, (stack_height,))
    np.testing.assert_equal(np.sum(z['output_projection']), np.array(12.))
    np.testing.assert_equal(
        np.all(z['raw_output'] == np.array([0., 1., 2.])[..., None, None]),
        np.array(True))


if __name__ == '__main__':
  absltest.main()
