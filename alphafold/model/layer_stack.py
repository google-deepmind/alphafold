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

"""Function to stack repeats of a layer function without shared parameters."""

import collections
import contextlib
import functools
import inspect
from typing import Any, Callable, Optional, Tuple, Union

import haiku as hk
import jax
import jax.numpy as jnp

LayerStackCarry = collections.namedtuple('LayerStackCarry', ['x', 'rng'])
LayerStackScanned = collections.namedtuple('LayerStackScanned',
                                           ['i', 'args_ys'])

# WrappedFn should take in arbitrarily nested `jnp.ndarray`, and return the
# exact same type. We cannot express this with `typing`. So we just use it
# to inform the user. In reality, the typing below will accept anything.
NestedArray = Any
WrappedFn = Callable[..., Union[NestedArray, Tuple[NestedArray]]]


def _check_no_varargs(f):
  if list(inspect.signature(
      f).parameters.values())[0].kind == inspect.Parameter.VAR_POSITIONAL:
    raise ValueError(
        'The function `f` should not have any `varargs` (that is *args) '
        'argument. Instead, it should only use explicit positional'
        'arguments.')


@contextlib.contextmanager
def nullcontext():
  yield


def maybe_with_rng(key):
  if key is not None:
    return hk.with_rng(key)
  else:
    return nullcontext()


def maybe_fold_in(key, data):
  if key is not None:
    return jax.random.fold_in(key, data)
  else:
    return None


class _LayerStack(hk.Module):
  """Module to compose parameterized functions, implemented as a scan."""

  def __init__(self,
               count: int,
               unroll: int,
               name: Optional[str] = None):
    """Iterate a function `f` `count` times, with non-shared parameters."""
    super().__init__(name=name)
    self._count = count
    self._unroll = unroll

  def __call__(self, x, *args_ys):
    count = self._count
    if hk.running_init():
      # At initialization time, we run just one layer but add an extra first
      # dimension to every initialized tensor, making sure to use different
      # random keys for different slices.
      def creator(next_creator, shape, dtype, init, context):
        del context

        def multi_init(shape, dtype):
          assert shape[0] == count
          key = hk.maybe_next_rng_key()

          def rng_context_init(slice_idx):
            slice_key = maybe_fold_in(key, slice_idx)
            with maybe_with_rng(slice_key):
              return init(shape[1:], dtype)

          return jax.vmap(rng_context_init)(jnp.arange(count))

        return next_creator((count,) + tuple(shape), dtype, multi_init)

      def getter(next_getter, value, context):
        trailing_dims = len(context.original_shape) + 1
        sliced_value = jax.lax.index_in_dim(
            value, index=0, axis=value.ndim - trailing_dims, keepdims=False)
        return next_getter(sliced_value)

      with hk.experimental.custom_creator(
          creator), hk.experimental.custom_getter(getter):
        if len(args_ys) == 1 and args_ys[0] is None:
          args0 = (None,)
        else:
          args0 = [
              jax.lax.dynamic_index_in_dim(ys, 0, keepdims=False)
              for ys in args_ys
          ]
        x, z = self._call_wrapped(x, *args0)
        if z is None:
          return x, z

        # Broadcast state to hold each layer state.
        def broadcast_state(layer_state):
          return jnp.broadcast_to(
              layer_state, [count,] + list(layer_state.shape))
        zs = jax.tree_util.tree_map(broadcast_state, z)
        return x, zs
    else:
      # Use scan during apply, threading through random seed so that it's
      # unique for each layer.
      def layer(carry: LayerStackCarry, scanned: LayerStackScanned):
        rng = carry.rng

        def getter(next_getter, value, context):
          # Getter slices the full param at the current loop index.
          trailing_dims = len(context.original_shape) + 1
          assert value.shape[value.ndim - trailing_dims] == count, (
              f'Attempting to use a parameter stack of size '
              f'{value.shape[value.ndim - trailing_dims]} for a LayerStack of '
              f'size {count}.')

          sliced_value = jax.lax.dynamic_index_in_dim(
              value, scanned.i, axis=value.ndim - trailing_dims, keepdims=False)
          return next_getter(sliced_value)

        with hk.experimental.custom_getter(getter):
          if rng is None:
            out_x, z = self._call_wrapped(carry.x, *scanned.args_ys)
          else:
            rng, rng_ = jax.random.split(rng)
            with hk.with_rng(rng_):
              out_x, z = self._call_wrapped(carry.x, *scanned.args_ys)
        return LayerStackCarry(x=out_x, rng=rng), z

      carry = LayerStackCarry(x=x, rng=hk.maybe_next_rng_key())
      scanned = LayerStackScanned(i=jnp.arange(count, dtype=jnp.int32),
                                  args_ys=args_ys)

      carry, zs = hk.scan(
          layer, carry, scanned, length=count, unroll=self._unroll)
      return carry.x, zs

  def _call_wrapped(self,
                    x: jnp.ndarray,
                    *args,
                    ) -> Tuple[jnp.ndarray, Optional[jnp.ndarray]]:
    raise NotImplementedError()


class _LayerStackNoState(_LayerStack):
  """_LayerStack impl with no per-layer state provided to the function."""

  def __init__(self,
               f: WrappedFn,
               count: int,
               unroll: int,
               name: Optional[str] = None):
    super().__init__(count=count, unroll=unroll, name=name)
    _check_no_varargs(f)
    self._f = f

  @hk.transparent
  def _call_wrapped(self, args, y):
    del y
    ret = self._f(*args)
    if len(args) == 1:
      # If the function takes a single argument, the wrapped function receives
      # a tuple of length 1, and therefore it must return a tuple of length 1.
      ret = (ret,)
    return ret, None


class _LayerStackWithState(_LayerStack):
  """_LayerStack impl with per-layer state provided to the function."""

  def __init__(self,
               f: WrappedFn,
               count: int,
               unroll: int,
               name: Optional[str] = None):
    super().__init__(count=count, unroll=unroll, name=name)
    self._f = f

  @hk.transparent
  def _call_wrapped(self, x, *args):
    return self._f(x, *args)


def layer_stack(num_layers: int,
                with_state=False,
                unroll: int = 1,
                name: Optional[str] = None):
  """Utility to wrap a Haiku function and recursively apply it to an input.

  A function is valid if it uses only explicit position parameters, and
  its return type matches its input type. The position parameters can be
  arbitrarily nested structures with `jnp.ndarray` at the leaf nodes. Note
  that kwargs are not supported, neither are functions with variable number
  of parameters (specified by `*args`).

  If `with_state=False` then the new, wrapped function can be understood as
  performing the following:
  ```
  for i in range(num_layers):
    x = f(x)
  return x
  ```

  And if `with_state=True`, assuming `f` takes two arguments on top of `x`:
  ```
  for i in range(num_layers):
    x, zs[i] = f(x, ys_0[i], ys_1[i])
  return x, zs
  ```
  The code using `layer_stack` for the above function would be:
  ```
  def f(x, y_0, y_1):
    ...
    return new_x, z
  x, zs = layer_stack.layer_stack(num_layers,
                                  with_state=True)(f)(x, ys_0, ys_1)
  ```

  Crucially, any parameters created inside `f` will not be shared across
  iterations.

  Args:
    num_layers: The number of times to iterate the wrapped function.
    with_state: Whether or not to pass per-layer state to the wrapped function.
    unroll: the unroll used by `scan`.
    name: Name of the Haiku context.

  Returns:
    Callable that will produce a layer stack when called with a valid function.
  """
  def iterate(f):
    if with_state:
      @functools.wraps(f)
      def wrapped(x, *args):
        for ys in args:
          assert ys.shape[0] == num_layers
        return _LayerStackWithState(
            f, num_layers, unroll=unroll, name=name)(x, *args)
    else:
      _check_no_varargs(f)
      @functools.wraps(f)
      def wrapped(*args):
        ret = _LayerStackNoState(
            f, num_layers, unroll=unroll, name=name)(args, None)[0]
        if len(args) == 1:
          # If the function takes a single argument, we must also return a
          # single value, and not a tuple of length 1.
          ret = ret[0]
        return ret

    return wrapped
  return iterate
