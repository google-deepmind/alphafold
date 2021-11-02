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
"""Class decorator to represent (nested) struct of arrays."""

import dataclasses

import jax


def get_item(instance, key):
  sliced = {}
  for field in get_array_fields(instance):
    num_trailing_dims = field.metadata.get('num_trailing_dims', 0)
    this_key = key
    if isinstance(key, tuple) and Ellipsis in this_key:
      this_key += (slice(None),) * num_trailing_dims
    sliced[field.name] = getattr(instance, field.name)[this_key]
  return dataclasses.replace(instance, **sliced)


@property
def get_shape(instance):
  """Returns Shape for given instance of dataclass."""
  first_field = dataclasses.fields(instance)[0]
  num_trailing_dims = first_field.metadata.get('num_trailing_dims', None)
  value = getattr(instance, first_field.name)
  if num_trailing_dims:
    return value.shape[:-num_trailing_dims]
  else:
    return value.shape


def get_len(instance):
  """Returns length for given instance of dataclass."""
  shape = instance.shape
  if shape:
    return shape[0]
  else:
    raise TypeError('len() of unsized object')  # Match jax.numpy behavior.


@property
def get_dtype(instance):
  """Returns Dtype for given instance of dataclass."""
  fields = dataclasses.fields(instance)
  sets_dtype = [
      field.name for field in fields if field.metadata.get('sets_dtype', False)
  ]
  if sets_dtype:
    assert len(sets_dtype) == 1, 'at most field can set dtype'
    field_value = getattr(instance, sets_dtype[0])
  elif instance.same_dtype:
    field_value = getattr(instance, fields[0].name)
  else:
    # Should this be Value Error?
    raise AttributeError('Trying to access Dtype on Struct of Array without'
                         'either "same_dtype" or field setting dtype')

  if hasattr(field_value, 'dtype'):
    return field_value.dtype
  else:
    # Should this be Value Error?
    raise AttributeError(f'field_value {field_value} does not have dtype')


def replace(instance, **kwargs):
  return dataclasses.replace(instance, **kwargs)


def post_init(instance):
  """Validate instance has same shapes & dtypes."""
  array_fields = get_array_fields(instance)
  arrays = list(get_array_fields(instance, return_values=True).values())
  first_field = array_fields[0]
  # These slightly weird constructions about checking whether the leaves are
  # actual arrays is since e.g. vmap internally relies on being able to
  # construct pytree's with object() as leaves, this would break the checking
  # as such we are only validating the object when the entries in the dataclass
  # Are arrays or other dataclasses of arrays.
  try:
    dtype = instance.dtype
  except AttributeError:
    dtype = None
  if dtype is not None:
    first_shape = instance.shape
    for array, field in zip(arrays, array_fields):
      field_shape = array.shape
      num_trailing_dims = field.metadata.get('num_trailing_dims', None)
      if num_trailing_dims:
        array_shape = array.shape
        field_shape = array_shape[:-num_trailing_dims]
        msg = (f'field {field} should have number of trailing dims'
               ' {num_trailing_dims}')
        assert len(array_shape) == len(first_shape) + num_trailing_dims, msg
      else:
        field_shape = array.shape

      shape_msg = (f"Stripped Shape {field_shape} of field {field} doesn't "
                   f"match shape {first_shape} of field {first_field}")
      assert field_shape == first_shape, shape_msg

      field_dtype = array.dtype

      allowed_metadata_dtypes = field.metadata.get('allowed_dtypes', [])
      if allowed_metadata_dtypes:
        msg = f'Dtype is {field_dtype} but must be in {allowed_metadata_dtypes}'
        assert field_dtype in allowed_metadata_dtypes, msg

      if 'dtype' in field.metadata:
        target_dtype = field.metadata['dtype']
      else:
        target_dtype = dtype

      msg = f'Dtype is {field_dtype} but must be {target_dtype}'
      assert field_dtype == target_dtype, msg


def flatten(instance):
  """Flatten Struct of Array instance."""
  array_likes = list(get_array_fields(instance, return_values=True).values())
  flat_array_likes = []
  inner_treedefs = []
  num_arrays = []
  for array_like in array_likes:
    flat_array_like, inner_treedef = jax.tree_flatten(array_like)
    inner_treedefs.append(inner_treedef)
    flat_array_likes += flat_array_like
    num_arrays.append(len(flat_array_like))
  metadata = get_metadata_fields(instance, return_values=True)
  metadata = type(instance).metadata_cls(**metadata)
  return flat_array_likes, (inner_treedefs, metadata, num_arrays)


def make_metadata_class(cls):
  metadata_fields = get_fields(cls,
                               lambda x: x.metadata.get('is_metadata', False))
  metadata_cls = dataclasses.make_dataclass(
      cls_name='Meta' + cls.__name__,
      fields=[(field.name, field.type, field) for field in metadata_fields],
      frozen=True,
      eq=True)
  return metadata_cls


def get_fields(cls_or_instance, filterfn, return_values=False):
  fields = dataclasses.fields(cls_or_instance)
  fields = [field for field in fields if filterfn(field)]
  if return_values:
    return {
        field.name: getattr(cls_or_instance, field.name) for field in fields
    }
  else:
    return fields


def get_array_fields(cls, return_values=False):
  return get_fields(
      cls,
      lambda x: not x.metadata.get('is_metadata', False),
      return_values=return_values)


def get_metadata_fields(cls, return_values=False):
  return get_fields(
      cls,
      lambda x: x.metadata.get('is_metadata', False),
      return_values=return_values)


class StructOfArray:
  """Class Decorator for Struct Of Arrays."""

  def __init__(self, same_dtype=True):
    self.same_dtype = same_dtype

  def __call__(self, cls):
    cls.__array_ufunc__ = None
    cls.replace = replace
    cls.same_dtype = self.same_dtype
    cls.dtype = get_dtype
    cls.shape = get_shape
    cls.__len__ = get_len
    cls.__getitem__ = get_item
    cls.__post_init__ = post_init
    new_cls = dataclasses.dataclass(cls, frozen=True, eq=False)  # pytype: disable=wrong-keyword-args
    # pytree claims to require metadata to be hashable, not sure why,
    # But making derived dataclass that can just hold metadata
    new_cls.metadata_cls = make_metadata_class(new_cls)

    def unflatten(aux, data):
      inner_treedefs, metadata, num_arrays = aux
      array_fields = [field.name for field in get_array_fields(new_cls)]
      value_dict = {}
      array_start = 0
      for num_array, inner_treedef, array_field in zip(num_arrays,
                                                       inner_treedefs,
                                                       array_fields):
        value_dict[array_field] = jax.tree_unflatten(
            inner_treedef, data[array_start:array_start + num_array])
        array_start += num_array
      metadata_fields = get_metadata_fields(new_cls)
      for field in metadata_fields:
        value_dict[field.name] = getattr(metadata, field.name)

      return new_cls(**value_dict)

    jax.tree_util.register_pytree_node(
        nodetype=new_cls, flatten_func=flatten, unflatten_func=unflatten)
    return new_cls
