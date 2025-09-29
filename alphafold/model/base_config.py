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

"""Config for the protein folding model and experiment."""

from collections.abc import Mapping
import contextlib
import copy
import dataclasses
import types
import typing
from typing import Any, ClassVar, Iterator, TypeVar


_T = TypeVar('_T')
_ConfigT = TypeVar('_ConfigT', bound='BaseConfig')


def _strip_optional(t: type[Any]) -> type[Any]:
  """Transforms type annotations of the form `T | None` to `T`."""
  if typing.get_origin(t) in (typing.Union, types.UnionType):
    args = set(typing.get_args(t)) - {types.NoneType}
    if len(args) == 1:
      return args.pop()
  return t


_NO_UPDATE = object()


class _Autocreate:

  def __init__(self, **defaults: Any):
    self.defaults = defaults


def autocreate(**defaults: Any) -> Any:
  """Marks a field as having a default factory derived from its type."""
  return _Autocreate(**defaults)


def _clone_field(
    field: dataclasses.Field[_T], new_default: _T
) -> dataclasses.Field[_T]:
  if new_default is _NO_UPDATE:
    return copy.copy(field)
  return dataclasses.field(
      default=new_default,
      init=True,
      kw_only=True,
      repr=field.repr,
      hash=field.hash,
      compare=field.compare,
      metadata=field.metadata,
  )


@typing.dataclass_transform()
class ConfigMeta(type):
  """Metaclass that synthesizes a __post_init__ that coerces dicts to Config subclass instances."""

  def __new__(mcs, name, bases, classdict):
    cls = super().__new__(mcs, name, bases, classdict)

    def _coercable_fields(self) -> Mapping[str, tuple[ConfigMeta, Any]]:
      type_hints = typing.get_type_hints(self.__class__)
      fields = dataclasses.fields(self.__class__)
      field_to_type_and_default = {
          field.name: (_strip_optional(type_hints[field.name]), field.default)
          for field in fields
      }
      coercable_fields = {
          f: t
          for f, t in field_to_type_and_default.items()
          if issubclass(type(t[0]), ConfigMeta)
      }
      return coercable_fields

    cls._coercable_fields = property(_coercable_fields)

    old_post_init = getattr(cls, '__post_init__', None)

    def _post_init(self) -> None:
      # Use get_type_hints instead of Field.type to ensure that forward
      # references are resolved.
      for field_name, (
          field_type,
          field_default,
      ) in self._coercable_fields.items():  # pylint: disable=protected-access
        field_value = getattr(self, field_name)
        if field_value is None:
          continue
        try:
          match field_value:
            case _Autocreate():
              # Construct from field defaults.
              setattr(self, field_name, field_type(**field_value.defaults))
            case Mapping():
              # Field value is not yet a `Config` instance; Assume we can create
              # one by splatting keys and values.
              args = {}
              # Apply default args first, if present.
              if isinstance(field_default, _Autocreate):
                args.update(field_default.defaults)
              args.update(field_value)
              setattr(self, field_name, field_type(**args))
            case _:
              pass
        except TypeError as e:
          raise TypeError(
              f'Failure while coercing field {field_name!r} of'
              f' {self.__class__.__qualname__}'
          ) from e
      if old_post_init:
        old_post_init(self)

    cls.__post_init__ = _post_init
    return dataclasses.dataclass(kw_only=True)(cls)


class BaseConfig(metaclass=ConfigMeta):
  """Config base class.

  Subclassing BaseConfig automatically makes the subclass a kw_only dataclass
  with a `__post_init__` that coerces Config-subclass field values from mappings
  to instances of the right type.
  """

  # Provided by dataclasses.make_dataclass
  __dataclass_fields__: ClassVar[dict[str, dataclasses.Field[Any]]]
  _is_frozen: ClassVar[bool] = dataclasses.field(
      default=False, init=False, repr=False
  )

  # Overridden by metaclass
  @property
  def _coercable_fields(self) -> Mapping[str, tuple[type['BaseConfig'], Any]]:
    return {}

  def as_dict(self, include_none: bool = True) -> Mapping[str, Any]:
    """Returns a dict representation of the config.

    Args:
      include_none: Whether to include fields with value None.
    """
    result = dataclasses.asdict(self)
    for field_name in self._coercable_fields:
      field_value = getattr(self, field_name, None)
      if isinstance(field_value, BaseConfig):
        result[field_name] = field_value.as_dict(include_none)
    return (
        result
        if include_none
        else {k: v for k, v in result.items() if v is not None}
    )

  def __setattr__(self, name: str, value: Any) -> None:
    if getattr(self, '_is_frozen', False) and name != '_is_frozen':
      # If we are frozen, raise an error
      raise dataclasses.FrozenInstanceError(
          f"Cannot assign to field '{name}'; instance is frozen."
      )

    # If not frozen, set the attribute normally
    super().__setattr__(name, value)

  def _toggle_freeze(self, frozen: bool) -> None:
    """Toggles the frozen state of the config and all subconfigs."""
    self._is_frozen = frozen
    for field_name in self._coercable_fields:
      field_value = getattr(self, field_name, None)
      if isinstance(field_value, BaseConfig):
        field_value._toggle_freeze(frozen)

  def freeze(self) -> None:
    """Freezes the config and all subconfigs to prevent further changes."""
    self._toggle_freeze(True)

  @contextlib.contextmanager
  def unfreeze(self: _ConfigT) -> Iterator[_ConfigT]:
    """A context manager to temporarily unfreeze the config."""
    was_frozen = self._is_frozen
    self._toggle_freeze(False)
    try:
      yield self
    finally:
      if was_frozen:
        self._toggle_freeze(True)
