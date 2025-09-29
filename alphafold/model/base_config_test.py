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

import dataclasses
from unittest import mock

from absl.testing import absltest
from alphafold.model import base_config
import jax


class InnerConfig(base_config.BaseConfig):
  a: int
  b: int = 10


class OuterConfig(base_config.BaseConfig):
  x: int
  z: InnerConfig
  optional_z: InnerConfig | None
  y: int = 11
  optional_z_default: InnerConfig | None = None
  z_requires_a: InnerConfig = base_config.autocreate()
  z_default: InnerConfig | None = base_config.autocreate(a=12)


class ModelConfigTest(absltest.TestCase):

  def _equal_at_path(self, path, a, b):
    self.assertEqual(a, b, f'trees differ at path {path}: {a} != {b}')

  def test_post_init_is_chained(self):
    post_init = mock.Mock()

    class Config(base_config.BaseConfig):
      x: int
      __post_init__ = post_init

    Config(x=1)
    post_init.assert_called_once()

  def test_config_is_dataclass(self):
    self.assertTrue(dataclasses.is_dataclass(OuterConfig))

  def test_nested_values_not_provided(self):
    with self.assertRaisesRegex(
        TypeError, r"Failure while coercing field 'z_requires_a' of OuterConfig"
    ):
      OuterConfig(x=5, z=InnerConfig(a=2), optional_z=None)

  def test_config_dict_escape_hatch(self):
    class Config(base_config.BaseConfig):
      x: int
      y: dict[str, int]

    conf = {'x': 1, 'y': {'z': 2}}
    conf2 = Config(**conf)
    self.assertIs(conf['y'], conf2.y)

  def test_config_from_dict(self):
    config = OuterConfig(**{
        'x': 5,
        'z': InnerConfig(a=2),
        'optional_z': InnerConfig(a=3),
        'optional_z_default': InnerConfig(a=4),
        'z_requires_a': InnerConfig(a=5, b=10),
    })
    expected = {
        'x': 5,
        'z': dict(a=2, b=10),
        'optional_z': dict(a=3, b=10),
        'y': 11,
        'optional_z_default': dict(a=4, b=10),
        'z_requires_a': dict(a=5, b=10),
        'z_default': dict(a=12, b=10),
    }
    jax.tree_util.tree_map_with_path(
        self._equal_at_path, config.as_dict(), expected
    )

  def test_create_config(self):
    config = OuterConfig(
        x=5,
        z=InnerConfig(a=2),
        optional_z=None,
        z_requires_a=InnerConfig(a=3),
        z_default=None,
    )
    expected = {
        'x': 5,
        'z': dict(a=2, b=10),
        'optional_z': None,
        'y': 11,
        'optional_z_default': None,
        'z_requires_a': dict(a=3, b=10),
        'z_default': None,
    }
    jax.tree_util.tree_map_with_path(
        self._equal_at_path, config.as_dict(), expected
    )

  def test_freeze(self):
    config = OuterConfig(
        x=5,
        z=InnerConfig(a=2),
        optional_z=None,
        z_requires_a=InnerConfig(a=3),
        z_default=None,
    )
    config.freeze()

    # Check that the config and all subconfigs are frozen.
    self.assertTrue(config._is_frozen)
    self.assertTrue(config.z._is_frozen)
    self.assertTrue(config.z_requires_a._is_frozen)

    # Check that we cannot modify the config.
    with self.assertRaises(dataclasses.FrozenInstanceError):
      config.x = 1
    with self.assertRaises(dataclasses.FrozenInstanceError):
      config.z.a = 1

  def test_unfreeze(self):
    config = OuterConfig(
        x=5,
        z=InnerConfig(a=2),
        optional_z=None,
        z_requires_a=InnerConfig(a=3),
        z_default=None,
    )
    config.freeze()

    # Check that we can modify the config within the unfrozen context.
    with config.unfreeze() as mutable_config:
      mutable_config.x = 1
      mutable_config.z.a = 1
    self.assertEqual(config.x, 1)
    self.assertEqual(config.z.a, 1)

    # Check that the config and all subconfigs are frozen again.
    self.assertTrue(config._is_frozen)
    self.assertTrue(config.z._is_frozen)
    self.assertTrue(config.z_requires_a._is_frozen)
    with self.assertRaises(dataclasses.FrozenInstanceError):
      config.x = 2
    with self.assertRaises(dataclasses.FrozenInstanceError):
      config.z.a = 2

    # Check that a config that was not frozen remains unfrozen.
    unfrozen_config = OuterConfig(
        x=5,
        z=InnerConfig(a=2),
        optional_z=None,
        z_requires_a=InnerConfig(a=3),
        z_default=None,
    )
    self.assertFalse(unfrozen_config._is_frozen)
    with unfrozen_config.unfreeze() as mutable_config:
      mutable_config.x = 1
    self.assertEqual(unfrozen_config.x, 1)
    self.assertFalse(unfrozen_config._is_frozen)


if __name__ == '__main__':
  absltest.main()
