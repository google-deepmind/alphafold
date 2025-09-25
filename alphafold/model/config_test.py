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

import json

from absl.testing import absltest
from absl.testing import parameterized
from alphafold.model import config


class ConfigTest(parameterized.TestCase):

  @parameterized.parameters(config.CONFIG_DIFFS.keys())
  def test_config_dict_and_dataclass_agree(self, model_name):
    """Ensures model_config() and get_model_config() return same values."""
    config_dict_json = json.dumps(config.model_config(model_name).to_dict())
    config_dataclass_json = json.dumps(
        config.get_model_config(model_name).as_dict(include_none=False)
    )
    self.assertJsonEqual(config_dict_json, config_dataclass_json)


if __name__ == '__main__':
  absltest.main()
