# Copyright 2025 DeepMind Technologies Limited
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

"""Fixture for pytest.

This is needed to parse the absl flags before running the test.
"""

import sys

from absl import flags
import pytest


@pytest.fixture(scope="session", autouse=True)
def initialize_absl_flags(request):
  del request
  # Parse any flags that make sense to absl as absl flags.
  flags.FLAGS(sys.argv, known_only=True)
