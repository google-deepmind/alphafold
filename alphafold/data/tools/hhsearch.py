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

"""Library to run HHsearch from Python."""

import glob
import os
import subprocess
from typing import Sequence

from absl import logging

from alphafold.data import parsers
from alphafold.data.tools import utils
# Internal import (7716).


class HHSearch:
  """Python wrapper of the HHsearch binary."""

  def __init__(self,
               *,
               binary_path: str,
               databases: Sequence[str],
               maxseq: int = 1_000_000,
               cpus: int = 8):
    """Initializes the Python HHsearch wrapper.

    Args:
      binary_path: The path to the HHsearch executable.
      databases: A sequence of HHsearch database paths. This should be the
        common prefix for the database files (i.e. up to but not including
        _hhm.ffindex etc.)
      maxseq: The maximum number of rows in an input alignment. Note that this
        parameter is only supported in HHBlits version 3.1 and higher.

    Raises:
      RuntimeError: If HHsearch binary not found within the path.
    """
    self.cpus = cpus
    self.binary_path = binary_path
    self.databases = databases
    self.maxseq = maxseq

    for database_path in self.databases:
      if not glob.glob(database_path + '_*'):
        logging.error('Could not find HHsearch database %s', database_path)
        raise ValueError(f'Could not find HHsearch database {database_path}')

  @property
  def output_format(self) -> str:
    return 'hhr'

  @property
  def input_format(self) -> str:
    return 'a3m'

  def query(self, a3m: str) -> str:
    """Queries the database using HHsearch using a given a3m."""
    with utils.tmpdir_manager() as query_tmp_dir:
      input_path = os.path.join(query_tmp_dir, 'query.a3m')
      hhr_path = os.path.join(query_tmp_dir, 'output.hhr')
      with open(input_path, 'w') as f:
        f.write(a3m)

      db_cmd = []
      for db_path in self.databases:
        db_cmd.append('-d')
        db_cmd.append(db_path)
      cmd = [self.binary_path,
             '-i', input_path,
             '-o', hhr_path,
             '-maxseq', str(self.maxseq),
             '-cpu', str(self.cpus)
             ] + db_cmd

      logging.info('Launching subprocess "%s"', ' '.join(cmd))
      process = subprocess.Popen(
          cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      with utils.timing('HHsearch query'):
        stdout, stderr = process.communicate()
        retcode = process.wait()

      if retcode:
        # Stderr is truncated to prevent proto size errors in Beam.
        raise RuntimeError(
            'HHSearch failed:\nstdout:\n%s\n\nstderr:\n%s\n' % (
                stdout.decode('utf-8'), stderr[:100_000].decode('utf-8')))

      with open(hhr_path) as f:
        hhr = f.read()
    return hhr

  def get_template_hits(self,
                        output_string: str,
                        input_sequence: str) -> Sequence[parsers.TemplateHit]:
    """Gets parsed template hits from the raw string output by the tool."""
    del input_sequence  # Used by hmmseach but not needed for hhsearch.
    return parsers.parse_hhr(output_string)
