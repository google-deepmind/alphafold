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

"""A Python wrapper for hmmsearch - search profile against a sequence db."""

import os
import subprocess
from typing import Optional, Sequence

from absl import logging
from alphafold.data import parsers
from alphafold.data.tools import hmmbuild
from alphafold.data.tools import utils
# Internal import (7716).


class Hmmsearch(object):
  """Python wrapper of the hmmsearch binary."""

  def __init__(self,
               *,
               binary_path: str,
               hmmbuild_binary_path: str,
               database_path: str,
               cpus: int = 8,
               flags: Optional[Sequence[str]] = None):
    """Initializes the Python hmmsearch wrapper.

    Args:
      binary_path: The path to the hmmsearch executable.
      hmmbuild_binary_path: The path to the hmmbuild executable. Used to build
        an hmm from an input a3m.
      database_path: The path to the hmmsearch database (FASTA format).
      flags: List of flags to be used by hmmsearch.

    Raises:
      RuntimeError: If hmmsearch binary not found within the path.
    """
    self.cpus = cpus
    self.binary_path = binary_path
    self.hmmbuild_runner = hmmbuild.Hmmbuild(binary_path=hmmbuild_binary_path)
    self.database_path = database_path
    if flags is None:
      # Default hmmsearch run settings.
      flags = ['--F1', '0.1',
               '--F2', '0.1',
               '--F3', '0.1',
               '--incE', '100',
               '-E', '100',
               '--domE', '100',
               '--incdomE', '100']
    self.flags = flags

    if not os.path.exists(self.database_path):
      logging.error('Could not find hmmsearch database %s', database_path)
      raise ValueError(f'Could not find hmmsearch database {database_path}')

  @property
  def output_format(self) -> str:
    return 'sto'

  @property
  def input_format(self) -> str:
    return 'sto'

  def query(self, msa_sto: str) -> str:
    """Queries the database using hmmsearch using a given stockholm msa."""
    hmm = self.hmmbuild_runner.build_profile_from_sto(msa_sto,
                                                      model_construction='hand')
    return self.query_with_hmm(hmm)

  def query_with_hmm(self, hmm: str) -> str:
    """Queries the database using hmmsearch using a given hmm."""
    with utils.tmpdir_manager() as query_tmp_dir:
      hmm_input_path = os.path.join(query_tmp_dir, 'query.hmm')
      out_path = os.path.join(query_tmp_dir, 'output.sto')
      with open(hmm_input_path, 'w') as f:
        f.write(hmm)

      cmd = [
          self.binary_path,
          '--noali',  # Don't include the alignment in stdout.
          '--cpu', str(self.cpus)
      ]
      # If adding flags, we have to do so before the output and input:
      if self.flags:
        cmd.extend(self.flags)
      cmd.extend([
          '-A', out_path,
          hmm_input_path,
          self.database_path,
      ])

      logging.info('Launching sub-process %s', cmd)
      process = subprocess.Popen(
          cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      with utils.timing(
          f'hmmsearch ({os.path.basename(self.database_path)}) query'):
        stdout, stderr = process.communicate()
        retcode = process.wait()

      if retcode:
        raise RuntimeError(
            'hmmsearch failed:\nstdout:\n%s\n\nstderr:\n%s\n' % (
                stdout.decode('utf-8'), stderr.decode('utf-8')))

      with open(out_path) as f:
        out_msa = f.read()

    return out_msa

  def get_template_hits(self,
                        output_string: str,
                        input_sequence: str) -> Sequence[parsers.TemplateHit]:
    """Gets parsed template hits from the raw string output by the tool."""
    a3m_string = parsers.convert_stockholm_to_a3m(output_string,
                                                  remove_first_row_gaps=False)
    template_hits = parsers.parse_hmmsearch_a3m(
        query_sequence=input_sequence,
        a3m_string=a3m_string,
        skip_first=False)
    return template_hits
