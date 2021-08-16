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

"""A Python wrapper for hmmbuild - construct HMM profiles from MSA."""

import os
import re
import subprocess

from absl import logging
from alphafold.data.tools import utils
# Internal import (7716).


class Hmmbuild(object):
  """Python wrapper of the hmmbuild binary."""

  def __init__(self,
               *,
               binary_path: str,
               singlemx: bool = False):
    """Initializes the Python hmmbuild wrapper.

    Args:
      binary_path: The path to the hmmbuild executable.
      singlemx: Whether to use --singlemx flag. If True, it forces HMMBuild to
        just use a common substitution score matrix.

    Raises:
      RuntimeError: If hmmbuild binary not found within the path.
    """
    self.binary_path = binary_path
    self.singlemx = singlemx

  def build_profile_from_sto(self, sto: str, model_construction='fast') -> str:
    """Builds a HHM for the aligned sequences given as an A3M string.

    Args:
      sto: A string with the aligned sequences in the Stockholm format.
      model_construction: Whether to use reference annotation in the msa to
        determine consensus columns ('hand') or default ('fast').

    Returns:
      A string with the profile in the HMM format.

    Raises:
      RuntimeError: If hmmbuild fails.
    """
    return self._build_profile(sto, model_construction=model_construction)

  def build_profile_from_a3m(self, a3m: str) -> str:
    """Builds a HHM for the aligned sequences given as an A3M string.

    Args:
      a3m: A string with the aligned sequences in the A3M format.

    Returns:
      A string with the profile in the HMM format.

    Raises:
      RuntimeError: If hmmbuild fails.
    """
    lines = []
    for line in a3m.splitlines():
      if not line.startswith('>'):
        line = re.sub('[a-z]+', '', line)  # Remove inserted residues.
      lines.append(line + '\n')
    msa = ''.join(lines)
    return self._build_profile(msa, model_construction='fast')

  def _build_profile(self, msa: str, model_construction: str = 'fast') -> str:
    """Builds a HMM for the aligned sequences given as an MSA string.

    Args:
      msa: A string with the aligned sequences, in A3M or STO format.
      model_construction: Whether to use reference annotation in the msa to
        determine consensus columns ('hand') or default ('fast').

    Returns:
      A string with the profile in the HMM format.

    Raises:
      RuntimeError: If hmmbuild fails.
      ValueError: If unspecified arguments are provided.
    """
    if model_construction not in {'hand', 'fast'}:
      raise ValueError(f'Invalid model_construction {model_construction} - only'
                       'hand and fast supported.')

    with utils.tmpdir_manager(base_dir='/tmp') as query_tmp_dir:
      input_query = os.path.join(query_tmp_dir, 'query.msa')
      output_hmm_path = os.path.join(query_tmp_dir, 'output.hmm')

      with open(input_query, 'w') as f:
        f.write(msa)

      cmd = [self.binary_path]
      # If adding flags, we have to do so before the output and input:

      if model_construction == 'hand':
        cmd.append(f'--{model_construction}')
      if self.singlemx:
        cmd.append('--singlemx')
      cmd.extend([
          '--amino',
          output_hmm_path,
          input_query,
      ])

      logging.info('Launching subprocess %s', cmd)
      process = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)

      with utils.timing('hmmbuild query'):
        stdout, stderr = process.communicate()
        retcode = process.wait()
        logging.info('hmmbuild stdout:\n%s\n\nstderr:\n%s\n',
                     stdout.decode('utf-8'), stderr.decode('utf-8'))

      if retcode:
        raise RuntimeError('hmmbuild failed\nstdout:\n%s\n\nstderr:\n%s\n'
                           % (stdout.decode('utf-8'), stderr.decode('utf-8')))

      with open(output_hmm_path, encoding='utf-8') as f:
        hmm = f.read()

    return hmm
