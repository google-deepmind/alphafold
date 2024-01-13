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

"""Library to run Jackhmmer from Python."""

from concurrent import futures
import glob
import os
import subprocess
from typing import Any, Callable, Mapping, Optional, Sequence
from urllib import request

from absl import logging

from alphafold.data import parsers
from alphafold.data.tools import utils
# Internal import (7716).


class Jackhmmer:
  """Python wrapper of the Jackhmmer binary."""

  def __init__(self,
               *,
               binary_path: str,
               database_path: str,
               n_cpu: int = 8,
               n_iter: int = 1,
               e_value: float = 0.0001,
               z_value: Optional[int] = None,
               get_tblout: bool = False,
               filter_f1: float = 0.0005,
               filter_f2: float = 0.00005,
               filter_f3: float = 0.0000005,
               incdom_e: Optional[float] = None,
               dom_e: Optional[float] = None,
               num_streamed_chunks: Optional[int] = None,
               streaming_callback: Optional[Callable[[int], None]] = None):
    """Initializes the Python Jackhmmer wrapper.

    Args:
      binary_path: The path to the jackhmmer executable.
      database_path: The path to the jackhmmer database (FASTA format).
      n_cpu: The number of CPUs to give Jackhmmer.
      n_iter: The number of Jackhmmer iterations.
      e_value: The E-value, see Jackhmmer docs for more details.
      z_value: The Z-value, see Jackhmmer docs for more details.
      get_tblout: Whether to save tblout string.
      filter_f1: MSV and biased composition pre-filter, set to >1.0 to turn off.
      filter_f2: Viterbi pre-filter, set to >1.0 to turn off.
      filter_f3: Forward pre-filter, set to >1.0 to turn off.
      incdom_e: Domain e-value criteria for inclusion of domains in MSA/next
        round.
      dom_e: Domain e-value criteria for inclusion in tblout.
      num_streamed_chunks: Number of database chunks to stream over.
      streaming_callback: Callback function run after each chunk iteration with
        the iteration number as argument.
    """
    self.binary_path = binary_path
    self.database_path = database_path
    self.num_streamed_chunks = num_streamed_chunks

    if not os.path.exists(self.database_path) and num_streamed_chunks is None:
      logging.error('Could not find Jackhmmer database %s', database_path)
      raise ValueError(f'Could not find Jackhmmer database {database_path}')

    self.n_cpu = n_cpu
    self.n_iter = n_iter
    self.e_value = e_value
    self.z_value = z_value
    self.filter_f1 = filter_f1
    self.filter_f2 = filter_f2
    self.filter_f3 = filter_f3
    self.incdom_e = incdom_e
    self.dom_e = dom_e
    self.get_tblout = get_tblout
    self.streaming_callback = streaming_callback

  def _query_chunk(self,
                   input_fasta_path: str,
                   database_path: str,
                   max_sequences: Optional[int] = None) -> Mapping[str, Any]:
    """Queries the database chunk using Jackhmmer."""
    with utils.tmpdir_manager() as query_tmp_dir:
      sto_path = os.path.join(query_tmp_dir, 'output.sto')

      # The F1/F2/F3 are the expected proportion to pass each of the filtering
      # stages (which get progressively more expensive), reducing these
      # speeds up the pipeline at the expensive of sensitivity.  They are
      # currently set very low to make querying Mgnify run in a reasonable
      # amount of time.
      cmd_flags = [
          # Don't pollute stdout with Jackhmmer output.
          '-o', '/dev/null',
          '-A', sto_path,
          '--noali',
          '--F1', str(self.filter_f1),
          '--F2', str(self.filter_f2),
          '--F3', str(self.filter_f3),
          '--incE', str(self.e_value),
          # Report only sequences with E-values <= x in per-sequence output.
          '-E', str(self.e_value),
          '--cpu', str(self.n_cpu),
          '-N', str(self.n_iter)
      ]
      if self.get_tblout:
        tblout_path = os.path.join(query_tmp_dir, 'tblout.txt')
        cmd_flags.extend(['--tblout', tblout_path])

      if self.z_value:
        cmd_flags.extend(['-Z', str(self.z_value)])

      if self.dom_e is not None:
        cmd_flags.extend(['--domE', str(self.dom_e)])

      if self.incdom_e is not None:
        cmd_flags.extend(['--incdomE', str(self.incdom_e)])

      cmd = [self.binary_path] + cmd_flags + [input_fasta_path,
                                              database_path]

      logging.info('Launching subprocess "%s"', ' '.join(cmd))
      process = subprocess.Popen(
          cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      with utils.timing(
          f'Jackhmmer ({os.path.basename(database_path)}) query'):
        _, stderr = process.communicate()
        retcode = process.wait()

      if retcode:
        raise RuntimeError(
            'Jackhmmer failed\nstderr:\n%s\n' % stderr.decode('utf-8'))

      # Get e-values for each target name
      tbl = ''
      if self.get_tblout:
        with open(tblout_path) as f:
          tbl = f.read()

      if max_sequences is None:
        with open(sto_path) as f:
          sto = f.read()
      else:
        sto = parsers.truncate_stockholm_msa(sto_path, max_sequences)

    raw_output = dict(
        sto=sto,
        tbl=tbl,
        stderr=stderr,
        n_iter=self.n_iter,
        e_value=self.e_value)

    return raw_output

  def query(self,
            input_fasta_path: str,
            max_sequences: Optional[int] = None) -> Sequence[Mapping[str, Any]]:
    """Queries the database using Jackhmmer."""
    return self.query_multiple([input_fasta_path], max_sequences)[0]

  def query_multiple(
      self,
      input_fasta_paths: Sequence[str],
      max_sequences: Optional[int] = None,
    ) -> Sequence[Sequence[Mapping[str, Any]]]:
    """Queries the database for multiple queries using Jackhmmer."""
    if self.num_streamed_chunks is None:
      single_chunk_results = []
      for input_fasta_path in input_fasta_paths:
        single_chunk_results.append([self._query_chunk(
            input_fasta_path, self.database_path, max_sequences)])
      return single_chunk_results

    db_basename = os.path.basename(self.database_path)
    db_remote_chunk = lambda db_idx: f'{self.database_path}.{db_idx}'
    db_local_chunk = lambda db_idx: f'/tmp/ramdisk/{db_basename}.{db_idx}'

    # Remove existing files to prevent OOM
    for f in glob.glob(db_local_chunk('[0-9]*')):
      try:
        os.remove(f)
      except OSError:
        print(f'OSError while deleting {f}')

    # Download the (i+1)-th chunk while Jackhmmer is running on the i-th chunk
    with futures.ThreadPoolExecutor(max_workers=2) as executor:
      chunked_outputs = [[] for _ in range(len(input_fasta_paths))]
      for i in range(1, self.num_streamed_chunks + 1):
        # Copy the chunk locally
        if i == 1:
          future = executor.submit(
              request.urlretrieve, db_remote_chunk(i), db_local_chunk(i))
        if i < self.num_streamed_chunks:
          next_future = executor.submit(
              request.urlretrieve, db_remote_chunk(i+1), db_local_chunk(i+1))

        # Run Jackhmmer with the chunk
        future.result()
        for fasta_index, input_fasta_path in enumerate(input_fasta_paths):
          chunked_outputs[fasta_index].append(self._query_chunk(
              input_fasta_path, db_local_chunk(i), max_sequences))
        # Remove the local copy of the chunk
        os.remove(db_local_chunk(i))
        # Do not set next_future for the last chunk so that this works even for
        # databases with only 1 chunk.
        if i < self.num_streamed_chunks:
          future = next_future
        if self.streaming_callback:
          self.streaming_callback(i)
    return chunked_outputs
