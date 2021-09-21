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

"""Singularity launch script for Alphafold Singularity image."""

import os
import signal
from pathlib import Path
from typing import Tuple
import tempfile

from absl import app
from absl import flags
from absl import logging
from spython.main import Client


#### USER CONFIGURATION ####

# Set to target of scripts/download_all_databases.sh
DOWNLOAD_DIR = 'SET ME'

# AlphaFold Singularity image.
singularity_image = Client.load('alphafold.sif')

# Path to a directory that will store the results.
if 'TMPDIR' in os.environ:
  output_dir = os.environ['TMPDIR']
else:
  output_dir = tempfile.mkdtemp(dir='/tmp', prefix='alphafold-')

# Names of models to use.
model_names = [
    'model_1',
    'model_2',
    'model_3',
    'model_4',
    'model_5',
]

# You can individually override the following paths if you have placed the
# data in locations other than the DOWNLOAD_DIR.

# Path to directory of supporting data, contains 'params' dir.
data_dir = DOWNLOAD_DIR

# Path to the Uniref90 database for use by JackHMMER.
uniref90_database_path = os.path.join(
    DOWNLOAD_DIR, 'uniref90', 'uniref90.fasta')

# Path to the MGnify database for use by JackHMMER.
mgnify_database_path = os.path.join(
    DOWNLOAD_DIR, 'mgnify', 'mgy_clusters_2018_12.fa')

# Path to the BFD database for use by HHblits.
bfd_database_path = os.path.join(
    DOWNLOAD_DIR, 'bfd', 'bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt')

# Path to the Small BFD database for use by JackHMMER.
small_bfd_database_path = os.path.join(
    DOWNLOAD_DIR, 'small_bfd', 'bfd-first_non_consensus_sequences.fasta')

# Path to the Uniclust30 database for use by HHblits.
uniclust30_database_path = os.path.join(
    DOWNLOAD_DIR, 'uniclust30', 'uniclust30_2018_08', 'uniclust30_2018_08')

# Path to the PDB70 database for use by HHsearch.
pdb70_database_path = os.path.join(DOWNLOAD_DIR, 'pdb70', 'pdb70')

# Path to a directory with template mmCIF structures, each named <pdb_id>.cif')
template_mmcif_dir = os.path.join(DOWNLOAD_DIR, 'pdb_mmcif', 'mmcif_files')

# Path to a file mapping obsolete PDB IDs to their replacements.
obsolete_pdbs_path = os.path.join(DOWNLOAD_DIR, 'pdb_mmcif', 'obsolete.dat')

#### END OF USER CONFIGURATION ####


flags.DEFINE_bool('use_gpu', True, 'Enable NVIDIA runtime to run with GPUs.')
flags.DEFINE_string('gpu_devices', 'all', 'Comma separated list of devices to '
                    'pass to NVIDIA_VISIBLE_DEVICES.')
flags.DEFINE_list('fasta_paths', None, 'Paths to FASTA files, each containing '
                  'one sequence. Paths should be separated by commas. '
                  'All FASTA paths must have a unique basename as the '
                  'basename is used to name the output directories for '
                  'each prediction.')
flags.DEFINE_string('max_template_date', None, 'Maximum template release date '
                    'to consider (ISO-8601 format - i.e. YYYY-MM-DD). '
                    'Important if folding historical test sets.')
flags.DEFINE_enum('preset', 'full_dbs',
                  ['reduced_dbs', 'full_dbs', 'casp14'],
                  'Choose preset model configuration - no ensembling and '
                  'smaller genetic database config (reduced_dbs), no '
                  'ensembling and full genetic database config  (full_dbs) or '
                  'full genetic database config and 8 model ensemblings '
                  '(casp14).')
flags.DEFINE_boolean('benchmark', False, 'Run multiple JAX model evaluations '
                     'to obtain a timing that excludes the compilation time, '
                     'which should be more indicative of the time required for '
                     'inferencing many proteins.')

FLAGS = flags.FLAGS

_ROOT_MOUNT_DIRECTORY = '/mnt/'


def _create_bind(bind_name: str, path: str) -> Tuple[str, str]:
  path = os.path.abspath(path)
  source_path = os.path.dirname(path) if bind_name != 'data_dir' else path
  target_path = os.path.join(_ROOT_MOUNT_DIRECTORY, bind_name)

  logging.info('Binding %s -> %s', source_path, target_path)

  # NOTE singularity binds are read-only by default
  if bind_name == 'data_dir':
    data_path = target_path
  else:
    data_path = f'{os.path.join(target_path, os.path.basename(path))}'

  return (f'{source_path}:{target_path}', f'{data_path}')


def main(argv):
  if len(argv) > 1:
    raise app.UsageError('Too many command-line arguments.')

  binds = []
  command_args = []

  # Mount each fasta path as a unique target directory.
  target_fasta_paths = []
  for i, fasta_path in enumerate(FLAGS.fasta_paths):
    bind, target_path = _create_bind(f'fasta_path_{i}', fasta_path)
    binds.append(bind)
    target_fasta_paths.append(target_path)
  command_args.append(f'--fasta_paths={",".join(target_fasta_paths)}')

  database_paths = [
      ('uniref90_database_path', uniref90_database_path),
      ('mgnify_database_path', mgnify_database_path),
      ('pdb70_database_path', pdb70_database_path),
      ('data_dir', data_dir),
      ('template_mmcif_dir', template_mmcif_dir),
      ('obsolete_pdbs_path', obsolete_pdbs_path),
  ]

  if FLAGS.preset == 'reduced_dbs':
    database_paths.append(('small_bfd_database_path', small_bfd_database_path))
  else:
    database_paths.extend([
      ('uniclust30_database_path', uniclust30_database_path),
      ('bfd_database_path', bfd_database_path),
    ])

  for name, path in database_paths:
    if path:
      bind, target_path = _create_bind(name, path)
      binds.append(bind)
      command_args.append(f'--{name}={target_path}')

  output_target_path = os.path.join(_ROOT_MOUNT_DIRECTORY, 'output')
  binds.append(f'{output_dir}:{output_target_path}')

  command_args.extend([
    f'--output_dir={output_target_path}',
    f'--model_names={",".join(model_names)}',
    f'--max_template_date={FLAGS.max_template_date}',
    f'--preset={FLAGS.preset}',
    f'--benchmark={FLAGS.benchmark}',
    '--logtostderr',
  ])

  options = [
    '--bind', f'{",".join(binds)}',
    '--env', 'TF_FORCE_UNIFIED_MEMORY=1',
    '--env', 'XLA_PYTHON_CLIENT_MEM_FRACTION=4.0'
  ]

  if FLAGS.use_gpu and FLAGS.gpu_devices:
    options.extend(['--env', f'NVIDIA_VISIBLE_DEVICES={FLAGS.gpu_devices}'])

  # result is a dict with keys "message" (value = all output as a single string), 
  # and "return_code" (value = integer return code)
  result = Client.run(
    singularity_image,
    command_args,
    nv=True if FLAGS.use_gpu else None,
    return_result=True,
    options=options
  )

  for line in result["message"].strip().split('\n'):
    print(f'{line}')

  return result['return_code']


if __name__ == '__main__':
  flags.mark_flags_as_required([
      'fasta_paths',
      'max_template_date',
  ])
  app.run(main)
