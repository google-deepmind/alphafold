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

"""Full AlphaFold protein structure prediction script."""
import enum
import json
import os
import pathlib
import pickle
import random
import shutil
import sys
import time
from typing import Any, Dict, Union
from tqdm import tqdm

from absl import app
from absl import flags
from absl import logging
from alphafold.common import confidence
from alphafold.common import protein
from alphafold.data import pipeline
from alphafold.data import pipeline_multimer
from alphafold.common import protein
from alphafold.relax import relax

from alphafold.model import model
from alphafold.model import config
from alphafold.model import data

import jax.numpy as jnp
import numpy as np

# Internal import (7716).

logging.set_verbosity(logging.INFO)


@enum.unique
class ModelsToRelax(enum.Enum):
  ALL = 0
  BEST = 1
  NONE = 2
@enum.unique
class ModelType(enum.Enum):
  MONOMER = 0
  MULTIMER = 1

flags.DEFINE_string('precomputed_msa', None, 'MSA to use for this run')

flags.DEFINE_string('data_dir', None, 'Path to directory of supporting data.')
flags.DEFINE_string('output_dir', None, 'Path to a directory that will '
                    'store the results.')
flags.DEFINE_integer('random_seed', None, 'The random seed for the data '
                     'pipeline. By default, this is randomly generated. Note '
                     'that even if this is set, Alphafold may still not be '
                     'deterministic, because processes like GPU inference are '
                     'nondeterministic.')
flags.DEFINE_enum_class('models_to_relax', ModelsToRelax.BEST, ModelsToRelax,
                        'The models to run the final relaxation step on. '
                        'If `all`, all models are relaxed, which may be time '
                        'consuming. If `best`, only the most confident model '
                        'is relaxed. If `none`, relaxation is not run. Turning '
                        'off relaxation might result in predictions with '
                        'distracting stereochemical violations but might help '
                        'in case you are having issues with the relaxation '
                        'stage.')
flags.DEFINE_boolean('use_gpu_relax', None, 'Whether to relax on GPU. '
                     'Relax on GPU can be much faster than CPU, so it is '
                     'recommended to enable if possible. GPUs must be available'
                     ' if this setting is enabled.')

flags.DEFINE_integer('num_multimer_predictions_per_model', 5, 'How many '
                     'predictions (each with a different random seed) will be '
                     'generated per model. E.g. if this is 2 and there are 5 '
                     'models then there will be 10 predictions per input. '
                     'Note: this FLAG only applies if model_preset=multimer')

FLAGS = flags.FLAGS

MAX_TEMPLATE_HITS = 20
RELAX_MAX_ITERATIONS = 0
RELAX_ENERGY_TOLERANCE = 2.39
RELAX_STIFFNESS = 10.0
RELAX_EXCLUDE_RESIDUES = []
RELAX_MAX_OUTER_ITERATIONS = 3


def _check_flag(flag_name: str,
                other_flag_name: str,
                should_be_set: bool):
  if should_be_set != bool(FLAGS[flag_name].value):
    verb = 'be' if should_be_set else 'not be'
    raise ValueError(f'{flag_name} must {verb} set when running with '
                     f'"--{other_flag_name}={FLAGS[other_flag_name].value}".')


def _jnp_to_np(output: Dict[str, Any]) -> Dict[str, Any]:
  """Recursively changes jax arrays to numpy arrays."""
  for k, v in output.items():
    if isinstance(v, dict):
      output[k] = _jnp_to_np(v)
    elif isinstance(v, jnp.ndarray):
      output[k] = np.array(v)
  return output


def _save_confidence_json_file(
    plddt: np.ndarray, output_dir: str, model_name: str
) -> None:
  confidence_json = confidence.confidence_json(plddt)

  # Save the confidence json.
  confidence_json_output_path = os.path.join(
      output_dir, f'confidence_{model_name}.json'
  )
  with open(confidence_json_output_path, 'w') as f:
    f.write(confidence_json)


def _save_pae_json_file(
    pae: np.ndarray, max_pae: float, output_dir: str, model_name: str
) -> None:
  """Check prediction result for PAE data and save to a JSON file if present.

  Args:
    pae: The n_res x n_res PAE array.
    max_pae: The maximum possible PAE value.
    output_dir: Directory to which files are saved.
    model_name: Name of a model.
  """
  pae_json = confidence.pae_json(pae, max_pae)

  # Save the PAE json.
  pae_json_output_path = os.path.join(output_dir, f'pae_{model_name}.json')
  with open(pae_json_output_path, 'w') as f:
    f.write(pae_json)


def predict_structure(
    precomputed_msa: str,
    output_dir_base : str,
    data_pipeline: Union[pipeline.DataPipeline, pipeline_multimer.DataPipeline],
    model_type_to_use = ModelType.MONOMER,
    multimer_model_max_num_recycles = 3,
    run_relax = False,
    relax_use_gpu  = False
    ):
  """Predicts structure using AlphaFold for the given sequence."""
  logging.info('Predicting %s', precomputed_msa)
  timings = {}
  msa_name = os.path.basename(precomputed_msa)
  output_dir = os.path.join(output_dir_base, msa_name)
  if not os.path.exists(output_dir):
    os.makedirs(output_dir)
  
    # Get features.
  t_0 = time.time()
  feature_dict = data_pipeline.process(precomputed_msa = precomputed_msa, output_dir=output_dir)
  timings['features'] = time.time() - t_0
  # Write out features as a pickled dictionary.
  features_output_path = os.path.join(output_dir, 'features.pkl')
  with open(features_output_path, 'wb') as f:
    pickle.dump(feature_dict, f, protocol=4)

  features_for_chain = {}
  features_for_chain[protein.PDB_CHAIN_IDS[0]] = feature_dict

# Do further feature post-processing depending on the model type.  
  np_example = features_for_chain[protein.PDB_CHAIN_IDS[0]]  
  model_names = config.MODEL_PRESETS['monomer'] + ('model_2_ptm',)

  plddts = {}
  ranking_confidences = {}
  pae_outputs = {}
  unrelaxed_proteins = {}
  TQDM_BAR_FORMAT = '{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]'
  with tqdm(total=len(model_names) + 1, bar_format=TQDM_BAR_FORMAT) as pbar:
    for model_name in model_names:
      pbar.set_description(f'Running {model_name}')

      cfg = config.model_config(model_name)

      if model_type_to_use == ModelType.MONOMER:
        cfg.data.eval.num_ensemble = 1
      elif model_type_to_use == ModelType.MULTIMER:
        cfg.model.num_ensemble_eval = 1

      if model_type_to_use == ModelType.MULTIMER:
        cfg.model.num_recycle = multimer_model_max_num_recycles
        cfg.model.recycle_early_stop_tolerance = 0.5

      params = data.get_model_haiku_params(model_name, './alphafold/data')
      model_runner = model.RunModel(cfg, params)
      processed_feature_dict = model_runner.process_features(np_example, random_seed=0)
      prediction = model_runner.predict(processed_feature_dict, random_seed=random.randrange(sys.maxsize))

      mean_plddt = prediction['plddt'].mean()

      if model_type_to_use == ModelType.MONOMER:
        if 'predicted_aligned_error' in prediction:
          pae_outputs[model_name] = (prediction['predicted_aligned_error'],
                                    prediction['max_predicted_aligned_error'])
        else:
          # Monomer models are sorted by mean pLDDT. Do not put monomer pTM models here as they
          # should never get selected.
          ranking_confidences[model_name] = prediction['ranking_confidence']
          plddts[model_name] = prediction['plddt']
      elif model_type_to_use == ModelType.MULTIMER:
        # Multimer models are sorted by pTM+ipTM.
        ranking_confidences[model_name] = prediction['ranking_confidence']
        plddts[model_name] = prediction['plddt']
        pae_outputs[model_name] = (prediction['predicted_aligned_error'],
                                  prediction['max_predicted_aligned_error'])

      # Set the b-factors to the per-residue plddt.
      final_atom_mask = prediction['structure_module']['final_atom_mask']
      b_factors = prediction['plddt'][:, None] * final_atom_mask
      unrelaxed_protein = protein.from_prediction(
          processed_feature_dict,
          prediction,
          b_factors=b_factors,
          remove_leading_feature_dimension=(
              model_type_to_use == ModelType.MONOMER))
      unrelaxed_proteins[model_name] = unrelaxed_protein

      # Delete unused outputs to save memory.
      del model_runner
      del params
      del prediction
      pbar.update(n=1)

    # --- AMBER relax the best model ---

    # Find the best model according to the mean pLDDT.
    best_model_name = max(ranking_confidences.keys(), key=lambda x: ranking_confidences[x])

    if run_relax:
      pbar.set_description(f'AMBER relaxation')
      amber_relaxer = relax.AmberRelaxation(
          max_iterations=0,
          tolerance=2.39,
          stiffness=10.0,
          exclude_residues=[],
          max_outer_iterations=3,
          use_gpu=relax_use_gpu)
      relaxed_pdb, _, _ = amber_relaxer.process(prot=unrelaxed_proteins[best_model_name])
    else:
      print('Warning: Running without the relaxation stage.')
      relaxed_pdb = protein.to_pdb(unrelaxed_proteins[best_model_name])
    pbar.update(n=1)  # Finished AMBER relax.

  # Write out the prediction
  pred_output_path = os.path.join(output_dir, 'selected_prediction.pdb')
  with open(pred_output_path, 'w') as f:
    f.write(relaxed_pdb)



def main(argv):
  if len(argv) > 1:
    raise app.UsageError('Too many command-line arguments.')
  
  monomer_data_pipeline = pipeline.DataPipeline()

  run_multimer_system = False
  if run_multimer_system:
    num_predictions_per_model = FLAGS.num_multimer_predictions_per_model
    data_pipeline = pipeline_multimer.DataPipeline(
        monomer_data_pipeline=monomer_data_pipeline,
        )
  else:
    num_predictions_per_model = 1
    data_pipeline = monomer_data_pipeline


  amber_relaxer = relax.AmberRelaxation(
      max_iterations=RELAX_MAX_ITERATIONS,
      tolerance=RELAX_ENERGY_TOLERANCE,
      stiffness=RELAX_STIFFNESS,
      exclude_residues=RELAX_EXCLUDE_RESIDUES,
      max_outer_iterations=RELAX_MAX_OUTER_ITERATIONS,
      use_gpu=FLAGS.use_gpu_relax)

  random_seed = FLAGS.random_seed
  if random_seed is None:
    random_seed = random.randrange(sys.maxsize // num_predictions_per_model)
  logging.info('Using random seed %d for the data pipeline', random_seed)

  # Predict structure for each of the sequences.
  predict_structure(
      precomputed_msa=FLAGS.precomputed_msa,
      output_dir_base=FLAGS.output_dir,
      data_pipeline=data_pipeline,
)


if __name__ == '__main__':
  flags.mark_flags_as_required([
      'precomputed_msa',
      'output_dir',
      'data_dir',
  ])

  app.run(main)
