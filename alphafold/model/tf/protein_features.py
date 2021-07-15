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

"""Contains descriptions of various protein features."""
import enum
from typing import Dict, Optional, Sequence, Tuple, Union

import tensorflow.compat.v1 as tf

from alphafold.common import residue_constants

# Type aliases.
FeaturesMetadata = Dict[str, Tuple[tf.dtypes.DType, Sequence[Union[str, int]]]]


class FeatureType(enum.Enum):
  ZERO_DIM = 0  # Shape [x]
  ONE_DIM = 1  # Shape [num_res, x]
  TWO_DIM = 2  # Shape [num_res, num_res, x]
  MSA = 3  # Shape [msa_length, num_res, x]


# Placeholder values that will be replaced with their true value at runtime.
NUM_RES = "num residues placeholder"
NUM_SEQ = "length msa placeholder"
NUM_TEMPLATES = "num templates placeholder"
# Sizes of the protein features, NUM_RES and NUM_SEQ are allowed as placeholders
# to be replaced with the number of residues and the number of sequences in the
# multiple sequence alignment, respectively.


FEATURES = {
    #### Static features of a protein sequence ####
    "aatype": (tf.float32, [NUM_RES, 21]),
    "between_segment_residues": (tf.int64, [NUM_RES, 1]),
    "deletion_matrix": (tf.float32, [NUM_SEQ, NUM_RES, 1]),
    "domain_name": (tf.string, [1]),
    "msa": (tf.int64, [NUM_SEQ, NUM_RES, 1]),
    "num_alignments": (tf.int64, [NUM_RES, 1]),
    "residue_index": (tf.int64, [NUM_RES, 1]),
    "seq_length": (tf.int64, [NUM_RES, 1]),
    "sequence": (tf.string, [1]),
    "all_atom_positions": (tf.float32,
                           [NUM_RES, residue_constants.atom_type_num, 3]),
    "all_atom_mask": (tf.int64, [NUM_RES, residue_constants.atom_type_num]),
    "resolution": (tf.float32, [1]),
    "template_domain_names": (tf.string, [NUM_TEMPLATES]),
    "template_sum_probs": (tf.float32, [NUM_TEMPLATES, 1]),
    "template_aatype": (tf.float32, [NUM_TEMPLATES, NUM_RES, 22]),
    "template_all_atom_positions": (tf.float32, [
        NUM_TEMPLATES, NUM_RES, residue_constants.atom_type_num, 3
    ]),
    "template_all_atom_masks": (tf.float32, [
        NUM_TEMPLATES, NUM_RES, residue_constants.atom_type_num, 1
    ]),
}

FEATURE_TYPES = {k: v[0] for k, v in FEATURES.items()}
FEATURE_SIZES = {k: v[1] for k, v in FEATURES.items()}


def register_feature(name: str,
                     type_: tf.dtypes.DType,
                     shape_: Tuple[Union[str, int]]):
  """Register extra features used in custom datasets."""
  FEATURES[name] = (type_, shape_)
  FEATURE_TYPES[name] = type_
  FEATURE_SIZES[name] = shape_


def shape(feature_name: str,
          num_residues: int,
          msa_length: int,
          num_templates: Optional[int] = None,
          features: Optional[FeaturesMetadata] = None):
  """Get the shape for the given feature name.

  This is near identical to _get_tf_shape_no_placeholders() but with 2
  differences:
  * This method does not calculate a single placeholder from the total number of
    elements (eg given <NUM_RES, 3> and size := 12, this won't deduce NUM_RES
    must be 4)
  * This method will work with tensors

  Args:
    feature_name: String identifier for the feature. If the feature name ends
      with "_unnormalized", theis suffix is stripped off.
    num_residues: The number of residues in the current domain - some elements
      of the shape can be dynamic and will be replaced by this value.
    msa_length: The number of sequences in the multiple sequence alignment, some
      elements of the shape can be dynamic and will be replaced by this value.
      If the number of alignments is unknown / not read, please pass None for
      msa_length.
    num_templates (optional): The number of templates in this tfexample.
    features: A feature_name to (tf_dtype, shape) lookup; defaults to FEATURES.

  Returns:
    List of ints representation the tensor size.

  Raises:
    ValueError: If a feature is requested but no concrete placeholder value is
        given.
  """
  features = features or FEATURES
  if feature_name.endswith("_unnormalized"):
    feature_name = feature_name[:-13]

  unused_dtype, raw_sizes = features[feature_name]
  replacements = {NUM_RES: num_residues,
                  NUM_SEQ: msa_length}

  if num_templates is not None:
    replacements[NUM_TEMPLATES] = num_templates

  sizes = [replacements.get(dimension, dimension) for dimension in raw_sizes]
  for dimension in sizes:
    if isinstance(dimension, str):
      raise ValueError("Could not parse %s (shape: %s) with values: %s" % (
          feature_name, raw_sizes, replacements))
  return sizes

