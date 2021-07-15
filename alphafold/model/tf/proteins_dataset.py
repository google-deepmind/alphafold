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

"""Datasets consisting of proteins."""
from typing import Dict, Mapping, Optional, Sequence

import numpy as np
import tensorflow.compat.v1 as tf

from alphafold.model.tf import protein_features

TensorDict = Dict[str, tf.Tensor]


def parse_tfexample(
    raw_data: bytes,
    features: protein_features.FeaturesMetadata,
    key: Optional[str] = None) -> Dict[str, tf.train.Feature]:
  """Read a single TF Example proto and return a subset of its features.

  Args:
    raw_data: A serialized tf.Example proto.
    features: A dictionary of features, mapping string feature names to a tuple
      (dtype, shape). This dictionary should be a subset of
      protein_features.FEATURES (or the dictionary itself for all features).
    key: Optional string with the SSTable key of that tf.Example. This will be
      added into features as a 'key' but only if requested in features.

  Returns:
    A dictionary of features mapping feature names to features. Only the given
    features are returned, all other ones are filtered out.
  """
  feature_map = {
      k: tf.io.FixedLenSequenceFeature(shape=(), dtype=v[0], allow_missing=True)
      for k, v in features.items()
  }
  parsed_features = tf.io.parse_single_example(raw_data, feature_map)
  reshaped_features = parse_reshape_logic(parsed_features, features, key=key)

  return reshaped_features


def _first(tensor: tf.Tensor) -> tf.Tensor:
  """Returns the 1st element - the input can be a tensor or a scalar."""
  return tf.reshape(tensor, shape=(-1,))[0]


def parse_reshape_logic(
    parsed_features: TensorDict,
    features: protein_features.FeaturesMetadata,
    key: Optional[str] = None) -> TensorDict:
  """Transforms parsed serial features to the correct shape."""
  # Find out what is the number of sequences and the number of alignments.
  num_residues = tf.cast(_first(parsed_features["seq_length"]), dtype=tf.int32)

  if "num_alignments" in parsed_features:
    num_msa = tf.cast(_first(parsed_features["num_alignments"]), dtype=tf.int32)
  else:
    num_msa = 0

  if "template_domain_names" in parsed_features:
    num_templates = tf.cast(
        tf.shape(parsed_features["template_domain_names"])[0], dtype=tf.int32)
  else:
    num_templates = 0

  if key is not None and "key" in features:
    parsed_features["key"] = [key]  # Expand dims from () to (1,).

  # Reshape the tensors according to the sequence length and num alignments.
  for k, v in parsed_features.items():
    new_shape = protein_features.shape(
        feature_name=k,
        num_residues=num_residues,
        msa_length=num_msa,
        num_templates=num_templates,
        features=features)
    new_shape_size = tf.constant(1, dtype=tf.int32)
    for dim in new_shape:
      new_shape_size *= tf.cast(dim, tf.int32)

    assert_equal = tf.assert_equal(
        tf.size(v), new_shape_size,
        name="assert_%s_shape_correct" % k,
        message="The size of feature %s (%s) could not be reshaped "
        "into %s" % (k, tf.size(v), new_shape))
    if "template" not in k:
      # Make sure the feature we are reshaping is not empty.
      assert_non_empty = tf.assert_greater(
          tf.size(v), 0, name="assert_%s_non_empty" % k,
          message="The feature %s is not set in the tf.Example. Either do not "
          "request the feature or use a tf.Example that has the "
          "feature set." % k)
      with tf.control_dependencies([assert_non_empty, assert_equal]):
        parsed_features[k] = tf.reshape(v, new_shape, name="reshape_%s" % k)
    else:
      with tf.control_dependencies([assert_equal]):
        parsed_features[k] = tf.reshape(v, new_shape, name="reshape_%s" % k)

  return parsed_features


def _make_features_metadata(
    feature_names: Sequence[str]) -> protein_features.FeaturesMetadata:
  """Makes a feature name to type and shape mapping from a list of names."""
  # Make sure these features are always read.
  required_features = ["aatype", "sequence", "seq_length"]
  feature_names = list(set(feature_names) | set(required_features))

  features_metadata = {name: protein_features.FEATURES[name]
                       for name in feature_names}
  return features_metadata


def create_tensor_dict(
    raw_data: bytes,
    features: Sequence[str],
    key: Optional[str] = None,
    ) -> TensorDict:
  """Creates a dictionary of tensor features.

  Args:
    raw_data: A serialized tf.Example proto.
    features: A list of strings of feature names to be returned in the dataset.
    key: Optional string with the SSTable key of that tf.Example. This will be
      added into features as a 'key' but only if requested in features.

  Returns:
    A dictionary of features mapping feature names to features. Only the given
    features are returned, all other ones are filtered out.
  """
  features_metadata = _make_features_metadata(features)
  return parse_tfexample(raw_data, features_metadata, key)


def np_to_tensor_dict(
    np_example: Mapping[str, np.ndarray],
    features: Sequence[str],
    ) -> TensorDict:
  """Creates dict of tensors from a dict of NumPy arrays.

  Args:
    np_example: A dict of NumPy feature arrays.
    features: A list of strings of feature names to be returned in the dataset.

  Returns:
    A dictionary of features mapping feature names to features. Only the given
    features are returned, all other ones are filtered out.
  """
  features_metadata = _make_features_metadata(features)
  tensor_dict = {k: tf.constant(v) for k, v in np_example.items()
                 if k in features_metadata}

  # Ensures shapes are as expected. Needed for setting size of empty features
  # e.g. when no template hits were found.
  tensor_dict = parse_reshape_logic(tensor_dict, features_metadata)
  return tensor_dict
