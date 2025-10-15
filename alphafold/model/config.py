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
"""Model config."""

import copy
import functools
from typing import Any, Dict, Final, Optional, Sequence

from alphafold.model import base_config
from alphafold.model.tf import shape_placeholders
import ml_collections


NUM_RES: Final = shape_placeholders.NUM_RES
NUM_MSA_SEQ: Final = shape_placeholders.NUM_MSA_SEQ
NUM_EXTRA_SEQ: Final = shape_placeholders.NUM_EXTRA_SEQ
NUM_TEMPLATES: Final = shape_placeholders.NUM_TEMPLATES

MODEL_PRESETS = {
    'monomer': (
        'model_1',
        'model_2',
        'model_3',
        'model_4',
        'model_5',
    ),
    'monomer_ptm': (
        'model_1_ptm',
        'model_2_ptm',
        'model_3_ptm',
        'model_4_ptm',
        'model_5_ptm',
    ),
    'multimer': (
        'model_1_multimer_v3',
        'model_2_multimer_v3',
        'model_3_multimer_v3',
        'model_4_multimer_v3',
        'model_5_multimer_v3',
    ),
}
MODEL_PRESETS['monomer_casp14'] = MODEL_PRESETS['monomer']

CONFIG_DIFFS = {
    'model_1': {
        # Jumper et al. (2021) Suppl. Table 5, Model 1.1.1
        'data.common.max_extra_msa': 5120,
        'data.common.reduce_msa_clusters_by_max_templates': True,
        'data.common.use_templates': True,
        'model.embeddings_and_evoformer.template.embed_torsion_angles': True,
        'model.embeddings_and_evoformer.template.enabled': True,
    },
    'model_2': {
        # Jumper et al. (2021) Suppl. Table 5, Model 1.1.2
        'data.common.reduce_msa_clusters_by_max_templates': True,
        'data.common.use_templates': True,
        'model.embeddings_and_evoformer.template.embed_torsion_angles': True,
        'model.embeddings_and_evoformer.template.enabled': True,
    },
    'model_3': {
        # Jumper et al. (2021) Suppl. Table 5, Model 1.2.1
        'data.common.max_extra_msa': 5120,
    },
    'model_4': {
        # Jumper et al. (2021) Suppl. Table 5, Model 1.2.2
        'data.common.max_extra_msa': 5120,
    },
    'model_5': {
        # Jumper et al. (2021) Suppl. Table 5, Model 1.2.3
    },
    # The following models are fine-tuned from the corresponding models above
    # with an additional predicted_aligned_error head that can produce
    # predicted TM-score (pTM) and predicted aligned errors.
    'model_1_ptm': {
        'data.common.max_extra_msa': 5120,
        'data.common.reduce_msa_clusters_by_max_templates': True,
        'data.common.use_templates': True,
        'model.embeddings_and_evoformer.template.embed_torsion_angles': True,
        'model.embeddings_and_evoformer.template.enabled': True,
        'model.heads.predicted_aligned_error.weight': 0.1,
    },
    'model_2_ptm': {
        'data.common.reduce_msa_clusters_by_max_templates': True,
        'data.common.use_templates': True,
        'model.embeddings_and_evoformer.template.embed_torsion_angles': True,
        'model.embeddings_and_evoformer.template.enabled': True,
        'model.heads.predicted_aligned_error.weight': 0.1,
    },
    'model_3_ptm': {
        'data.common.max_extra_msa': 5120,
        'model.heads.predicted_aligned_error.weight': 0.1,
    },
    'model_4_ptm': {
        'data.common.max_extra_msa': 5120,
        'model.heads.predicted_aligned_error.weight': 0.1,
    },
    'model_5_ptm': {'model.heads.predicted_aligned_error.weight': 0.1},
    'model_1_multimer_v3': {},
    'model_2_multimer_v3': {},
    'model_3_multimer_v3': {},
    'model_4_multimer_v3': {
        'model.embeddings_and_evoformer.num_extra_msa': 1152
    },
    'model_5_multimer_v3': {
        'model.embeddings_and_evoformer.num_extra_msa': 1152
    },
}
# Key differences between multimer v1/v2 and v3, mostly due to numerical
# optimisations in the TriangleMultiplication module.
common_updates = {
    'model.embeddings_and_evoformer.num_msa': 252,
    'model.embeddings_and_evoformer.num_extra_msa': 1152,
    'model.embeddings_and_evoformer.evoformer.triangle_multiplication_incoming.fuse_projection_weights': False,  # pylint: disable=line-too-long
    'model.embeddings_and_evoformer.evoformer.triangle_multiplication_outgoing.fuse_projection_weights': False,  # pylint: disable=line-too-long
    'model.embeddings_and_evoformer.template.template_pair_stack.triangle_multiplication_incoming.fuse_projection_weights': False,  # pylint: disable=line-too-long
    'model.embeddings_and_evoformer.template.template_pair_stack.triangle_multiplication_outgoing.fuse_projection_weights': False,  # pylint: disable=line-too-long
}
CONFIG_DIFFS.update(
    {f'model_{i}_multimer': common_updates for i in range(1, 6)}
)
CONFIG_DIFFS.update(
    {f'model_{i}_multimer_v2': common_updates for i in range(1, 6)}
)

CONFIG = ml_collections.ConfigDict({
    'data': {
        'common': {
            'masked_msa': {
                'profile_prob': 0.1,
                'same_prob': 0.1,
                'uniform_prob': 0.1,
            },
            'max_extra_msa': 1024,
            'msa_cluster_features': True,
            'num_recycle': 3,
            'reduce_msa_clusters_by_max_templates': False,
            'resample_msa_in_recycling': True,
            'template_features': [
                'template_all_atom_positions',
                'template_sum_probs',
                'template_aatype',
                'template_all_atom_masks',
                'template_domain_names',
            ],
            'unsupervised_features': [
                'aatype',
                'residue_index',
                'sequence',
                'msa',
                'domain_name',
                'num_alignments',
                'seq_length',
                'between_segment_residues',
                'deletion_matrix',
            ],
            'use_templates': False,
        },
        'eval': {
            'feat': {
                'aatype': [NUM_RES],
                'all_atom_mask': [NUM_RES, None],
                'all_atom_positions': [NUM_RES, None, None],
                'alt_chi_angles': [NUM_RES, None],
                'atom14_alt_gt_exists': [NUM_RES, None],
                'atom14_alt_gt_positions': [NUM_RES, None, None],
                'atom14_atom_exists': [NUM_RES, None],
                'atom14_atom_is_ambiguous': [NUM_RES, None],
                'atom14_gt_exists': [NUM_RES, None],
                'atom14_gt_positions': [NUM_RES, None, None],
                'atom37_atom_exists': [NUM_RES, None],
                'backbone_affine_mask': [NUM_RES],
                'backbone_affine_tensor': [NUM_RES, None],
                'bert_mask': [NUM_MSA_SEQ, NUM_RES],
                'chi_angles': [NUM_RES, None],
                'chi_mask': [NUM_RES, None],
                'extra_deletion_value': [NUM_EXTRA_SEQ, NUM_RES],
                'extra_has_deletion': [NUM_EXTRA_SEQ, NUM_RES],
                'extra_msa': [NUM_EXTRA_SEQ, NUM_RES],
                'extra_msa_mask': [NUM_EXTRA_SEQ, NUM_RES],
                'extra_msa_row_mask': [NUM_EXTRA_SEQ],
                'is_distillation': [],
                'msa_feat': [NUM_MSA_SEQ, NUM_RES, None],
                'msa_mask': [NUM_MSA_SEQ, NUM_RES],
                'msa_row_mask': [NUM_MSA_SEQ],
                'pseudo_beta': [NUM_RES, None],
                'pseudo_beta_mask': [NUM_RES],
                'random_crop_to_size_seed': [None],
                'residue_index': [NUM_RES],
                'residx_atom14_to_atom37': [NUM_RES, None],
                'residx_atom37_to_atom14': [NUM_RES, None],
                'resolution': [],
                'rigidgroups_alt_gt_frames': [NUM_RES, None, None],
                'rigidgroups_group_exists': [NUM_RES, None],
                'rigidgroups_group_is_ambiguous': [NUM_RES, None],
                'rigidgroups_gt_exists': [NUM_RES, None],
                'rigidgroups_gt_frames': [NUM_RES, None, None],
                'seq_length': [],
                'seq_mask': [NUM_RES],
                'target_feat': [NUM_RES, None],
                'template_aatype': [NUM_TEMPLATES, NUM_RES],
                'template_all_atom_masks': [NUM_TEMPLATES, NUM_RES, None],
                'template_all_atom_positions': [
                    NUM_TEMPLATES,
                    NUM_RES,
                    None,
                    None,
                ],
                'template_backbone_affine_mask': [NUM_TEMPLATES, NUM_RES],
                'template_backbone_affine_tensor': [
                    NUM_TEMPLATES,
                    NUM_RES,
                    None,
                ],
                'template_mask': [NUM_TEMPLATES],
                'template_pseudo_beta': [NUM_TEMPLATES, NUM_RES, None],
                'template_pseudo_beta_mask': [NUM_TEMPLATES, NUM_RES],
                'template_sum_probs': [NUM_TEMPLATES, None],
                'true_msa': [NUM_MSA_SEQ, NUM_RES],
            },
            'fixed_size': True,
            'subsample_templates': False,  # We want top templates.
            'masked_msa_replace_fraction': 0.15,
            'max_msa_clusters': 512,
            'max_templates': 4,
            'num_ensemble': 1,
        },
    },
    'model': {
        'embeddings_and_evoformer': {
            'evoformer_num_block': 48,
            'evoformer': {
                'msa_row_attention_with_pair_bias': {
                    'dropout_rate': 0.15,
                    'gating': True,
                    'num_head': 8,
                    'orientation': 'per_row',
                    'shared_dropout': True,
                },
                'msa_column_attention': {
                    'dropout_rate': 0.0,
                    'gating': True,
                    'num_head': 8,
                    'orientation': 'per_column',
                    'shared_dropout': True,
                },
                'msa_transition': {
                    'dropout_rate': 0.0,
                    'num_intermediate_factor': 4,
                    'orientation': 'per_row',
                    'shared_dropout': True,
                },
                'outer_product_mean': {
                    'first': False,
                    'chunk_size': 128,
                    'dropout_rate': 0.0,
                    'num_outer_channel': 32,
                    'orientation': 'per_row',
                    'shared_dropout': True,
                },
                'triangle_attention_starting_node': {
                    'dropout_rate': 0.25,
                    'gating': True,
                    'num_head': 4,
                    'orientation': 'per_row',
                    'shared_dropout': True,
                },
                'triangle_attention_ending_node': {
                    'dropout_rate': 0.25,
                    'gating': True,
                    'num_head': 4,
                    'orientation': 'per_column',
                    'shared_dropout': True,
                },
                'triangle_multiplication_outgoing': {
                    'dropout_rate': 0.25,
                    'equation': 'ikc,jkc->ijc',
                    'num_intermediate_channel': 128,
                    'orientation': 'per_row',
                    'shared_dropout': True,
                    'fuse_projection_weights': False,
                },
                'triangle_multiplication_incoming': {
                    'dropout_rate': 0.25,
                    'equation': 'kjc,kic->ijc',
                    'num_intermediate_channel': 128,
                    'orientation': 'per_row',
                    'shared_dropout': True,
                    'fuse_projection_weights': False,
                },
                'pair_transition': {
                    'dropout_rate': 0.0,
                    'num_intermediate_factor': 4,
                    'orientation': 'per_row',
                    'shared_dropout': True,
                },
            },
            'extra_msa_channel': 64,
            'extra_msa_stack_num_block': 4,
            'max_relative_feature': 32,
            'msa_channel': 256,
            'pair_channel': 128,
            'prev_pos': {'min_bin': 3.25, 'max_bin': 20.75, 'num_bins': 15},
            'recycle_features': True,
            'recycle_pos': True,
            'seq_channel': 384,
            'template': {
                'attention': {
                    'gating': False,
                    'key_dim': 64,
                    'num_head': 4,
                    'value_dim': 64,
                },
                'dgram_features': {
                    'min_bin': 3.25,
                    'max_bin': 50.75,
                    'num_bins': 39,
                },
                'embed_torsion_angles': False,
                'enabled': False,
                'template_pair_stack': {
                    'num_block': 2,
                    'triangle_attention_starting_node': {
                        'dropout_rate': 0.25,
                        'gating': True,
                        'key_dim': 64,
                        'num_head': 4,
                        'orientation': 'per_row',
                        'shared_dropout': True,
                        'value_dim': 64,
                    },
                    'triangle_attention_ending_node': {
                        'dropout_rate': 0.25,
                        'gating': True,
                        'key_dim': 64,
                        'num_head': 4,
                        'orientation': 'per_column',
                        'shared_dropout': True,
                        'value_dim': 64,
                    },
                    'triangle_multiplication_outgoing': {
                        'dropout_rate': 0.25,
                        'equation': 'ikc,jkc->ijc',
                        'num_intermediate_channel': 64,
                        'orientation': 'per_row',
                        'shared_dropout': True,
                        'fuse_projection_weights': False,
                    },
                    'triangle_multiplication_incoming': {
                        'dropout_rate': 0.25,
                        'equation': 'kjc,kic->ijc',
                        'num_intermediate_channel': 64,
                        'orientation': 'per_row',
                        'shared_dropout': True,
                        'fuse_projection_weights': False,
                    },
                    'pair_transition': {
                        'dropout_rate': 0.0,
                        'num_intermediate_factor': 2,
                        'orientation': 'per_row',
                        'shared_dropout': True,
                    },
                },
                'max_templates': 4,
                'subbatch_size': 128,
                'use_template_unit_vector': False,
            },
        },
        'global_config': {
            'deterministic': False,
            'multimer_mode': False,
            'subbatch_size': 4,
            'use_remat': False,
            'zero_init': True,
            'eval_dropout': False,
        },
        'heads': {
            'distogram': {
                'first_break': 2.3125,
                'last_break': 21.6875,
                'num_bins': 64,
                'weight': 0.3,
            },
            'predicted_aligned_error': {
                # `num_bins - 1` bins uniformly space the
                # [0, max_error_bin A] range.
                # The final bin covers [max_error_bin A, +infty]
                # 31A gives bins with 0.5A width.
                'max_error_bin': 31.0,
                'num_bins': 64,
                'num_channels': 128,
                'filter_by_resolution': True,
                'min_resolution': 0.1,
                'max_resolution': 3.0,
                'weight': 0.0,
            },
            'experimentally_resolved': {
                'filter_by_resolution': True,
                'max_resolution': 3.0,
                'min_resolution': 0.1,
                'weight': 0.01,
            },
            'structure_module': {
                'num_layer': 8,
                'fape': {
                    'clamp_distance': 10.0,
                    'clamp_type': 'relu',
                    'loss_unit_distance': 10.0,
                },
                'angle_norm_weight': 0.01,
                'chi_weight': 0.5,
                'clash_overlap_tolerance': 1.5,
                'compute_in_graph_metrics': True,
                'dropout': 0.1,
                'num_channel': 384,
                'num_head': 12,
                'num_layer_in_transition': 3,
                'num_point_qk': 4,
                'num_point_v': 8,
                'num_scalar_qk': 16,
                'num_scalar_v': 16,
                'position_scale': 10.0,
                'sidechain': {
                    'atom_clamp_distance': 10.0,
                    'num_channel': 128,
                    'num_residual_block': 2,
                    'weight_frac': 0.5,
                    'length_scale': 10.0,
                },
                'structural_violation_loss_weight': 1.0,
                'violation_tolerance_factor': 12.0,
                'weight': 1.0,
            },
            'predicted_lddt': {
                'filter_by_resolution': True,
                'max_resolution': 3.0,
                'min_resolution': 0.1,
                'num_bins': 50,
                'num_channels': 128,
                'weight': 0.01,
            },
            'masked_msa': {'num_output': 23, 'weight': 2.0},
        },
        'num_recycle': 3,
        'resample_msa_in_recycling': True,
    },
})


CONFIG_MULTIMER = ml_collections.ConfigDict({
    'model': {
        'embeddings_and_evoformer': {
            'evoformer_num_block': 48,
            'evoformer': {
                'msa_column_attention': {
                    'dropout_rate': 0.0,
                    'gating': True,
                    'num_head': 8,
                    'orientation': 'per_column',
                    'shared_dropout': True,
                },
                'msa_row_attention_with_pair_bias': {
                    'dropout_rate': 0.15,
                    'gating': True,
                    'num_head': 8,
                    'orientation': 'per_row',
                    'shared_dropout': True,
                },
                'msa_transition': {
                    'dropout_rate': 0.0,
                    'num_intermediate_factor': 4,
                    'orientation': 'per_row',
                    'shared_dropout': True,
                },
                'outer_product_mean': {
                    'chunk_size': 128,
                    'dropout_rate': 0.0,
                    'first': True,
                    'num_outer_channel': 32,
                    'orientation': 'per_row',
                    'shared_dropout': True,
                },
                'pair_transition': {
                    'dropout_rate': 0.0,
                    'num_intermediate_factor': 4,
                    'orientation': 'per_row',
                    'shared_dropout': True,
                },
                'triangle_attention_ending_node': {
                    'dropout_rate': 0.25,
                    'gating': True,
                    'num_head': 4,
                    'orientation': 'per_column',
                    'shared_dropout': True,
                },
                'triangle_attention_starting_node': {
                    'dropout_rate': 0.25,
                    'gating': True,
                    'num_head': 4,
                    'orientation': 'per_row',
                    'shared_dropout': True,
                },
                'triangle_multiplication_incoming': {
                    'dropout_rate': 0.25,
                    'equation': 'kjc,kic->ijc',
                    'num_intermediate_channel': 128,
                    'orientation': 'per_row',
                    'shared_dropout': True,
                    'fuse_projection_weights': True,
                },
                'triangle_multiplication_outgoing': {
                    'dropout_rate': 0.25,
                    'equation': 'ikc,jkc->ijc',
                    'num_intermediate_channel': 128,
                    'orientation': 'per_row',
                    'shared_dropout': True,
                    'fuse_projection_weights': True,
                },
            },
            'extra_msa_channel': 64,
            'extra_msa_stack_num_block': 4,
            'num_msa': 508,
            'num_extra_msa': 2048,
            'masked_msa': {
                'profile_prob': 0.1,
                'replace_fraction': 0.15,
                'same_prob': 0.1,
                'uniform_prob': 0.1,
            },
            'use_chain_relative': True,
            'max_relative_chain': 2,
            'max_relative_idx': 32,
            'seq_channel': 384,
            'msa_channel': 256,
            'pair_channel': 128,
            'prev_pos': {'max_bin': 20.75, 'min_bin': 3.25, 'num_bins': 15},
            'recycle_features': True,
            'recycle_pos': True,
            'template': {
                'attention': {'gating': False, 'num_head': 4},
                'dgram_features': {
                    'max_bin': 50.75,
                    'min_bin': 3.25,
                    'num_bins': 39,
                },
                'enabled': True,
                'max_templates': 4,
                'num_channels': 64,
                'subbatch_size': 128,
                'template_pair_stack': {
                    'num_block': 2,
                    'pair_transition': {
                        'dropout_rate': 0.0,
                        'num_intermediate_factor': 2,
                        'orientation': 'per_row',
                        'shared_dropout': True,
                    },
                    'triangle_attention_ending_node': {
                        'dropout_rate': 0.25,
                        'gating': True,
                        'num_head': 4,
                        'orientation': 'per_column',
                        'shared_dropout': True,
                    },
                    'triangle_attention_starting_node': {
                        'dropout_rate': 0.25,
                        'gating': True,
                        'num_head': 4,
                        'orientation': 'per_row',
                        'shared_dropout': True,
                    },
                    'triangle_multiplication_incoming': {
                        'dropout_rate': 0.25,
                        'equation': 'kjc,kic->ijc',
                        'num_intermediate_channel': 64,
                        'orientation': 'per_row',
                        'shared_dropout': True,
                        'fuse_projection_weights': True,
                    },
                    'triangle_multiplication_outgoing': {
                        'dropout_rate': 0.25,
                        'equation': 'ikc,jkc->ijc',
                        'num_intermediate_channel': 64,
                        'orientation': 'per_row',
                        'shared_dropout': True,
                        'fuse_projection_weights': True,
                    },
                },
            },
        },
        'global_config': {
            'bfloat16': True,
            'bfloat16_output': False,
            'deterministic': False,
            'multimer_mode': True,
            'subbatch_size': 4,
            'use_remat': False,
            'zero_init': True,
            'eval_dropout': False,
        },
        'heads': {
            'distogram': {
                'first_break': 2.3125,
                'last_break': 21.6875,
                'num_bins': 64,
                'weight': 0.3,
            },
            'experimentally_resolved': {
                'filter_by_resolution': True,
                'max_resolution': 3.0,
                'min_resolution': 0.1,
                'weight': 0.01,
            },
            'masked_msa': {'weight': 2.0},
            'predicted_aligned_error': {
                'filter_by_resolution': True,
                'max_error_bin': 31.0,
                'max_resolution': 3.0,
                'min_resolution': 0.1,
                'num_bins': 64,
                'num_channels': 128,
                'weight': 0.1,
            },
            'predicted_lddt': {
                'filter_by_resolution': True,
                'max_resolution': 3.0,
                'min_resolution': 0.1,
                'num_bins': 50,
                'num_channels': 128,
                'weight': 0.01,
            },
            'structure_module': {
                'angle_norm_weight': 0.01,
                'chi_weight': 0.5,
                'clash_overlap_tolerance': 1.5,
                'dropout': 0.1,
                'interface_fape': {
                    'atom_clamp_distance': 1000.0,
                    'loss_unit_distance': 20.0,
                },
                'intra_chain_fape': {
                    'atom_clamp_distance': 10.0,
                    'loss_unit_distance': 10.0,
                },
                'num_channel': 384,
                'num_head': 12,
                'num_layer': 8,
                'num_layer_in_transition': 3,
                'num_point_qk': 4,
                'num_point_v': 8,
                'num_scalar_qk': 16,
                'num_scalar_v': 16,
                'position_scale': 20.0,
                'sidechain': {
                    'atom_clamp_distance': 10.0,
                    'loss_unit_distance': 10.0,
                    'num_channel': 128,
                    'num_residual_block': 2,
                    'weight_frac': 0.5,
                },
                'structural_violation_loss_weight': 1.0,
                'violation_tolerance_factor': 12.0,
                'weight': 1.0,
            },
        },
        'num_ensemble_eval': 1,
        'num_recycle': 20,
        # A negative value indicates that no early stopping will occur, i.e.
        # the model will always run `num_recycle` number of recycling
        # iterations.  A positive value will enable early stopping if the
        # difference in pairwise distances is less than the tolerance between
        # recycling steps.
        'recycle_early_stop_tolerance': 0.5,
        'resample_msa_in_recycling': True,
    }
})


class MaskedMsa(base_config.BaseConfig):
  profile_prob: float
  same_prob: float
  uniform_prob: float
  replace_fraction: Optional[float] = None


class CommonData(base_config.BaseConfig):
  masked_msa: MaskedMsa
  max_extra_msa: int
  msa_cluster_features: bool
  num_recycle: int
  reduce_msa_clusters_by_max_templates: bool
  resample_msa_in_recycling: bool
  template_features: Sequence[str]
  unsupervised_features: Sequence[str]
  use_templates: bool


class EvalData(base_config.BaseConfig):
  feat: Dict[str, Any]
  fixed_size: bool
  subsample_templates: bool
  masked_msa_replace_fraction: float
  max_msa_clusters: int
  max_templates: int
  num_ensemble: int
  crop_size: Optional[int] = None


class Data(base_config.BaseConfig):
  common: CommonData
  eval: EvalData


class MsaRowAttentionWithPairBias(base_config.BaseConfig):
  dropout_rate: float
  gating: bool
  num_head: int
  orientation: str
  shared_dropout: bool


class MsaColumnAttention(base_config.BaseConfig):
  dropout_rate: float
  gating: bool
  num_head: int
  orientation: str
  shared_dropout: bool


class MsaTransition(base_config.BaseConfig):
  dropout_rate: float
  num_intermediate_factor: int
  orientation: str
  shared_dropout: bool


class OuterProductMean(base_config.BaseConfig):
  first: bool
  chunk_size: int
  dropout_rate: float
  num_outer_channel: int
  orientation: str
  shared_dropout: bool


class TriangleAttention(base_config.BaseConfig):
  dropout_rate: float
  gating: bool
  num_head: int
  orientation: str
  shared_dropout: bool


class TriangleMultiplication(base_config.BaseConfig):
  dropout_rate: float
  equation: str
  num_intermediate_channel: int
  orientation: str
  shared_dropout: bool
  fuse_projection_weights: bool


class PairTransition(base_config.BaseConfig):
  dropout_rate: float
  num_intermediate_factor: int
  orientation: str
  shared_dropout: bool


class Evoformer(base_config.BaseConfig):
  msa_row_attention_with_pair_bias: MsaRowAttentionWithPairBias
  msa_column_attention: MsaColumnAttention
  msa_transition: MsaTransition
  outer_product_mean: OuterProductMean
  triangle_attention_starting_node: TriangleAttention
  triangle_attention_ending_node: TriangleAttention
  triangle_multiplication_outgoing: TriangleMultiplication
  triangle_multiplication_incoming: TriangleMultiplication
  pair_transition: PairTransition


class TemplateAttention(base_config.BaseConfig):
  gating: bool
  num_head: int
  key_dim: Optional[int] = None
  value_dim: Optional[int] = None


class DgramFeatures(base_config.BaseConfig):
  min_bin: float
  max_bin: float
  num_bins: int


class TemplatePairStackAttention(base_config.BaseConfig):
  dropout_rate: float
  gating: bool
  num_head: int
  orientation: str
  shared_dropout: bool
  key_dim: Optional[int] = None
  value_dim: Optional[int] = None


class TemplatePairStackTriangleMultiplication(base_config.BaseConfig):
  dropout_rate: float
  equation: str
  num_intermediate_channel: int
  orientation: str
  shared_dropout: bool
  fuse_projection_weights: bool


class TemplatePairStackTransition(base_config.BaseConfig):
  dropout_rate: float
  num_intermediate_factor: int
  orientation: str
  shared_dropout: bool


class TemplatePairStack(base_config.BaseConfig):
  num_block: int
  triangle_attention_starting_node: TemplatePairStackAttention
  triangle_attention_ending_node: TemplatePairStackAttention
  triangle_multiplication_outgoing: TemplatePairStackTriangleMultiplication
  triangle_multiplication_incoming: TemplatePairStackTriangleMultiplication
  pair_transition: TemplatePairStackTransition


class Template(base_config.BaseConfig):
  attention: TemplateAttention
  dgram_features: DgramFeatures
  enabled: bool
  template_pair_stack: TemplatePairStack
  max_templates: int
  subbatch_size: int
  use_template_unit_vector: Optional[bool] = None
  embed_torsion_angles: Optional[bool] = None
  num_channels: Optional[int] = None


class PrevPos(base_config.BaseConfig):
  min_bin: float
  max_bin: float
  num_bins: int


class EmbeddingsAndEvoformer(base_config.BaseConfig):
  """Config for the embeddings and evoformer."""

  evoformer_num_block: int
  evoformer: Evoformer
  extra_msa_channel: int
  extra_msa_stack_num_block: int
  msa_channel: int
  pair_channel: int
  prev_pos: PrevPos
  recycle_features: bool
  recycle_pos: bool
  seq_channel: int
  template: Template
  max_relative_feature: Optional[int] = None
  num_msa: Optional[int] = None
  num_extra_msa: Optional[int] = None
  masked_msa: Optional[MaskedMsa] = None
  use_chain_relative: Optional[bool] = None
  max_relative_chain: Optional[int] = None
  max_relative_idx: Optional[int] = None


class GlobalConfig(base_config.BaseConfig):
  deterministic: bool
  multimer_mode: bool
  subbatch_size: int
  use_remat: bool
  zero_init: bool
  eval_dropout: bool
  bfloat16: Optional[bool] = None
  bfloat16_output: Optional[bool] = None


class DistogramHead(base_config.BaseConfig):
  first_break: float
  last_break: float
  num_bins: int
  weight: float


class PredictedAlignedErrorHead(base_config.BaseConfig):
  max_error_bin: float
  num_bins: int
  num_channels: int
  filter_by_resolution: bool
  min_resolution: float
  max_resolution: float
  weight: float


class ExperimentallyResolvedHead(base_config.BaseConfig):
  filter_by_resolution: bool
  max_resolution: float
  min_resolution: float
  weight: float


class Fape(base_config.BaseConfig):
  clamp_distance: float
  clamp_type: str
  loss_unit_distance: float


class Sidechain(base_config.BaseConfig):
  atom_clamp_distance: float
  num_channel: int
  num_residual_block: int
  weight_frac: float
  length_scale: Optional[float] = None
  loss_unit_distance: Optional[float] = None


class StructureModuleHead(base_config.BaseConfig):
  """Config for the structure module head."""

  num_layer: int
  angle_norm_weight: float
  chi_weight: float
  clash_overlap_tolerance: float
  dropout: float
  num_channel: int
  num_head: int
  num_layer_in_transition: int
  num_point_qk: int
  num_point_v: int
  num_scalar_qk: int
  num_scalar_v: int
  position_scale: float
  sidechain: Sidechain
  structural_violation_loss_weight: float
  violation_tolerance_factor: float
  weight: float
  fape: Optional[Fape] = None
  compute_in_graph_metrics: Optional[bool] = None
  interface_fape: Optional[Dict[str, float]] = None
  intra_chain_fape: Optional[Dict[str, float]] = None


class PredictedLDDTHead(base_config.BaseConfig):
  filter_by_resolution: bool
  max_resolution: float
  min_resolution: float
  num_bins: int
  num_channels: int
  weight: float


class MaskedMSAHead(base_config.BaseConfig):
  weight: float
  num_output: Optional[int] = None


class Heads(base_config.BaseConfig):
  distogram: DistogramHead
  predicted_aligned_error: PredictedAlignedErrorHead
  experimentally_resolved: ExperimentallyResolvedHead
  structure_module: StructureModuleHead
  predicted_lddt: PredictedLDDTHead
  masked_msa: MaskedMSAHead


class Model(base_config.BaseConfig):
  embeddings_and_evoformer: EmbeddingsAndEvoformer
  global_config: GlobalConfig
  heads: Heads
  num_recycle: int
  resample_msa_in_recycling: bool
  num_ensemble_eval: Optional[int] = None
  recycle_early_stop_tolerance: Optional[float] = None


class AlphaFoldConfig(base_config.BaseConfig):
  model: Model
  data: Optional[Data] = None


def model_config(name: str) -> ml_collections.ConfigDict:
  """Get the ConfigDict of a CASP14 model."""

  if name not in CONFIG_DIFFS:
    raise ValueError(f'Invalid model name {name}.')
  if 'multimer' in name:
    cfg = copy.deepcopy(CONFIG_MULTIMER)
  else:
    cfg = copy.deepcopy(CONFIG)
  cfg.update_from_flattened_dict(CONFIG_DIFFS[name])
  return cfg


@functools.lru_cache  #  Cache configs per model to avoid re-initializing.
def get_model_config(name: str, frozen: bool = True) -> AlphaFoldConfig:
  """Get the Config DataClass of a CASP14 model."""
  if name not in CONFIG_DIFFS:
    raise ValueError(f'Invalid model name {name}.')
  cfg = (
      AlphaFoldConfig(**(CONFIG_MULTIMER.to_dict()))
      if 'multimer' in name
      else AlphaFoldConfig(**(CONFIG.to_dict()))
  )
  apply_diff_op = CONFIG_DIFF_OPS[name]
  apply_diff_op(cfg)
  if frozen:
    cfg.freeze()
  return cfg


def _apply_model_1_diff(cfg: AlphaFoldConfig) -> None:
  if cfg.data:
    cfg.data.common.max_extra_msa = 5120
    cfg.data.common.reduce_msa_clusters_by_max_templates = True
    cfg.data.common.use_templates = True
  cfg.model.embeddings_and_evoformer.template.embed_torsion_angles = True
  cfg.model.embeddings_and_evoformer.template.enabled = True


def _apply_model_2_diff(cfg: AlphaFoldConfig) -> None:
  if cfg.data:
    cfg.data.common.reduce_msa_clusters_by_max_templates = True
    cfg.data.common.use_templates = True
  cfg.model.embeddings_and_evoformer.template.embed_torsion_angles = True
  cfg.model.embeddings_and_evoformer.template.enabled = True


def _apply_model_3_diff(cfg: AlphaFoldConfig) -> None:
  if cfg.data:
    cfg.data.common.max_extra_msa = 5120


def _apply_model_4_diff(cfg: AlphaFoldConfig) -> None:
  if cfg.data:
    cfg.data.common.max_extra_msa = 5120


def _apply_model_5_diff(cfg: AlphaFoldConfig) -> None:  # pylint: disable=unused-argument
  pass


def _apply_model_1_ptm_diff(cfg: AlphaFoldConfig) -> None:
  if cfg.data:
    cfg.data.common.max_extra_msa = 5120
    cfg.data.common.reduce_msa_clusters_by_max_templates = True
    cfg.data.common.use_templates = True
  cfg.model.embeddings_and_evoformer.template.embed_torsion_angles = True
  cfg.model.embeddings_and_evoformer.template.enabled = True
  cfg.model.heads.predicted_aligned_error.weight = 0.1


def _apply_model_2_ptm_diff(cfg: AlphaFoldConfig) -> None:
  if cfg.data:
    cfg.data.common.reduce_msa_clusters_by_max_templates = True
    cfg.data.common.use_templates = True
  cfg.model.embeddings_and_evoformer.template.embed_torsion_angles = True
  cfg.model.embeddings_and_evoformer.template.enabled = True
  cfg.model.heads.predicted_aligned_error.weight = 0.1


def _apply_model_3_ptm_diff(cfg: AlphaFoldConfig) -> None:
  if cfg.data:
    cfg.data.common.max_extra_msa = 5120
  cfg.model.heads.predicted_aligned_error.weight = 0.1


def _apply_model_4_ptm_diff(cfg: AlphaFoldConfig) -> None:
  if cfg.data:
    cfg.data.common.max_extra_msa = 5120
  cfg.model.heads.predicted_aligned_error.weight = 0.1


def _apply_model_5_ptm_diff(cfg: AlphaFoldConfig) -> None:
  cfg.model.heads.predicted_aligned_error.weight = 0.1


def _apply_model_1_multimer_v3_diff(cfg: AlphaFoldConfig) -> None:  # pylint: disable=unused-argument
  pass


def _apply_model_2_multimer_v3_diff(cfg: AlphaFoldConfig) -> None:  # pylint: disable=unused-argument
  pass


def _apply_model_3_multimer_v3_diff(cfg: AlphaFoldConfig) -> None:  # pylint: disable=unused-argument
  pass


def _apply_model_4_multimer_v3_diff(cfg: AlphaFoldConfig) -> None:
  cfg.model.embeddings_and_evoformer.num_extra_msa = 1152


def _apply_model_5_multimer_v3_diff(cfg: AlphaFoldConfig) -> None:
  cfg.model.embeddings_and_evoformer.num_extra_msa = 1152


def _common_updates(cfg: AlphaFoldConfig) -> None:
  """Applies common updates to the AlphaFold config."""
  cfg.model.embeddings_and_evoformer.num_msa = 252
  cfg.model.embeddings_and_evoformer.num_extra_msa = 1152
  cfg.model.embeddings_and_evoformer.evoformer.triangle_multiplication_incoming.fuse_projection_weights = (
      False
  )
  cfg.model.embeddings_and_evoformer.evoformer.triangle_multiplication_outgoing.fuse_projection_weights = (
      False
  )
  cfg.model.embeddings_and_evoformer.template.template_pair_stack.triangle_multiplication_incoming.fuse_projection_weights = (
      False
  )
  cfg.model.embeddings_and_evoformer.template.template_pair_stack.triangle_multiplication_outgoing.fuse_projection_weights = (
      False
  )


CONFIG_DIFF_OPS = {
    'model_1': _apply_model_1_diff,
    'model_2': _apply_model_2_diff,
    'model_3': _apply_model_3_diff,
    'model_4': _apply_model_4_diff,
    'model_5': _apply_model_5_diff,
    'model_1_ptm': _apply_model_1_ptm_diff,
    'model_2_ptm': _apply_model_2_ptm_diff,
    'model_3_ptm': _apply_model_3_ptm_diff,
    'model_4_ptm': _apply_model_4_ptm_diff,
    'model_5_ptm': _apply_model_5_ptm_diff,
    'model_1_multimer_v3': _apply_model_1_multimer_v3_diff,
    'model_2_multimer_v3': _apply_model_2_multimer_v3_diff,
    'model_3_multimer_v3': _apply_model_3_multimer_v3_diff,
    'model_4_multimer_v3': _apply_model_4_multimer_v3_diff,
    'model_5_multimer_v3': _apply_model_5_multimer_v3_diff,
}

CONFIG_DIFF_OPS.update(
    {f'model_{i}_multimer': _common_updates for i in range(1, 6)}
)
CONFIG_DIFF_OPS.update(
    {f'model_{i}_multimer_v2': _common_updates for i in range(1, 6)}
)
