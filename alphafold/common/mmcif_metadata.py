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

"""mmCIF metadata."""

from typing import Mapping, Sequence
from alphafold import version
import numpy as np


_DISCLAIMER = """ALPHAFOLD DATA, COPYRIGHT (2021) DEEPMIND TECHNOLOGIES LIMITED.
THE INFORMATION PROVIDED IS THEORETICAL MODELLING ONLY AND CAUTION SHOULD BE
EXERCISED IN ITS USE. IT IS PROVIDED "AS-IS" WITHOUT ANY WARRANTY OF ANY KIND,
WHETHER EXPRESSED OR IMPLIED. NO WARRANTY IS GIVEN THAT USE OF THE INFORMATION
SHALL NOT INFRINGE THE RIGHTS OF ANY THIRD PARTY. DISCLAIMER: THE INFORMATION IS
NOT INTENDED TO BE A SUBSTITUTE FOR PROFESSIONAL MEDICAL ADVICE, DIAGNOSIS, OR
TREATMENT, AND DOES NOT CONSTITUTE MEDICAL OR OTHER PROFESSIONAL ADVICE. IT IS
AVAILABLE FOR ACADEMIC AND COMMERCIAL PURPOSES, UNDER CC-BY 4.0 LICENCE."""

# Authors of the Nature methods paper we reference in the mmCIF.
_MMCIF_PAPER_AUTHORS = (
    'Jumper, John',
    'Evans, Richard',
    'Pritzel, Alexander',
    'Green, Tim',
    'Figurnov, Michael',
    'Ronneberger, Olaf',
    'Tunyasuvunakool, Kathryn',
    'Bates, Russ',
    'Zidek, Augustin',
    'Potapenko, Anna',
    'Bridgland, Alex',
    'Meyer, Clemens',
    'Kohl, Simon A. A.',
    'Ballard, Andrew J.',
    'Cowie, Andrew',
    'Romera-Paredes, Bernardino',
    'Nikolov, Stanislav',
    'Jain, Rishub',
    'Adler, Jonas',
    'Back, Trevor',
    'Petersen, Stig',
    'Reiman, David',
    'Clancy, Ellen',
    'Zielinski, Michal',
    'Steinegger, Martin',
    'Pacholska, Michalina',
    'Berghammer, Tamas',
    'Silver, David',
    'Vinyals, Oriol',
    'Senior, Andrew W.',
    'Kavukcuoglu, Koray',
    'Kohli, Pushmeet',
    'Hassabis, Demis',
)

# Authors of the mmCIF - we set them to be equal to the authors of the paper.
_MMCIF_AUTHORS = _MMCIF_PAPER_AUTHORS


def add_metadata_to_mmcif(
    old_cif: Mapping[str, Sequence[str]], model_type: str
) -> Mapping[str, Sequence[str]]:
  """Adds AlphaFold metadata in the given mmCIF."""
  cif = {}

  # ModelCIF conformation dictionary.
  cif['_audit_conform.dict_name'] = ['mmcif_ma.dic']
  cif['_audit_conform.dict_version'] = ['1.3.9']
  cif['_audit_conform.dict_location'] = [
      'https://raw.githubusercontent.com/ihmwg/ModelCIF/master/dist/'
      'mmcif_ma.dic'
  ]

  # License and disclaimer.
  cif['_pdbx_data_usage.id'] = ['1', '2']
  cif['_pdbx_data_usage.type'] = ['license', 'disclaimer']
  cif['_pdbx_data_usage.details'] = [
      'Data in this file is available under a CC-BY-4.0 license.',
      _DISCLAIMER,
  ]
  cif['_pdbx_data_usage.url'] = [
      'https://creativecommons.org/licenses/by/4.0/',
      '?',
  ]
  cif['_pdbx_data_usage.name'] = ['CC-BY-4.0', '?']

  # Structure author details.
  cif['_audit_author.name'] = []
  cif['_audit_author.pdbx_ordinal'] = []
  for author_index, author_name in enumerate(_MMCIF_AUTHORS, start=1):
    cif['_audit_author.name'].append(author_name)
    cif['_audit_author.pdbx_ordinal'].append(str(author_index))

  # Paper author details.
  cif['_citation_author.citation_id'] = []
  cif['_citation_author.name'] = []
  cif['_citation_author.ordinal'] = []
  for author_index, author_name in enumerate(_MMCIF_PAPER_AUTHORS, start=1):
    cif['_citation_author.citation_id'].append('primary')
    cif['_citation_author.name'].append(author_name)
    cif['_citation_author.ordinal'].append(str(author_index))

  # Paper citation details.
  cif['_citation.id'] = ['primary']
  cif['_citation.title'] = [
      'Highly accurate protein structure prediction with AlphaFold'
  ]
  cif['_citation.journal_full'] = ['Nature']
  cif['_citation.journal_volume'] = ['596']
  cif['_citation.page_first'] = ['583']
  cif['_citation.page_last'] = ['589']
  cif['_citation.year'] = ['2021']
  cif['_citation.journal_id_ASTM'] = ['NATUAS']
  cif['_citation.country'] = ['UK']
  cif['_citation.journal_id_ISSN'] = ['0028-0836']
  cif['_citation.journal_id_CSD'] = ['0006']
  cif['_citation.book_publisher'] = ['?']
  cif['_citation.pdbx_database_id_PubMed'] = ['34265844']
  cif['_citation.pdbx_database_id_DOI'] = ['10.1038/s41586-021-03819-2']

  # Type of data in the dataset including data used in the model generation.
  cif['_ma_data.id'] = ['1']
  cif['_ma_data.name'] = ['Model']
  cif['_ma_data.content_type'] = ['model coordinates']

  # Description of number of instances for each entity.
  cif['_ma_target_entity_instance.asym_id'] = old_cif['_struct_asym.id']
  cif['_ma_target_entity_instance.entity_id'] = old_cif[
      '_struct_asym.entity_id'
  ]
  cif['_ma_target_entity_instance.details'] = ['.'] * len(
      cif['_ma_target_entity_instance.entity_id']
  )

  # Details about the target entities.
  cif['_ma_target_entity.entity_id'] = cif[
      '_ma_target_entity_instance.entity_id'
  ]
  cif['_ma_target_entity.data_id'] = ['1'] * len(
      cif['_ma_target_entity.entity_id']
  )
  cif['_ma_target_entity.origin'] = ['.'] * len(
      cif['_ma_target_entity.entity_id']
  )

  # Details of the models being deposited.
  cif['_ma_model_list.ordinal_id'] = ['1']
  cif['_ma_model_list.model_id'] = ['1']
  cif['_ma_model_list.model_group_id'] = ['1']
  cif['_ma_model_list.model_name'] = ['Top ranked model']

  cif['_ma_model_list.model_group_name'] = [
      f'AlphaFold {model_type} v{version.__version__} model'
  ]
  cif['_ma_model_list.data_id'] = ['1']
  cif['_ma_model_list.model_type'] = ['Ab initio model']

  # Software used.
  cif['_software.pdbx_ordinal'] = ['1']
  cif['_software.name'] = ['AlphaFold']
  cif['_software.version'] = [f'v{version.__version__}']
  cif['_software.type'] = ['package']
  cif['_software.description'] = ['Structure prediction']
  cif['_software.classification'] = ['other']
  cif['_software.date'] = ['?']

  # Collection of software into groups.
  cif['_ma_software_group.ordinal_id'] = ['1']
  cif['_ma_software_group.group_id'] = ['1']
  cif['_ma_software_group.software_id'] = ['1']

  # Method description to conform with ModelCIF.
  cif['_ma_protocol_step.ordinal_id'] = ['1', '2', '3']
  cif['_ma_protocol_step.protocol_id'] = ['1', '1', '1']
  cif['_ma_protocol_step.step_id'] = ['1', '2', '3']
  cif['_ma_protocol_step.method_type'] = [
      'coevolution MSA',
      'template search',
      'modeling',
  ]

  # Details of the metrics use to assess model confidence.
  cif['_ma_qa_metric.id'] = ['1', '2']
  cif['_ma_qa_metric.name'] = ['pLDDT', 'pLDDT']
  # Accepted values are distance, energy, normalised score, other, zscore.
  cif['_ma_qa_metric.type'] = ['pLDDT', 'pLDDT']
  cif['_ma_qa_metric.mode'] = ['global', 'local']
  cif['_ma_qa_metric.software_group_id'] = ['1', '1']

  # Global model confidence metric value.
  cif['_ma_qa_metric_global.ordinal_id'] = ['1']
  cif['_ma_qa_metric_global.model_id'] = ['1']
  cif['_ma_qa_metric_global.metric_id'] = ['1']
  global_plddt = np.mean(
      [float(v) for v in old_cif['_atom_site.B_iso_or_equiv']]
  )
  cif['_ma_qa_metric_global.metric_value'] = [f'{global_plddt:.2f}']

  cif['_atom_type.symbol'] = sorted(set(old_cif['_atom_site.type_symbol']))

  return cif
