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

"""Tests for notebook_utils."""
import io

from absl.testing import absltest
from absl.testing import parameterized
from alphafold.data import parsers
from alphafold.data import templates
from alphafold.notebooks import notebook_utils

import mock
import numpy as np


ONLY_QUERY_HIT = {
    'sto': (
        '# STOCKHOLM 1.0\n'
        '#=GF ID query-i1\n'
        'query   MAAHKGAEHHHKAAEHHEQAAKHHHAAAEHHEKGEHEQAAHHADTAYAHHKHAEEH\n'
        '//\n'),
    'tbl': '',
    'stderr': b'',
    'n_iter': 1,
    'e_value': 0.0001}

# pylint: disable=line-too-long
MULTI_SEQUENCE_HIT_1 = {
    'sto': (
        '# STOCKHOLM 1.0\n'
        '#=GF ID query-i1\n'
        '#=GS ERR1700680_4602609/41-109  DE [subseq from] ERR1700680_4602609\n'
        '#=GS ERR1019366_5760491/40-105  DE [subseq from] ERR1019366_5760491\n'
        '#=GS SRR5580704_12853319/61-125 DE [subseq from] SRR5580704_12853319\n'
        'query                              MAAHKGAEHHHKAAEHHEQAAKHHHAAAEHHEKGEHEQAAHHADTAYAHHKHAEEHAAQAAKHDAEHHAPKPH\n'
        'ERR1700680_4602609/41-109          --INKGAEYHKKAAEHHELAAKHHREAAKHHEAGSHEKAAHHSEIAAGHGLTAVHHTEEATK-HHPEEHTEK--\n'
        'ERR1019366_5760491/40-105          ---RSGAQHHDAAAQHYEEAARHHRMAAKQYQASHHEKAAHYAQLAYAHHMYAEQHAAEAAK-AHAKNHG----\n'
        'SRR5580704_12853319/61-125         ----PAADHHMKAAEHHEEAAKHHRAAAEHHTAGDHQKAGHHAHVANGHHVNAVHHAEEASK-HHATDHS----\n'
        '//\n'),
    'tbl': (
        'ERR1700680_4602609   -          query                -            7.7e-09   47.7  33.8   1.1e-08   47.2  33.8   1.2   1   0   0   1   1   1   1 -\n'
        'ERR1019366_5760491   -          query                -            1.7e-08   46.6  33.1   2.5e-08   46.1  33.1   1.3   1   0   0   1   1   1   1 -\n'
        'SRR5580704_12853319  -          query                -            1.1e-07   44.0  41.6     2e-07   43.1  41.6   1.4   1   0   0   1   1   1   1 -\n'),
        'stderr': b'',
        'n_iter': 1,
        'e_value': 0.0001}

MULTI_SEQUENCE_HIT_2 = {
    'sto': (
        '# STOCKHOLM 1.0\n'
        '#=GF ID query-i1\n'
        '#=GS ERR1700719_3476944/70-137  DE [subseq from] ERR1700719_3476944\n'
        '#=GS ERR1700761_4254522/72-138  DE [subseq from] ERR1700761_4254522\n'
        '#=GS SRR5438477_9761204/64-132  DE [subseq from] SRR5438477_9761204\n'
        'query                              MAAHKGAEHHHKAAEHHEQAAKHHHAAAEHHEKGEHEQAAHHADTAYAHHKHAEEHAAQAAKHDAEHHAPKPH\n'
        'ERR1700719_3476944/70-137          ---KQAAEHHHQAAEHHEHAARHHREAAKHHEAGDHESAAHHAHTAQGHLHQATHHASEAAKLHVEHHGQK--\n'
        'ERR1700761_4254522/72-138          ----QASEHHNLAAEHHEHAARHHRDAAKHHKAGDHEKAAHHAHVAHGHHLHATHHATEAAKHHVEAHGEK--\n'
        'SRR5438477_9761204/64-132          MPKHEGAEHHKKAAEHNEHAARHHKEAARHHEEGSHEKVGHHAHIAHGHHLHATHHAEEAAKTHSNQHE----\n'
        '//\n'),
    'tbl': (
        'ERR1700719_3476944   -          query                -              2e-07   43.2  47.5   3.5e-07   42.4  47.5   1.4   1   0   0   1   1   1   1 -\n'
        'ERR1700761_4254522   -          query                -            6.1e-07   41.6  48.1   8.1e-07   41.3  48.1   1.2   1   0   0   1   1   1   1 -\n'
        'SRR5438477_9761204   -          query                -            1.8e-06   40.2  46.9   2.3e-06   39.8  46.9   1.2   1   0   0   1   1   1   1 -\n'),
        'stderr': b'',
        'n_iter': 1,
        'e_value': 0.0001}
# pylint: enable=line-too-long


class NotebookUtilsTest(parameterized.TestCase):

  @parameterized.parameters(
      ('DeepMind', 'DEEPMIND'), ('A ', 'A'), ('\tA', 'A'), (' A\t\n', 'A'),
      ('ACDEFGHIKLMNPQRSTVWY', 'ACDEFGHIKLMNPQRSTVWY'))
  def test_clean_and_validate_sequence_ok(self, sequence, exp_clean):
    clean = notebook_utils.clean_and_validate_sequence(
        sequence, min_length=1, max_length=100)
    self.assertEqual(clean, exp_clean)

  @parameterized.named_parameters(
      ('too_short', 'AA', 'too short'),
      ('too_long', 'AAAAAAAAAA', 'too long'),
      ('bad_amino_acids_B', 'BBBB', 'non-amino acid'),
      ('bad_amino_acids_J', 'JJJJ', 'non-amino acid'),
      ('bad_amino_acids_O', 'OOOO', 'non-amino acid'),
      ('bad_amino_acids_U', 'UUUU', 'non-amino acid'),
      ('bad_amino_acids_X', 'XXXX', 'non-amino acid'),
      ('bad_amino_acids_Z', 'ZZZZ', 'non-amino acid'))
  def test_clean_and_validate_sequence_bad(self, sequence, exp_error):
    with self.assertRaisesRegex(ValueError, f'.*{exp_error}.*'):
      notebook_utils.clean_and_validate_sequence(
          sequence, min_length=4, max_length=8)

  @parameterized.parameters(
      (['A', '', '', ' ', '\t', ' \t\n', '', ''], ['A'],
       notebook_utils.ModelType.MONOMER),
      (['', 'A'], ['A'],
       notebook_utils.ModelType.MONOMER),
      (['A', 'C ', ''], ['A', 'C'],
       notebook_utils.ModelType.MULTIMER),
      (['', 'A', '', 'C '], ['A', 'C'],
       notebook_utils.ModelType.MULTIMER))
  def test_validate_input_ok(
      self, input_sequences, exp_sequences, exp_model_type):
    sequences, model_type = notebook_utils.validate_input(
        input_sequences=input_sequences,
        min_length=1, max_length=100, max_multimer_length=100)
    self.assertSequenceEqual(sequences, exp_sequences)
    self.assertEqual(model_type, exp_model_type)

  @parameterized.named_parameters(
      ('no_input_sequence', ['', '\t', '\n'], 'No input amino acid sequence'),
      ('too_long_single', ['AAAAAAAAA', 'AAAA'], 'Input sequence is too long'),
      ('too_long_multimer', ['AAAA', 'AAAAA'], 'The total length of multimer'))
  def test_validate_input_bad(self, input_sequences, exp_error):
    with self.assertRaisesRegex(ValueError, f'.*{exp_error}.*'):
      notebook_utils.validate_input(
          input_sequences=input_sequences,
          min_length=4, max_length=8, max_multimer_length=6)

  def test_merge_chunked_msa_no_hits(self):
    results = [ONLY_QUERY_HIT, ONLY_QUERY_HIT]
    merged_msa = notebook_utils.merge_chunked_msa(
        results=results)
    self.assertSequenceEqual(
        merged_msa.sequences,
        ('MAAHKGAEHHHKAAEHHEQAAKHHHAAAEHHEKGEHEQAAHHADTAYAHHKHAEEH',))
    self.assertSequenceEqual(merged_msa.deletion_matrix, ([0] * 56,))

  def test_merge_chunked_msa(self):
    results = [MULTI_SEQUENCE_HIT_1, MULTI_SEQUENCE_HIT_2]
    merged_msa = notebook_utils.merge_chunked_msa(
        results=results)
    self.assertLen(merged_msa.sequences, 7)
    # The 1st one is the query.
    self.assertEqual(
        merged_msa.sequences[0],
        'MAAHKGAEHHHKAAEHHEQAAKHHHAAAEHHEKGEHEQAAHHADTAYAHHKHAEEHAAQAAKHDAEHHAP'
        'KPH')
    # The 2nd one is the one with the lowest e-value: ERR1700680_4602609.
    self.assertEqual(
        merged_msa.sequences[1],
        '--INKGAEYHKKAAEHHELAAKHHREAAKHHEAGSHEKAAHHSEIAAGHGLTAVHHTEEATK-HHPEEHT'
        'EK-')
    # The last one is the one with the largest e-value: SRR5438477_9761204.
    self.assertEqual(
        merged_msa.sequences[-1],
        'MPKHEGAEHHKKAAEHNEHAARHHKEAARHHEEGSHEKVGHHAHIAHGHHLHATHHAEEAAKTHSNQHE-'
        '---')
    self.assertLen(merged_msa.deletion_matrix, 7)

  @mock.patch('sys.stdout', new_callable=io.StringIO)
  def test_show_msa_info(self, mocked_stdout):
    single_chain_msas = [
        parsers.Msa(sequences=['A', 'B', 'C', 'C'],
                    deletion_matrix=[None] * 4,
                    descriptions=[''] * 4),
        parsers.Msa(sequences=['A', 'A', 'A', 'D'],
                    deletion_matrix=[None] * 4,
                    descriptions=[''] * 4)
    ]
    notebook_utils.show_msa_info(
        single_chain_msas=single_chain_msas, sequence_index=1)
    self.assertEqual(mocked_stdout.getvalue(),
                     '\n4 unique sequences found in total for sequence 1\n\n')

  @parameterized.named_parameters(
      ('some_templates', 4), ('no_templates', 0))
  def test_empty_placeholder_template_features(self, num_templates):
    template_features = notebook_utils.empty_placeholder_template_features(
        num_templates=num_templates, num_res=16)
    self.assertCountEqual(template_features.keys(),
                          templates.TEMPLATE_FEATURES.keys())
    self.assertSameElements(
        [v.shape[0] for v in template_features.values()], [num_templates])
    self.assertSequenceEqual(
        [t.dtype for t in template_features.values()],
        [np.array([], dtype=templates.TEMPLATE_FEATURES[feat_name]).dtype
         for feat_name in template_features])

  def test_get_pae_json(self):
    pae = np.array([[0.01, 13.12345], [20.0987, 0.0]])
    pae_json = notebook_utils.get_pae_json(pae=pae, max_pae=31.75)
    self.assertEqual(
        pae_json,
        '[{"residue1":[1,1,2,2],"residue2":[1,2,1,2],"distance":'
        '[0.0,13.1,20.1,0.0],"max_predicted_aligned_error":31.75}]')


if __name__ == '__main__':
  absltest.main()
