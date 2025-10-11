# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# You may not use this file except in compliance with the License.
# You may obtain a copy at:
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.



import json
import os
import tempfile
import glob
from unittest import mock
from absl.testing import absltest, parameterized
import numpy as np
import run_alphafold


TEST_DATA_DIR = 'alphafold/common/testdata/'


class RunAlphafoldTest(parameterized.TestCase):
    """End-to-end AlphaFold structure prediction test suite."""

    def setUp(self):
        """Configure mocks and reusable test setup."""
        super().setUp()
        self.data_pipeline_mock = mock.Mock()
        self.model_runner_mock = mock.Mock()
        self.amber_relaxer_mock = mock.Mock()
        self._configure_mocks()

    def _configure_mocks(self):
        """Initialize default mock return values."""
        self.data_pipeline_mock.process.return_value = {}
        self.model_runner_mock.process_features.return_value = {
            'aatype': np.zeros((12, 10), dtype=np.int32),
            'residue_index': np.tile(np.arange(10, dtype=np.int32)[None], (12, 1)),
        }
        self.model_runner_mock.predict.return_value = {
            'structure_module': {
                'final_atom_positions': np.zeros((10, 37, 3)),
                'final_atom_mask': np.ones((10, 37)),
            },
            'predicted_lddt': {'logits': np.ones((10, 50))},
            'plddt': np.ones(10) * 42,
            'ranking_confidence': 90,
            'ptm': np.array(0.0),
            'aligned_confidence_probs': np.zeros((10, 10, 50)),
            'predicted_aligned_error': np.zeros((10, 10)),
            'max_predicted_aligned_error': np.array(0.0),
        }
        self.model_runner_mock.multimer_mode = False

    def _create_fasta(self, out_dir: str) -> str:
        """Helper to create a small FASTA file."""
        fasta_path = os.path.join(out_dir, 'target.fasta')
        with open(fasta_path, 'wt') as f:
            f.write('>A\nAAAAAAAAAAAAA')
        return fasta_path

    def _read_test_pdb(self, pdb_filename: str) -> str:
        """Helper to load a test PDB structure."""
        with open(
            os.path.join(
                absltest.get_default_test_srcdir(),
                TEST_DATA_DIR,
                pdb_filename
            )
        ) as f:
            return f.read()

    @parameterized.named_parameters(
        ('relax_monomer', run_alphafold.ModelsToRelax.ALL, 'Monomer'),
        ('no_relax_monomer', run_alphafold.ModelsToRelax.NONE, 'Monomer'),
        ('relax_multimer', run_alphafold.ModelsToRelax.ALL, 'Multimer'),
    )
    def test_end_to_end(self, models_to_relax, model_type):
        """Ensures AlphaFold runs from FASTA → PDB end-to-end correctly."""
        out_dir = tempfile.mkdtemp()
        fasta_path = self._create_fasta(out_dir)
        pdb_string = self._read_test_pdb('glucagon.pdb')

        self.amber_relaxer_mock.process.return_value = (
            pdb_string, None, [1.0, 0.0, 0.0]
        )

        fasta_name = 'test'

        # Run AlphaFold prediction
        run_alphafold.predict_structure(
            fasta_path=fasta_path,
            fasta_name=fasta_name,
            output_dir_base=out_dir,
            data_pipeline=self.data_pipeline_mock,
            model_runners={'model1': self.model_runner_mock},
            amber_relaxer=self.amber_relaxer_mock,
            benchmark=False,
            random_seed=0,
            models_to_relax=models_to_relax,
            model_type=model_type,
        )

        # Validate output directory structure
        base_files = os.listdir(out_dir)
        self.assertIn('target.fasta', base_files, msg="FASTA file missing in output.")
        self.assertIn('test', base_files, msg="Model output folder missing.")

        target_dir = os.path.join(out_dir, 'test')
        output_files = os.listdir(target_dir)

        expected_files = [
            'confidence_model1.json',
            'features.pkl',
            'msas',
            'pae_model1.json',
            'ranked_0.cif',
            'ranked_0.pdb',
            'ranking_debug.json',
            'result_model1.pkl',
            'timings.json',
            'unrelaxed_model1.cif',
            'unrelaxed_model1.pdb',
        ]

        if models_to_relax == run_alphafold.ModelsToRelax.ALL:
            expected_files += [
                'relaxed_model1.cif',
                'relaxed_model1.pdb',
                'relax_metrics.json',
            ]

            with open(os.path.join(target_dir, 'relax_metrics.json')) as f:
                relax_metrics = json.load(f)

            self.assertDictEqual(
                {
                    'model1': {
                        'remaining_violations': [1.0, 0.0, 0.0],
                        'remaining_violations_count': 1.0,
                    }
                },
                relax_metrics,
                msg="Relaxation metrics content mismatch.",
            )

        # File validation
        for f in expected_files:
            file_path = os.path.join(target_dir, f)
            self.assertTrue(os.path.exists(file_path), msg=f"Missing: {f}")
            # Ensure non-empty files (except msas folder)
            if not os.path.isdir(file_path):
                self.assertGreater(os.path.getsize(file_path), 0, msg=f"Empty file: {f}")

        # Validate B-factor correctness (pLDDT stored)
        pdb_path = os.path.join(target_dir, 'unrelaxed_model1.pdb')
        with open(pdb_path) as f:
            atom_lines = [line for line in f if line.startswith('ATOM')]
        self.assertTrue(atom_lines, msg="No ATOM records found in PDB file.")
        for line in atom_lines:
            self.assertEqual(line[61:66], '42.00', msg="Incorrect pLDDT value in B-factor.")

    def test_invalid_fasta_path_raises_error(self):
        """Verify that an invalid FASTA path raises a FileNotFoundError."""
        with self.assertRaises(FileNotFoundError):
            run_alphafold.predict_structure(
                fasta_path='invalid_path.fasta',
                fasta_name='invalid',
                output_dir_base='.',
                data_pipeline=self.data_pipeline_mock,
                model_runners={'model1': self.model_runner_mock},
                amber_relaxer=self.amber_relaxer_mock,
                benchmark=False,
                random_seed=0,
                models_to_relax=run_alphafold.ModelsToRelax.NONE,
                model_type='Monomer',
            )


if __name__ == '__main__':
    absltest.main()
