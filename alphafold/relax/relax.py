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

"""Amber relaxation."""
from typing import Any, Dict, Sequence, Tuple
from alphafold.common import protein
from alphafold.relax import amber_minimize
from alphafold.relax import utils
import numpy as np


class AmberRelaxation(object):
  """Amber relaxation."""

  def __init__(self,
               *,
               max_iterations: int,
               tolerance: float,
               stiffness: float,
               exclude_residues: Sequence[int],
               max_outer_iterations: int):
    """Initialize Amber Relaxer.

    Args:
      max_iterations: Maximum number of L-BFGS iterations. 0 means no max.
      tolerance: kcal/mol, the energy tolerance of L-BFGS.
      stiffness: kcal/mol A**2, spring constant of heavy atom restraining
        potential.
      exclude_residues: Residues to exclude from per-atom restraining.
        Zero-indexed.
      max_outer_iterations: Maximum number of violation-informed relax
       iterations. A value of 1 will run the non-iterative procedure used in
       CASP14. Use 20 so that >95% of the bad cases are relaxed. Relax finishes
       as soon as there are no violations, hence in most cases this causes no
       slowdown. In the worst case we do 20 outer iterations.
    """

    self._max_iterations = max_iterations
    self._tolerance = tolerance
    self._stiffness = stiffness
    self._exclude_residues = exclude_residues
    self._max_outer_iterations = max_outer_iterations

  def process(self, *,
              prot: protein.Protein) -> Tuple[str, Dict[str, Any], np.ndarray]:
    """Runs Amber relax on a prediction, adds hydrogens, returns PDB string."""
    out = amber_minimize.run_pipeline(
        prot=prot, max_iterations=self._max_iterations,
        tolerance=self._tolerance, stiffness=self._stiffness,
        exclude_residues=self._exclude_residues,
        max_outer_iterations=self._max_outer_iterations)
    min_pos = out['pos']
    start_pos = out['posinit']
    rmsd = np.sqrt(np.sum((start_pos - min_pos)**2) / start_pos.shape[0])
    debug_data = {
        'initial_energy': out['einit'],
        'final_energy': out['efinal'],
        'attempts': out['min_attempts'],
        'rmsd': rmsd
    }
    pdb_str = amber_minimize.clean_protein(prot)
    min_pdb = utils.overwrite_pdb_coordinates(pdb_str, min_pos)
    min_pdb = utils.overwrite_b_factors(min_pdb, prot.b_factors)
    utils.assert_equal_nonterminal_atom_types(
        protein.from_pdb_string(min_pdb).atom_mask,
        prot.atom_mask)
    violations = out['structural_violations'][
        'total_per_residue_violations_mask']
    return min_pdb, debug_data, violations
