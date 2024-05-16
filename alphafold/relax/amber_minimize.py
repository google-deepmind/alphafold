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

"""Restrained Amber Minimization of a structure."""

import io
import time
from typing import Collection, Optional, Sequence

from absl import logging
from alphafold.common import protein
from alphafold.common import residue_constants
from alphafold.model import folding
from alphafold.relax import cleanup
from alphafold.relax import utils
import ml_collections
import numpy as np
import jax
import openmm
from openmm import unit
from openmm import app as openmm_app
from openmm.app.internal.pdbstructure import PdbStructure


ENERGY = unit.kilocalories_per_mole
LENGTH = unit.angstroms


def will_restrain(atom: openmm_app.Atom, rset: str) -> bool:
  """Returns True if the atom will be restrained by the given restraint set."""

  if rset == "non_hydrogen":
    return atom.element.name != "hydrogen"
  elif rset == "c_alpha":
    return atom.name == "CA"


def _add_restraints(
    system: openmm.System,
    reference_pdb: openmm_app.PDBFile,
    stiffness: unit.Unit,
    rset: str,
    exclude_residues: Sequence[int]):
  """Adds a harmonic potential that restrains the system to a structure."""
  assert rset in ["non_hydrogen", "c_alpha"]

  force = openmm.CustomExternalForce(
      "0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
  force.addGlobalParameter("k", stiffness)
  for p in ["x0", "y0", "z0"]:
    force.addPerParticleParameter(p)

  for i, atom in enumerate(reference_pdb.topology.atoms()):
    if atom.residue.index in exclude_residues:
      continue
    if will_restrain(atom, rset):
      force.addParticle(i, reference_pdb.positions[i])
  logging.info("Restraining %d / %d particles.",
               force.getNumParticles(), system.getNumParticles())
  system.addForce(force)


def _openmm_minimize(
    pdb_str: str,
    max_iterations: int,
    tolerance: unit.Unit,
    stiffness: unit.Unit,
    restraint_set: str,
    exclude_residues: Sequence[int],
    use_gpu: bool):
  """Minimize energy via openmm."""

  pdb_file = io.StringIO(pdb_str)
  pdb = openmm_app.PDBFile(pdb_file)

  force_field = openmm_app.ForceField("amber99sb.xml")
  constraints = openmm_app.HBonds
  system = force_field.createSystem(
      pdb.topology, constraints=constraints)
  if stiffness > 0 * ENERGY / (LENGTH**2):
    _add_restraints(system, pdb, stiffness, restraint_set, exclude_residues)

  integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)
  platform = openmm.Platform.getPlatformByName("CUDA" if use_gpu else "CPU")
  simulation = openmm_app.Simulation(
      pdb.topology, system, integrator, platform)
  simulation.context.setPositions(pdb.positions)

  ret = {}
  state = simulation.context.getState(getEnergy=True, getPositions=True)
  ret["einit"] = state.getPotentialEnergy().value_in_unit(ENERGY)
  ret["posinit"] = state.getPositions(asNumpy=True).value_in_unit(LENGTH)
  simulation.minimizeEnergy(maxIterations=max_iterations,
                            tolerance=tolerance)
  state = simulation.context.getState(getEnergy=True, getPositions=True)
  ret["efinal"] = state.getPotentialEnergy().value_in_unit(ENERGY)
  ret["pos"] = state.getPositions(asNumpy=True).value_in_unit(LENGTH)
  ret["min_pdb"] = _get_pdb_string(simulation.topology, state.getPositions())
  return ret


def _get_pdb_string(topology: openmm_app.Topology, positions: unit.Quantity):
  """Returns a pdb string provided OpenMM topology and positions."""
  with io.StringIO() as f:
    openmm_app.PDBFile.writeFile(topology, positions, f)
    return f.getvalue()


def _check_cleaned_atoms(pdb_cleaned_string: str, pdb_ref_string: str):
  """Checks that no atom positions have been altered by cleaning."""
  cleaned = openmm_app.PDBFile(io.StringIO(pdb_cleaned_string))
  reference = openmm_app.PDBFile(io.StringIO(pdb_ref_string))

  cl_xyz = np.array(cleaned.getPositions().value_in_unit(LENGTH))
  ref_xyz = np.array(reference.getPositions().value_in_unit(LENGTH))

  for ref_res, cl_res in zip(reference.topology.residues(),
                             cleaned.topology.residues()):
    assert ref_res.name == cl_res.name
    for rat in ref_res.atoms():
      for cat in cl_res.atoms():
        if cat.name == rat.name:
          if not np.array_equal(cl_xyz[cat.index], ref_xyz[rat.index]):
            raise ValueError(f"Coordinates of cleaned atom {cat} do not match "
                             f"coordinates of reference atom {rat}.")


def _check_residues_are_well_defined(prot: protein.Protein):
  """Checks that all residues contain non-empty atom sets."""
  if (prot.atom_mask.sum(axis=-1) == 0).any():
    raise ValueError("Amber minimization can only be performed on proteins with"
                     " well-defined residues. This protein contains at least"
                     " one residue with no atoms.")


def _check_atom_mask_is_ideal(prot):
  """Sanity-check the atom mask is ideal, up to a possible OXT."""
  atom_mask = prot.atom_mask
  ideal_atom_mask = protein.ideal_atom_mask(prot)
  utils.assert_equal_nonterminal_atom_types(atom_mask, ideal_atom_mask)


def clean_protein(
    prot: protein.Protein,
    checks: bool = True):
  """Adds missing atoms to Protein instance.

  Args:
    prot: A `protein.Protein` instance.
    checks: A `bool` specifying whether to add additional checks to the cleaning
      process.

  Returns:
    pdb_string: A string of the cleaned protein.
  """
  _check_atom_mask_is_ideal(prot)

  # Clean pdb.
  prot_pdb_string = protein.to_pdb(prot)
  pdb_file = io.StringIO(prot_pdb_string)
  alterations_info = {}
  fixed_pdb = cleanup.fix_pdb(pdb_file, alterations_info)
  fixed_pdb_file = io.StringIO(fixed_pdb)
  pdb_structure = PdbStructure(fixed_pdb_file)
  cleanup.clean_structure(pdb_structure, alterations_info)

  logging.info("alterations info: %s", alterations_info)

  # Write pdb file of cleaned structure.
  as_file = openmm_app.PDBFile(pdb_structure)
  pdb_string = _get_pdb_string(as_file.getTopology(), as_file.getPositions())
  if checks:
    _check_cleaned_atoms(pdb_string, prot_pdb_string)
  return pdb_string


def make_atom14_positions(prot):
  """Constructs denser atom positions (14 dimensions instead of 37)."""
  restype_atom14_to_atom37 = []  # mapping (restype, atom14) --> atom37
  restype_atom37_to_atom14 = []  # mapping (restype, atom37) --> atom14
  restype_atom14_mask = []

  for rt in residue_constants.restypes:
    atom_names = residue_constants.restype_name_to_atom14_names[
        residue_constants.restype_1to3[rt]]

    restype_atom14_to_atom37.append([
        (residue_constants.atom_order[name] if name else 0)
        for name in atom_names
    ])

    atom_name_to_idx14 = {name: i for i, name in enumerate(atom_names)}
    restype_atom37_to_atom14.append([
        (atom_name_to_idx14[name] if name in atom_name_to_idx14 else 0)
        for name in residue_constants.atom_types
    ])

    restype_atom14_mask.append([(1. if name else 0.) for name in atom_names])

  # Add dummy mapping for restype 'UNK'.
  restype_atom14_to_atom37.append([0] * 14)
  restype_atom37_to_atom14.append([0] * 37)
  restype_atom14_mask.append([0.] * 14)

  restype_atom14_to_atom37 = np.array(restype_atom14_to_atom37, dtype=np.int32)
  restype_atom37_to_atom14 = np.array(restype_atom37_to_atom14, dtype=np.int32)
  restype_atom14_mask = np.array(restype_atom14_mask, dtype=np.float32)

  # Create the mapping for (residx, atom14) --> atom37, i.e. an array
  # with shape (num_res, 14) containing the atom37 indices for this protein.
  residx_atom14_to_atom37 = restype_atom14_to_atom37[prot["aatype"]]
  residx_atom14_mask = restype_atom14_mask[prot["aatype"]]

  # Create a mask for known ground truth positions.
  residx_atom14_gt_mask = residx_atom14_mask * np.take_along_axis(
      prot["all_atom_mask"], residx_atom14_to_atom37, axis=1).astype(np.float32)

  # Gather the ground truth positions.
  residx_atom14_gt_positions = residx_atom14_gt_mask[:, :, None] * (
      np.take_along_axis(prot["all_atom_positions"],
                         residx_atom14_to_atom37[..., None],
                         axis=1))

  prot["atom14_atom_exists"] = residx_atom14_mask
  prot["atom14_gt_exists"] = residx_atom14_gt_mask
  prot["atom14_gt_positions"] = residx_atom14_gt_positions

  prot["residx_atom14_to_atom37"] = residx_atom14_to_atom37

  # Create the gather indices for mapping back.
  residx_atom37_to_atom14 = restype_atom37_to_atom14[prot["aatype"]]
  prot["residx_atom37_to_atom14"] = residx_atom37_to_atom14

  # Create the corresponding mask.
  restype_atom37_mask = np.zeros([21, 37], dtype=np.float32)
  for restype, restype_letter in enumerate(residue_constants.restypes):
    restype_name = residue_constants.restype_1to3[restype_letter]
    atom_names = residue_constants.residue_atoms[restype_name]
    for atom_name in atom_names:
      atom_type = residue_constants.atom_order[atom_name]
      restype_atom37_mask[restype, atom_type] = 1

  residx_atom37_mask = restype_atom37_mask[prot["aatype"]]
  prot["atom37_atom_exists"] = residx_atom37_mask

  # As the atom naming is ambiguous for 7 of the 20 amino acids, provide
  # alternative ground truth coordinates where the naming is swapped
  restype_3 = [
      residue_constants.restype_1to3[res] for res in residue_constants.restypes
  ]
  restype_3 += ["UNK"]

  # Matrices for renaming ambiguous atoms.
  all_matrices = {res: np.eye(14, dtype=np.float32) for res in restype_3}
  for resname, swap in residue_constants.residue_atom_renaming_swaps.items():
    correspondences = np.arange(14)
    for source_atom_swap, target_atom_swap in swap.items():
      source_index = residue_constants.restype_name_to_atom14_names[
          resname].index(source_atom_swap)
      target_index = residue_constants.restype_name_to_atom14_names[
          resname].index(target_atom_swap)
      correspondences[source_index] = target_index
      correspondences[target_index] = source_index
      renaming_matrix = np.zeros((14, 14), dtype=np.float32)
      for index, correspondence in enumerate(correspondences):
        renaming_matrix[index, correspondence] = 1.
    all_matrices[resname] = renaming_matrix.astype(np.float32)
  renaming_matrices = np.stack([all_matrices[restype] for restype in restype_3])

  # Pick the transformation matrices for the given residue sequence
  # shape (num_res, 14, 14).
  renaming_transform = renaming_matrices[prot["aatype"]]

  # Apply it to the ground truth positions. shape (num_res, 14, 3).
  alternative_gt_positions = np.einsum("rac,rab->rbc",
                                       residx_atom14_gt_positions,
                                       renaming_transform)
  prot["atom14_alt_gt_positions"] = alternative_gt_positions

  # Create the mask for the alternative ground truth (differs from the
  # ground truth mask, if only one of the atoms in an ambiguous pair has a
  # ground truth position).
  alternative_gt_mask = np.einsum("ra,rab->rb",
                                  residx_atom14_gt_mask,
                                  renaming_transform)

  prot["atom14_alt_gt_exists"] = alternative_gt_mask

  # Create an ambiguous atoms mask.  shape: (21, 14).
  restype_atom14_is_ambiguous = np.zeros((21, 14), dtype=np.float32)
  for resname, swap in residue_constants.residue_atom_renaming_swaps.items():
    for atom_name1, atom_name2 in swap.items():
      restype = residue_constants.restype_order[
          residue_constants.restype_3to1[resname]]
      atom_idx1 = residue_constants.restype_name_to_atom14_names[resname].index(
          atom_name1)
      atom_idx2 = residue_constants.restype_name_to_atom14_names[resname].index(
          atom_name2)
      restype_atom14_is_ambiguous[restype, atom_idx1] = 1
      restype_atom14_is_ambiguous[restype, atom_idx2] = 1

  # From this create an ambiguous_mask for the given sequence.
  prot["atom14_atom_is_ambiguous"] = (
      restype_atom14_is_ambiguous[prot["aatype"]])

  return prot


def find_violations(prot_np: protein.Protein):
  """Analyzes a protein and returns structural violation information.

  Args:
    prot_np: A protein.

  Returns:
    violations: A `dict` of structure components with structural violations.
    violation_metrics: A `dict` of violation metrics.
  """
  batch = {
      "aatype": prot_np.aatype,
      "all_atom_positions": prot_np.atom_positions.astype(np.float32),
      "all_atom_mask": prot_np.atom_mask.astype(np.float32),
      "residue_index": prot_np.residue_index,
  }

  batch["seq_mask"] = np.ones_like(batch["aatype"], np.float32)
  batch = make_atom14_positions(batch)

  violations = folding.find_structural_violations(
      batch=batch,
      atom14_pred_positions=batch["atom14_gt_positions"],
      config=ml_collections.ConfigDict(
          {"violation_tolerance_factor": 12,  # Taken from model config.
           "clash_overlap_tolerance": 1.5,  # Taken from model config.
          }))
  violation_metrics = folding.compute_violation_metrics(
      batch=batch,
      atom14_pred_positions=batch["atom14_gt_positions"],
      violations=violations,
  )

  return violations, violation_metrics


def get_violation_metrics(prot: protein.Protein):
  """Computes violation and alignment metrics."""
  structural_violations, struct_metrics = find_violations(prot)
  violation_idx = np.flatnonzero(
      structural_violations["total_per_residue_violations_mask"])

  struct_metrics["residue_violations"] = violation_idx
  struct_metrics["num_residue_violations"] = len(violation_idx)
  struct_metrics["structural_violations"] = structural_violations
  return struct_metrics


def _run_one_iteration(
    *,
    pdb_string: str,
    max_iterations: int,
    tolerance: float,
    stiffness: float,
    restraint_set: str,
    max_attempts: int,
    use_gpu: bool,
    exclude_residues: Optional[Collection[int]] = None):
  """Runs the minimization pipeline.

  Args:
    pdb_string: A pdb string.
    max_iterations: An `int` specifying the maximum number of L-BFGS iterations.
    A value of 0 specifies no limit.
    tolerance: kcal/mol, the energy tolerance of L-BFGS.
    stiffness: kcal/mol A**2, spring constant of heavy atom restraining
      potential.
    restraint_set: The set of atoms to restrain.
    max_attempts: The maximum number of minimization attempts.
    use_gpu: Whether to run on GPU.
    exclude_residues: An optional list of zero-indexed residues to exclude from
        restraints.

  Returns:
    A `dict` of minimization info.
  """
  exclude_residues = exclude_residues or []

  # Assign physical dimensions.
  tolerance = tolerance * ENERGY
  stiffness = stiffness * ENERGY / (LENGTH**2)

  start = time.time()
  minimized = False
  attempts = 0
  while not minimized and attempts < max_attempts:
    attempts += 1
    try:
      logging.info("Minimizing protein, attempt %d of %d.",
                   attempts, max_attempts)
      ret = _openmm_minimize(
          pdb_string, max_iterations=max_iterations,
          tolerance=tolerance, stiffness=stiffness,
          restraint_set=restraint_set,
          exclude_residues=exclude_residues,
          use_gpu=use_gpu)
      minimized = True
    except Exception as e:  # pylint: disable=broad-except
      logging.info(e)
  if not minimized:
    raise ValueError(f"Minimization failed after {max_attempts} attempts.")
  ret["opt_time"] = time.time() - start
  ret["min_attempts"] = attempts
  return ret


def run_pipeline(
    prot: protein.Protein,
    stiffness: float,
    use_gpu: bool,
    max_outer_iterations: int = 1,
    place_hydrogens_every_iteration: bool = True,
    max_iterations: int = 0,
    tolerance: float = 2.39,
    restraint_set: str = "non_hydrogen",
    max_attempts: int = 100,
    checks: bool = True,
    exclude_residues: Optional[Sequence[int]] = None):
  """Run iterative amber relax.

  Successive relax iterations are performed until all violations have been
  resolved. Each iteration involves a restrained Amber minimization, with
  restraint exclusions determined by violation-participating residues.

  Args:
    prot: A protein to be relaxed.
    stiffness: kcal/mol A**2, the restraint stiffness.
    use_gpu: Whether to run on GPU.
    max_outer_iterations: The maximum number of iterative minimization.
    place_hydrogens_every_iteration: Whether hydrogens are re-initialized
        prior to every minimization.
    max_iterations: An `int` specifying the maximum number of L-BFGS steps
        per relax iteration. A value of 0 specifies no limit.
    tolerance: kcal/mol, the energy tolerance of L-BFGS.
        The default value is the OpenMM default.
    restraint_set: The set of atoms to restrain.
    max_attempts: The maximum number of minimization attempts per iteration.
    checks: Whether to perform cleaning checks.
    exclude_residues: An optional list of zero-indexed residues to exclude from
        restraints.

  Returns:
    out: A dictionary of output values.
  """

  # `protein.to_pdb` will strip any poorly-defined residues so we need to
  # perform this check before `clean_protein`.
  _check_residues_are_well_defined(prot)
  pdb_string = clean_protein(prot, checks=checks)

  exclude_residues = exclude_residues or []
  exclude_residues = set(exclude_residues)
  violations = np.inf
  iteration = 0

  while violations > 0 and iteration < max_outer_iterations:
    ret = _run_one_iteration(
        pdb_string=pdb_string,
        exclude_residues=exclude_residues,
        max_iterations=max_iterations,
        tolerance=tolerance,
        stiffness=stiffness,
        restraint_set=restraint_set,
        max_attempts=max_attempts,
        use_gpu=use_gpu)
    prot = protein.from_pdb_string(ret["min_pdb"])
    if place_hydrogens_every_iteration:
      pdb_string = clean_protein(prot, checks=True)
    else:
      pdb_string = ret["min_pdb"]
    # Calculation of violations can cause CUDA errors for some JAX versions.
    with jax.default_device(jax.local_devices(backend="cpu")[0]):
      ret.update(get_violation_metrics(prot))
    ret.update({
        "num_exclusions": len(exclude_residues),
        "iteration": iteration,
    })
    violations = ret["violations_per_residue"]
    exclude_residues = exclude_residues.union(ret["residue_violations"])

    logging.info("Iteration completed: Einit %.2f Efinal %.2f Time %.2f s "
                 "num residue violations %d num residue exclusions %d ",
                 ret["einit"], ret["efinal"], ret["opt_time"],
                 ret["num_residue_violations"], ret["num_exclusions"])
    iteration += 1
  return ret
