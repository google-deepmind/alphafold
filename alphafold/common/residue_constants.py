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

"""Constants used in AlphaFold."""

import collections
import functools
import os
from typing import Final, List, Mapping, Tuple

import numpy as np
import tree

# Internal import (35fd).


# Distance from one CA to next CA [trans configuration: omega = 180].
ca_ca = 3.80209737096

# Format: The list for each AA type contains chi1, chi2, chi3, chi4 in
# this order (or a relevant subset from chi1 onwards). ALA and GLY don't have
# chi angles so their chi angle lists are empty.
chi_angles_atoms = {
    'ALA': [],
    # Chi5 in arginine is always 0 +- 5 degrees, so ignore it.
    'ARG': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD'],
            ['CB', 'CG', 'CD', 'NE'], ['CG', 'CD', 'NE', 'CZ']],
    'ASN': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'OD1']],
    'ASP': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'OD1']],
    'CYS': [['N', 'CA', 'CB', 'SG']],
    'GLN': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD'],
            ['CB', 'CG', 'CD', 'OE1']],
    'GLU': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD'],
            ['CB', 'CG', 'CD', 'OE1']],
    'GLY': [],
    'HIS': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'ND1']],
    'ILE': [['N', 'CA', 'CB', 'CG1'], ['CA', 'CB', 'CG1', 'CD1']],
    'LEU': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD1']],
    'LYS': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD'],
            ['CB', 'CG', 'CD', 'CE'], ['CG', 'CD', 'CE', 'NZ']],
    'MET': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'SD'],
            ['CB', 'CG', 'SD', 'CE']],
    'PHE': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD1']],
    'PRO': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD']],
    'SER': [['N', 'CA', 'CB', 'OG']],
    'THR': [['N', 'CA', 'CB', 'OG1']],
    'TRP': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD1']],
    'TYR': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD1']],
    'VAL': [['N', 'CA', 'CB', 'CG1']],
}

# If chi angles given in fixed-length array, this matrix determines how to mask
# them for each AA type. The order is as per restype_order (see below).
chi_angles_mask = [
    [0.0, 0.0, 0.0, 0.0],  # ALA
    [1.0, 1.0, 1.0, 1.0],  # ARG
    [1.0, 1.0, 0.0, 0.0],  # ASN
    [1.0, 1.0, 0.0, 0.0],  # ASP
    [1.0, 0.0, 0.0, 0.0],  # CYS
    [1.0, 1.0, 1.0, 0.0],  # GLN
    [1.0, 1.0, 1.0, 0.0],  # GLU
    [0.0, 0.0, 0.0, 0.0],  # GLY
    [1.0, 1.0, 0.0, 0.0],  # HIS
    [1.0, 1.0, 0.0, 0.0],  # ILE
    [1.0, 1.0, 0.0, 0.0],  # LEU
    [1.0, 1.0, 1.0, 1.0],  # LYS
    [1.0, 1.0, 1.0, 0.0],  # MET
    [1.0, 1.0, 0.0, 0.0],  # PHE
    [1.0, 1.0, 0.0, 0.0],  # PRO
    [1.0, 0.0, 0.0, 0.0],  # SER
    [1.0, 0.0, 0.0, 0.0],  # THR
    [1.0, 1.0, 0.0, 0.0],  # TRP
    [1.0, 1.0, 0.0, 0.0],  # TYR
    [1.0, 0.0, 0.0, 0.0],  # VAL
]

# The following chi angles are pi periodic: they can be rotated by a multiple
# of pi without affecting the structure.
chi_pi_periodic = [
    [0.0, 0.0, 0.0, 0.0],  # ALA
    [0.0, 0.0, 0.0, 0.0],  # ARG
    [0.0, 0.0, 0.0, 0.0],  # ASN
    [0.0, 1.0, 0.0, 0.0],  # ASP
    [0.0, 0.0, 0.0, 0.0],  # CYS
    [0.0, 0.0, 0.0, 0.0],  # GLN
    [0.0, 0.0, 1.0, 0.0],  # GLU
    [0.0, 0.0, 0.0, 0.0],  # GLY
    [0.0, 0.0, 0.0, 0.0],  # HIS
    [0.0, 0.0, 0.0, 0.0],  # ILE
    [0.0, 0.0, 0.0, 0.0],  # LEU
    [0.0, 0.0, 0.0, 0.0],  # LYS
    [0.0, 0.0, 0.0, 0.0],  # MET
    [0.0, 1.0, 0.0, 0.0],  # PHE
    [0.0, 0.0, 0.0, 0.0],  # PRO
    [0.0, 0.0, 0.0, 0.0],  # SER
    [0.0, 0.0, 0.0, 0.0],  # THR
    [0.0, 0.0, 0.0, 0.0],  # TRP
    [0.0, 1.0, 0.0, 0.0],  # TYR
    [0.0, 0.0, 0.0, 0.0],  # VAL
    [0.0, 0.0, 0.0, 0.0],  # UNK
]

# Atoms positions relative to the 8 rigid groups, defined by the pre-omega, phi,
# psi and chi angles:
# 0: 'backbone group',
# 1: 'pre-omega-group', (empty)
# 2: 'phi-group', (currently empty, because it defines only hydrogens)
# 3: 'psi-group',
# 4,5,6,7: 'chi1,2,3,4-group'
# The atom positions are relative to the axis-end-atom of the corresponding
# rotation axis. The x-axis is in direction of the rotation axis, and the y-axis
# is defined such that the dihedral-angle-defining atom (the last entry in
# chi_angles_atoms above) is in the xy-plane (with a positive y-coordinate).
# format: [atomname, group_idx, rel_position]
rigid_group_atom_positions = {
    'ALA': [
        ['N', 0, (-0.525, 1.363, 0.000)],
        ['CA', 0, (0.000, 0.000, 0.000)],
        ['C', 0, (1.526, -0.000, -0.000)],
        ['CB', 0, (-0.529, -0.774, -1.205)],
        ['O', 3, (0.627, 1.062, 0.000)],
    ],
    'ARG': [
        ['N', 0, (-0.524, 1.362, -0.000)],
        ['CA', 0, (0.000, 0.000, 0.000)],
        ['C', 0, (1.525, -0.000, -0.000)],
        ['CB', 0, (-0.524, -0.778, -1.209)],
        ['O', 3, (0.626, 1.062, 0.000)],
        ['CG', 4, (0.616, 1.390, -0.000)],
        ['CD', 5, (0.564, 1.414, 0.000)],
        ['NE', 6, (0.539, 1.357, -0.000)],
        ['NH1', 7, (0.206, 2.301, 0.000)],
        ['NH2', 7, (2.078, 0.978, -0.000)],
        ['CZ', 7, (0.758, 1.093, -0.000)],
    ],
    'ASN': [
        ['N', 0, (-0.536, 1.357, 0.000)],
        ['CA', 0, (0.000, 0.000, 0.000)],
        ['C', 0, (1.526, -0.000, -0.000)],
        ['CB', 0, (-0.531, -0.787, -1.200)],
        ['O', 3, (0.625, 1.062, 0.000)],
        ['CG', 4, (0.584, 1.399, 0.000)],
        ['ND2', 5, (0.593, -1.188, 0.001)],
        ['OD1', 5, (0.633, 1.059, 0.000)],
    ],
    'ASP': [
        ['N', 0, (-0.525, 1.362, -0.000)],
        ['CA', 0, (0.000, 0.000, 0.000)],
        ['C', 0, (1.527, 0.000, -0.000)],
        ['CB', 0, (-0.526, -0.778, -1.208)],
        ['O', 3, (0.626, 1.062, -0.000)],
        ['CG', 4, (0.593, 1.398, -0.000)],
        ['OD1', 5, (0.610, 1.091, 0.000)],
        ['OD2', 5, (0.592, -1.101, -0.003)],
    ],
    'CYS': [
        ['N', 0, (-0.522, 1.362, -0.000)],
        ['CA', 0, (0.000, 0.000, 0.000)],
        ['C', 0, (1.524, 0.000, 0.000)],
        ['CB', 0, (-0.519, -0.773, -1.212)],
        ['O', 3, (0.625, 1.062, -0.000)],
        ['SG', 4, (0.728, 1.653, 0.000)],
    ],
    'GLN': [
        ['N', 0, (-0.526, 1.361, -0.000)],
        ['CA', 0, (0.000, 0.000, 0.000)],
        ['C', 0, (1.526, 0.000, 0.000)],
        ['CB', 0, (-0.525, -0.779, -1.207)],
        ['O', 3, (0.626, 1.062, -0.000)],
        ['CG', 4, (0.615, 1.393, 0.000)],
        ['CD', 5, (0.587, 1.399, -0.000)],
        ['NE2', 6, (0.593, -1.189, -0.001)],
        ['OE1', 6, (0.634, 1.060, 0.000)],
    ],
    'GLU': [
        ['N', 0, (-0.528, 1.361, 0.000)],
        ['CA', 0, (0.000, 0.000, 0.000)],
        ['C', 0, (1.526, -0.000, -0.000)],
        ['CB', 0, (-0.526, -0.781, -1.207)],
        ['O', 3, (0.626, 1.062, 0.000)],
        ['CG', 4, (0.615, 1.392, 0.000)],
        ['CD', 5, (0.600, 1.397, 0.000)],
        ['OE1', 6, (0.607, 1.095, -0.000)],
        ['OE2', 6, (0.589, -1.104, -0.001)],
    ],
    'GLY': [
        ['N', 0, (-0.572, 1.337, 0.000)],
        ['CA', 0, (0.000, 0.000, 0.000)],
        ['C', 0, (1.517, -0.000, -0.000)],
        ['O', 3, (0.626, 1.062, -0.000)],
    ],
    'HIS': [
        ['N', 0, (-0.527, 1.360, 0.000)],
        ['CA', 0, (0.000, 0.000, 0.000)],
        ['C', 0, (1.525, 0.000, 0.000)],
        ['CB', 0, (-0.525, -0.778, -1.208)],
        ['O', 3, (0.625, 1.063, 0.000)],
        ['CG', 4, (0.600, 1.370, -0.000)],
        ['CD2', 5, (0.889, -1.021, 0.003)],
        ['ND1', 5, (0.744, 1.160, -0.000)],
        ['CE1', 5, (2.030, 0.851, 0.002)],
        ['NE2', 5, (2.145, -0.466, 0.004)],
    ],
    'ILE': [
        ['N', 0, (-0.493, 1.373, -0.000)],
        ['CA', 0, (0.000, 0.000, 0.000)],
        ['C', 0, (1.527, -0.000, -0.000)],
        ['CB', 0, (-0.536, -0.793, -1.213)],
        ['O', 3, (0.627, 1.062, -0.000)],
        ['CG1', 4, (0.534, 1.437, -0.000)],
        ['CG2', 4, (0.540, -0.785, -1.199)],
        ['CD1', 5, (0.619, 1.391, 0.000)],
    ],
    'LEU': [
        ['N', 0, (-0.520, 1.363, 0.000)],
        ['CA', 0, (0.000, 0.000, 0.000)],
        ['C', 0, (1.525, -0.000, -0.000)],
        ['CB', 0, (-0.522, -0.773, -1.214)],
        ['O', 3, (0.625, 1.063, -0.000)],
        ['CG', 4, (0.678, 1.371, 0.000)],
        ['CD1', 5, (0.530, 1.430, -0.000)],
        ['CD2', 5, (0.535, -0.774, 1.200)],
    ],
    'LYS': [
        ['N', 0, (-0.526, 1.362, -0.000)],
        ['CA', 0, (0.000, 0.000, 0.000)],
        ['C', 0, (1.526, 0.000, 0.000)],
        ['CB', 0, (-0.524, -0.778, -1.208)],
        ['O', 3, (0.626, 1.062, -0.000)],
        ['CG', 4, (0.619, 1.390, 0.000)],
        ['CD', 5, (0.559, 1.417, 0.000)],
        ['CE', 6, (0.560, 1.416, 0.000)],
        ['NZ', 7, (0.554, 1.387, 0.000)],
    ],
    'MET': [
        ['N', 0, (-0.521, 1.364, -0.000)],
        ['CA', 0, (0.000, 0.000, 0.000)],
        ['C', 0, (1.525, 0.000, 0.000)],
        ['CB', 0, (-0.523, -0.776, -1.210)],
        ['O', 3, (0.625, 1.062, -0.000)],
        ['CG', 4, (0.613, 1.391, -0.000)],
        ['SD', 5, (0.703, 1.695, 0.000)],
        ['CE', 6, (0.320, 1.786, -0.000)],
    ],
    'PHE': [
        ['N', 0, (-0.518, 1.363, 0.000)],
        ['CA', 0, (0.000, 0.000, 0.000)],
        ['C', 0, (1.524, 0.000, -0.000)],
        ['CB', 0, (-0.525, -0.776, -1.212)],
        ['O', 3, (0.626, 1.062, -0.000)],
        ['CG', 4, (0.607, 1.377, 0.000)],
        ['CD1', 5, (0.709, 1.195, -0.000)],
        ['CD2', 5, (0.706, -1.196, 0.000)],
        ['CE1', 5, (2.102, 1.198, -0.000)],
        ['CE2', 5, (2.098, -1.201, -0.000)],
        ['CZ', 5, (2.794, -0.003, -0.001)],
    ],
    'PRO': [
        ['N', 0, (-0.566, 1.351, -0.000)],
        ['CA', 0, (0.000, 0.000, 0.000)],
        ['C', 0, (1.527, -0.000, 0.000)],
        ['CB', 0, (-0.546, -0.611, -1.293)],
        ['O', 3, (0.621, 1.066, 0.000)],
        ['CG', 4, (0.382, 1.445, 0.0)],
        # ['CD', 5, (0.427, 1.440, 0.0)],
        ['CD', 5, (0.477, 1.424, 0.0)],  # manually made angle 2 degrees larger
    ],
    'SER': [
        ['N', 0, (-0.529, 1.360, -0.000)],
        ['CA', 0, (0.000, 0.000, 0.000)],
        ['C', 0, (1.525, -0.000, -0.000)],
        ['CB', 0, (-0.518, -0.777, -1.211)],
        ['O', 3, (0.626, 1.062, -0.000)],
        ['OG', 4, (0.503, 1.325, 0.000)],
    ],
    'THR': [
        ['N', 0, (-0.517, 1.364, 0.000)],
        ['CA', 0, (0.000, 0.000, 0.000)],
        ['C', 0, (1.526, 0.000, -0.000)],
        ['CB', 0, (-0.516, -0.793, -1.215)],
        ['O', 3, (0.626, 1.062, 0.000)],
        ['CG2', 4, (0.550, -0.718, -1.228)],
        ['OG1', 4, (0.472, 1.353, 0.000)],
    ],
    'TRP': [
        ['N', 0, (-0.521, 1.363, 0.000)],
        ['CA', 0, (0.000, 0.000, 0.000)],
        ['C', 0, (1.525, -0.000, 0.000)],
        ['CB', 0, (-0.523, -0.776, -1.212)],
        ['O', 3, (0.627, 1.062, 0.000)],
        ['CG', 4, (0.609, 1.370, -0.000)],
        ['CD1', 5, (0.824, 1.091, 0.000)],
        ['CD2', 5, (0.854, -1.148, -0.005)],
        ['CE2', 5, (2.186, -0.678, -0.007)],
        ['CE3', 5, (0.622, -2.530, -0.007)],
        ['NE1', 5, (2.140, 0.690, -0.004)],
        ['CH2', 5, (3.028, -2.890, -0.013)],
        ['CZ2', 5, (3.283, -1.543, -0.011)],
        ['CZ3', 5, (1.715, -3.389, -0.011)],
    ],
    'TYR': [
        ['N', 0, (-0.522, 1.362, 0.000)],
        ['CA', 0, (0.000, 0.000, 0.000)],
        ['C', 0, (1.524, -0.000, -0.000)],
        ['CB', 0, (-0.522, -0.776, -1.213)],
        ['O', 3, (0.627, 1.062, -0.000)],
        ['CG', 4, (0.607, 1.382, -0.000)],
        ['CD1', 5, (0.716, 1.195, -0.000)],
        ['CD2', 5, (0.713, -1.194, -0.001)],
        ['CE1', 5, (2.107, 1.200, -0.002)],
        ['CE2', 5, (2.104, -1.201, -0.003)],
        ['OH', 5, (4.168, -0.002, -0.005)],
        ['CZ', 5, (2.791, -0.001, -0.003)],
    ],
    'VAL': [
        ['N', 0, (-0.494, 1.373, -0.000)],
        ['CA', 0, (0.000, 0.000, 0.000)],
        ['C', 0, (1.527, -0.000, -0.000)],
        ['CB', 0, (-0.533, -0.795, -1.213)],
        ['O', 3, (0.627, 1.062, -0.000)],
        ['CG1', 4, (0.540, 1.429, -0.000)],
        ['CG2', 4, (0.533, -0.776, 1.203)],
    ],
}

# A list of atoms (excluding hydrogen) for each AA type. PDB naming convention.
residue_atoms = {
    'ALA': ['C', 'CA', 'CB', 'N', 'O'],
    'ARG': ['C', 'CA', 'CB', 'CG', 'CD', 'CZ', 'N', 'NE', 'O', 'NH1', 'NH2'],
    'ASP': ['C', 'CA', 'CB', 'CG', 'N', 'O', 'OD1', 'OD2'],
    'ASN': ['C', 'CA', 'CB', 'CG', 'N', 'ND2', 'O', 'OD1'],
    'CYS': ['C', 'CA', 'CB', 'N', 'O', 'SG'],
    'GLU': ['C', 'CA', 'CB', 'CG', 'CD', 'N', 'O', 'OE1', 'OE2'],
    'GLN': ['C', 'CA', 'CB', 'CG', 'CD', 'N', 'NE2', 'O', 'OE1'],
    'GLY': ['C', 'CA', 'N', 'O'],
    'HIS': ['C', 'CA', 'CB', 'CG', 'CD2', 'CE1', 'N', 'ND1', 'NE2', 'O'],
    'ILE': ['C', 'CA', 'CB', 'CG1', 'CG2', 'CD1', 'N', 'O'],
    'LEU': ['C', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'N', 'O'],
    'LYS': ['C', 'CA', 'CB', 'CG', 'CD', 'CE', 'N', 'NZ', 'O'],
    'MET': ['C', 'CA', 'CB', 'CG', 'CE', 'N', 'O', 'SD'],
    'PHE': ['C', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'N', 'O'],
    'PRO': ['C', 'CA', 'CB', 'CG', 'CD', 'N', 'O'],
    'SER': ['C', 'CA', 'CB', 'N', 'O', 'OG'],
    'THR': ['C', 'CA', 'CB', 'CG2', 'N', 'O', 'OG1'],
    'TRP': ['C', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3',
            'CH2', 'N', 'NE1', 'O'],
    'TYR': ['C', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'N', 'O',
            'OH'],
    'VAL': ['C', 'CA', 'CB', 'CG1', 'CG2', 'N', 'O']
}

# Naming swaps for ambiguous atom names.
# Due to symmetries in the amino acids the naming of atoms is ambiguous in
# 4 of the 20 amino acids.
# (The LDDT paper lists 7 amino acids as ambiguous, but the naming ambiguities
# in LEU, VAL and ARG can be resolved by using the 3d constellations of
# the 'ambiguous' atoms and their neighbours)
residue_atom_renaming_swaps = {
    'ASP': {'OD1': 'OD2'},
    'GLU': {'OE1': 'OE2'},
    'PHE': {'CD1': 'CD2', 'CE1': 'CE2'},
    'TYR': {'CD1': 'CD2', 'CE1': 'CE2'},
}

# Van der Waals radii [Angstroem] of the atoms (from Wikipedia)
van_der_waals_radius = {
    'C': 1.7,
    'N': 1.55,
    'O': 1.52,
    'S': 1.8,
}

Bond = collections.namedtuple(
    'Bond', ['atom1_name', 'atom2_name', 'length', 'stddev'])
BondAngle = collections.namedtuple(
    'BondAngle',
    ['atom1_name', 'atom2_name', 'atom3name', 'angle_rad', 'stddev'])


@functools.lru_cache(maxsize=None)
def load_stereo_chemical_props() -> Tuple[Mapping[str, List[Bond]],
                                          Mapping[str, List[Bond]],
                                          Mapping[str, List[BondAngle]]]:
  """Load stereo_chemical_props.txt into a nice structure.

  Load literature values for bond lengths and bond angles and translate
  bond angles into the length of the opposite edge of the triangle
  ("residue_virtual_bonds").

  Returns:
    residue_bonds: Dict that maps resname -> list of Bond tuples.
    residue_virtual_bonds: Dict that maps resname -> list of Bond tuples.
    residue_bond_angles: Dict that maps resname -> list of BondAngle tuples.
  """
  stereo_chemical_props_path = os.path.join(
      os.path.dirname(os.path.abspath(__file__)), 'stereo_chemical_props.txt'
  )
  with open(stereo_chemical_props_path, 'rt') as f:
    stereo_chemical_props = f.read()
  lines_iter = iter(stereo_chemical_props.splitlines())
  # Load bond lengths.
  residue_bonds = {}
  next(lines_iter)  # Skip header line.
  for line in lines_iter:
    if line.strip() == '-':
      break
    bond, resname, length, stddev = line.split()
    atom1, atom2 = bond.split('-')
    if resname not in residue_bonds:
      residue_bonds[resname] = []
    residue_bonds[resname].append(
        Bond(atom1, atom2, float(length), float(stddev)))
  residue_bonds['UNK'] = []

  # Load bond angles.
  residue_bond_angles = {}
  next(lines_iter)  # Skip empty line.
  next(lines_iter)  # Skip header line.
  for line in lines_iter:
    if line.strip() == '-':
      break
    bond, resname, angle_degree, stddev_degree = line.split()
    atom1, atom2, atom3 = bond.split('-')
    if resname not in residue_bond_angles:
      residue_bond_angles[resname] = []
    residue_bond_angles[resname].append(
        BondAngle(atom1, atom2, atom3,
                  float(angle_degree) / 180. * np.pi,
                  float(stddev_degree) / 180. * np.pi))
  residue_bond_angles['UNK'] = []

  def make_bond_key(atom1_name, atom2_name):
    """Unique key to lookup bonds."""
    return '-'.join(sorted([atom1_name, atom2_name]))

  # Translate bond angles into distances ("virtual bonds").
  residue_virtual_bonds = {}
  for resname, bond_angles in residue_bond_angles.items():
    # Create a fast lookup dict for bond lengths.
    bond_cache = {}
    for b in residue_bonds[resname]:
      bond_cache[make_bond_key(b.atom1_name, b.atom2_name)] = b
    residue_virtual_bonds[resname] = []
    for ba in bond_angles:
      bond1 = bond_cache[make_bond_key(ba.atom1_name, ba.atom2_name)]
      bond2 = bond_cache[make_bond_key(ba.atom2_name, ba.atom3name)]

      # Compute distance between atom1 and atom3 using the law of cosines
      # c^2 = a^2 + b^2 - 2ab*cos(gamma).
      gamma = ba.angle_rad
      length = np.sqrt(bond1.length**2 + bond2.length**2
                       - 2 * bond1.length * bond2.length * np.cos(gamma))

      # Propagation of uncertainty assuming uncorrelated errors.
      dl_outer = 0.5 / length
      dl_dgamma = (2 * bond1.length * bond2.length * np.sin(gamma)) * dl_outer
      dl_db1 = (2 * bond1.length - 2 * bond2.length * np.cos(gamma)) * dl_outer
      dl_db2 = (2 * bond2.length - 2 * bond1.length * np.cos(gamma)) * dl_outer
      stddev = np.sqrt((dl_dgamma * ba.stddev)**2 +
                       (dl_db1 * bond1.stddev)**2 +
                       (dl_db2 * bond2.stddev)**2)
      residue_virtual_bonds[resname].append(
          Bond(ba.atom1_name, ba.atom3name, length, stddev))

  return (residue_bonds,
          residue_virtual_bonds,
          residue_bond_angles)


# Between-residue bond lengths for general bonds (first element) and for Proline
# (second element).
between_res_bond_length_c_n = [1.329, 1.341]
between_res_bond_length_stddev_c_n = [0.014, 0.016]

# Between-residue cos_angles.
between_res_cos_angles_c_n_ca = [-0.5203, 0.0353]  # degrees: 121.352 +- 2.315
between_res_cos_angles_ca_c_n = [-0.4473, 0.0311]  # degrees: 116.568 +- 1.995

# This mapping is used when we need to store atom data in a format that requires
# fixed atom data size for every residue (e.g. a numpy array).
atom_types = [
    'N', 'CA', 'C', 'CB', 'O', 'CG', 'CG1', 'CG2', 'OG', 'OG1', 'SG', 'CD',
    'CD1', 'CD2', 'ND1', 'ND2', 'OD1', 'OD2', 'SD', 'CE', 'CE1', 'CE2', 'CE3',
    'NE', 'NE1', 'NE2', 'OE1', 'OE2', 'CH2', 'NH1', 'NH2', 'OH', 'CZ', 'CZ2',
    'CZ3', 'NZ', 'OXT'
]
atom_order = {atom_type: i for i, atom_type in enumerate(atom_types)}
atom_type_num = len(atom_types)  # := 37.


# A compact atom encoding with 14 columns
# pylint: disable=line-too-long
# pylint: disable=bad-whitespace
restype_name_to_atom14_names = {
    'ALA': ['N', 'CA', 'C', 'O', 'CB', '',    '',    '',    '',    '',    '',    '',    '',    ''],
    'ARG': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD',  'NE',  'CZ',  'NH1', 'NH2', '',    '',    ''],
    'ASN': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'OD1', 'ND2', '',    '',    '',    '',    '',    ''],
    'ASP': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'OD1', 'OD2', '',    '',    '',    '',    '',    ''],
    'CYS': ['N', 'CA', 'C', 'O', 'CB', 'SG',  '',    '',    '',    '',    '',    '',    '',    ''],
    'GLN': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD',  'OE1', 'NE2', '',    '',    '',    '',    ''],
    'GLU': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD',  'OE1', 'OE2', '',    '',    '',    '',    ''],
    'GLY': ['N', 'CA', 'C', 'O', '',   '',    '',    '',    '',    '',    '',    '',    '',    ''],
    'HIS': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'ND1', 'CD2', 'CE1', 'NE2', '',    '',    '',    ''],
    'ILE': ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1', '',    '',    '',    '',    '',    ''],
    'LEU': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD1', 'CD2', '',    '',    '',    '',    '',    ''],
    'LYS': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD',  'CE',  'NZ',  '',    '',    '',    '',    ''],
    'MET': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'SD',  'CE',  '',    '',    '',    '',    '',    ''],
    'PHE': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD1', 'CD2', 'CE1', 'CE2', 'CZ',  '',    '',    ''],
    'PRO': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD',  '',    '',    '',    '',    '',    '',    ''],
    'SER': ['N', 'CA', 'C', 'O', 'CB', 'OG',  '',    '',    '',    '',    '',    '',    '',    ''],
    'THR': ['N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2', '',    '',    '',    '',    '',    '',    ''],
    'TRP': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
    'TYR': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD1', 'CD2', 'CE1', 'CE2', 'CZ',  'OH',  '',    ''],
    'VAL': ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', '',    '',    '',    '',    '',    '',    ''],
    'UNK': ['',  '',   '',  '',  '',   '',    '',    '',    '',    '',    '',    '',    '',    ''],

}
# pylint: enable=line-too-long
# pylint: enable=bad-whitespace


# This is the standard residue order when coding AA type as a number.
# Reproduce it by taking 3-letter AA codes and sorting them alphabetically.
restypes = [
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P',
    'S', 'T', 'W', 'Y', 'V'
]
restype_order = {restype: i for i, restype in enumerate(restypes)}
restype_num = len(restypes)  # := 20.
unk_restype_index = restype_num  # Catch-all index for unknown restypes.

restypes_with_x = restypes + ['X']
restype_order_with_x = {restype: i for i, restype in enumerate(restypes_with_x)}


def sequence_to_onehot(
    sequence: str,
    mapping: Mapping[str, int],
    map_unknown_to_x: bool = False) -> np.ndarray:
  """Maps the given sequence into a one-hot encoded matrix.

  Args:
    sequence: An amino acid sequence.
    mapping: A dictionary mapping amino acids to integers.
    map_unknown_to_x: If True, any amino acid that is not in the mapping will be
      mapped to the unknown amino acid 'X'. If the mapping doesn't contain
      amino acid 'X', an error will be thrown. If False, any amino acid not in
      the mapping will throw an error.

  Returns:
    A numpy array of shape (seq_len, num_unique_aas) with one-hot encoding of
    the sequence.

  Raises:
    ValueError: If the mapping doesn't contain values from 0 to
      num_unique_aas - 1 without any gaps.
  """
  num_entries = max(mapping.values()) + 1

  if sorted(set(mapping.values())) != list(range(num_entries)):
    raise ValueError('The mapping must have values from 0 to num_unique_aas-1 '
                     'without any gaps. Got: %s' % sorted(mapping.values()))

  one_hot_arr = np.zeros((len(sequence), num_entries), dtype=np.int32)

  for aa_index, aa_type in enumerate(sequence):
    if map_unknown_to_x:
      if aa_type.isalpha() and aa_type.isupper():
        aa_id = mapping.get(aa_type, mapping['X'])
      else:
        raise ValueError(f'Invalid character in the sequence: {aa_type}')
    else:
      aa_id = mapping[aa_type]
    one_hot_arr[aa_index, aa_id] = 1

  return one_hot_arr


restype_1to3 = {
    'A': 'ALA',
    'R': 'ARG',
    'N': 'ASN',
    'D': 'ASP',
    'C': 'CYS',
    'Q': 'GLN',
    'E': 'GLU',
    'G': 'GLY',
    'H': 'HIS',
    'I': 'ILE',
    'L': 'LEU',
    'K': 'LYS',
    'M': 'MET',
    'F': 'PHE',
    'P': 'PRO',
    'S': 'SER',
    'T': 'THR',
    'W': 'TRP',
    'Y': 'TYR',
    'V': 'VAL',
}

PROTEIN_CHAIN: Final[str] = 'polypeptide(L)'
POLYMER_CHAIN: Final[str] = 'polymer'


def atom_id_to_type(atom_id: str) -> str:
  """Convert atom ID to atom type, works only for standard protein residues.

  Args:
    atom_id: Atom ID to be converted.

  Returns:
    String corresponding to atom type.

  Raises:
    ValueError: If atom ID not recognized.
  """

  if atom_id.startswith('C'):
    return 'C'
  elif atom_id.startswith('N'):
    return 'N'
  elif atom_id.startswith('O'):
    return 'O'
  elif atom_id.startswith('H'):
    return 'H'
  elif atom_id.startswith('S'):
    return 'S'
  raise ValueError('Atom ID not recognized.')


# NB: restype_3to1 differs from Bio.PDB.protein_letters_3to1 by being a simple
# 1-to-1 mapping of 3 letter names to one letter names. The latter contains
# many more, and less common, three letter names as keys and maps many of these
# to the same one letter name (including 'X' and 'U' which we don't use here).
restype_3to1 = {v: k for k, v in restype_1to3.items()}

# Define a restype name for all unknown residues.
unk_restype = 'UNK'

resnames = [restype_1to3[r] for r in restypes] + [unk_restype]
resname_to_idx = {resname: i for i, resname in enumerate(resnames)}


# The mapping here uses hhblits convention, so that B is mapped to D, J and O
# are mapped to X, U is mapped to C, and Z is mapped to E. Other than that the
# remaining 20 amino acids are kept in alphabetical order.
# There are 2 non-amino acid codes, X (representing any amino acid) and
# "-" representing a missing amino acid in an alignment.  The id for these
# codes is put at the end (20 and 21) so that they can easily be ignored if
# desired.
HHBLITS_AA_TO_ID = {
    'A': 0,
    'B': 2,
    'C': 1,
    'D': 2,
    'E': 3,
    'F': 4,
    'G': 5,
    'H': 6,
    'I': 7,
    'J': 20,
    'K': 8,
    'L': 9,
    'M': 10,
    'N': 11,
    'O': 20,
    'P': 12,
    'Q': 13,
    'R': 14,
    'S': 15,
    'T': 16,
    'U': 1,
    'V': 17,
    'W': 18,
    'X': 20,
    'Y': 19,
    'Z': 3,
    '-': 21,
}

# Partial inversion of HHBLITS_AA_TO_ID.
ID_TO_HHBLITS_AA = {
    0: 'A',
    1: 'C',  # Also U.
    2: 'D',  # Also B.
    3: 'E',  # Also Z.
    4: 'F',
    5: 'G',
    6: 'H',
    7: 'I',
    8: 'K',
    9: 'L',
    10: 'M',
    11: 'N',
    12: 'P',
    13: 'Q',
    14: 'R',
    15: 'S',
    16: 'T',
    17: 'V',
    18: 'W',
    19: 'Y',
    20: 'X',  # Includes J and O.
    21: '-',
}

restypes_with_x_and_gap = restypes + ['X', '-']
MAP_HHBLITS_AATYPE_TO_OUR_AATYPE = tuple(
    restypes_with_x_and_gap.index(ID_TO_HHBLITS_AA[i])
    for i in range(len(restypes_with_x_and_gap)))


def _make_standard_atom_mask() -> np.ndarray:
  """Returns [num_res_types, num_atom_types] mask array."""
  # +1 to account for unknown (all 0s).
  mask = np.zeros([restype_num + 1, atom_type_num], dtype=np.int32)
  for restype, restype_letter in enumerate(restypes):
    restype_name = restype_1to3[restype_letter]
    atom_names = residue_atoms[restype_name]
    for atom_name in atom_names:
      atom_type = atom_order[atom_name]
      mask[restype, atom_type] = 1
  return mask


STANDARD_ATOM_MASK = _make_standard_atom_mask()


# A one hot representation for the first and second atoms defining the axis
# of rotation for each chi-angle in each residue.
def chi_angle_atom(atom_index: int) -> np.ndarray:
  """Define chi-angle rigid groups via one-hot representations."""
  chi_angles_index = {}
  one_hots = []

  for k, v in chi_angles_atoms.items():
    indices = [atom_types.index(s[atom_index]) for s in v]
    indices.extend([-1]*(4-len(indices)))
    chi_angles_index[k] = indices

  for r in restypes:
    res3 = restype_1to3[r]
    one_hot = np.eye(atom_type_num)[chi_angles_index[res3]]
    one_hots.append(one_hot)

  one_hots.append(np.zeros([4, atom_type_num]))  # Add zeros for residue `X`.
  one_hot = np.stack(one_hots, axis=0)
  one_hot = np.transpose(one_hot, [0, 2, 1])

  return one_hot

chi_atom_1_one_hot = chi_angle_atom(1)
chi_atom_2_one_hot = chi_angle_atom(2)

# An array like chi_angles_atoms but using indices rather than names.
chi_angles_atom_indices = [chi_angles_atoms[restype_1to3[r]] for r in restypes]
chi_angles_atom_indices = tree.map_structure(
    lambda atom_name: atom_order[atom_name], chi_angles_atom_indices)
chi_angles_atom_indices = np.array([
    chi_atoms + ([[0, 0, 0, 0]] * (4 - len(chi_atoms)))
    for chi_atoms in chi_angles_atom_indices])

# Mapping from (res_name, atom_name) pairs to the atom's chi group index
# and atom index within that group.
chi_groups_for_atom = collections.defaultdict(list)
for res_name, chi_angle_atoms_for_res in chi_angles_atoms.items():
  for chi_group_i, chi_group in enumerate(chi_angle_atoms_for_res):
    for atom_i, atom in enumerate(chi_group):
      chi_groups_for_atom[(res_name, atom)].append((chi_group_i, atom_i))
chi_groups_for_atom = dict(chi_groups_for_atom)


def _make_rigid_transformation_4x4(ex, ey, translation):
  """Create a rigid 4x4 transformation matrix from two axes and transl."""
  # Normalize ex.
  ex_normalized = ex / np.linalg.norm(ex)

  # make ey perpendicular to ex
  ey_normalized = ey - np.dot(ey, ex_normalized) * ex_normalized
  ey_normalized /= np.linalg.norm(ey_normalized)

  # compute ez as cross product
  eznorm = np.cross(ex_normalized, ey_normalized)
  m = np.stack([ex_normalized, ey_normalized, eznorm, translation]).transpose()
  m = np.concatenate([m, [[0., 0., 0., 1.]]], axis=0)
  return m


# create an array with (restype, atomtype) --> rigid_group_idx
# and an array with (restype, atomtype, coord) for the atom positions
# and compute affine transformation matrices (4,4) from one rigid group to the
# previous group
restype_atom37_to_rigid_group = np.zeros([21, 37], dtype=int)
restype_atom37_mask = np.zeros([21, 37], dtype=np.float32)
restype_atom37_rigid_group_positions = np.zeros([21, 37, 3], dtype=np.float32)
restype_atom14_to_rigid_group = np.zeros([21, 14], dtype=int)
restype_atom14_mask = np.zeros([21, 14], dtype=np.float32)
restype_atom14_rigid_group_positions = np.zeros([21, 14, 3], dtype=np.float32)
restype_rigid_group_default_frame = np.zeros([21, 8, 4, 4], dtype=np.float32)


def _make_rigid_group_constants():
  """Fill the arrays above."""
  for restype, restype_letter in enumerate(restypes):
    resname = restype_1to3[restype_letter]
    for atomname, group_idx, atom_position in rigid_group_atom_positions[
        resname]:
      atomtype = atom_order[atomname]
      restype_atom37_to_rigid_group[restype, atomtype] = group_idx
      restype_atom37_mask[restype, atomtype] = 1
      restype_atom37_rigid_group_positions[restype, atomtype, :] = atom_position

      atom14idx = restype_name_to_atom14_names[resname].index(atomname)
      restype_atom14_to_rigid_group[restype, atom14idx] = group_idx
      restype_atom14_mask[restype, atom14idx] = 1
      restype_atom14_rigid_group_positions[restype,
                                           atom14idx, :] = atom_position

  for restype, restype_letter in enumerate(restypes):
    resname = restype_1to3[restype_letter]
    atom_positions = {name: np.array(pos) for name, _, pos
                      in rigid_group_atom_positions[resname]}

    # backbone to backbone is the identity transform
    restype_rigid_group_default_frame[restype, 0, :, :] = np.eye(4)

    # pre-omega-frame to backbone (currently dummy identity matrix)
    restype_rigid_group_default_frame[restype, 1, :, :] = np.eye(4)

    # phi-frame to backbone
    mat = _make_rigid_transformation_4x4(
        ex=atom_positions['N'] - atom_positions['CA'],
        ey=np.array([1., 0., 0.]),
        translation=atom_positions['N'])
    restype_rigid_group_default_frame[restype, 2, :, :] = mat

    # psi-frame to backbone
    mat = _make_rigid_transformation_4x4(
        ex=atom_positions['C'] - atom_positions['CA'],
        ey=atom_positions['CA'] - atom_positions['N'],
        translation=atom_positions['C'])
    restype_rigid_group_default_frame[restype, 3, :, :] = mat

    # chi1-frame to backbone
    if chi_angles_mask[restype][0]:
      base_atom_names = chi_angles_atoms[resname][0]
      base_atom_positions = [atom_positions[name] for name in base_atom_names]
      mat = _make_rigid_transformation_4x4(
          ex=base_atom_positions[2] - base_atom_positions[1],
          ey=base_atom_positions[0] - base_atom_positions[1],
          translation=base_atom_positions[2])
      restype_rigid_group_default_frame[restype, 4, :, :] = mat

    # chi2-frame to chi1-frame
    # chi3-frame to chi2-frame
    # chi4-frame to chi3-frame
    # luckily all rotation axes for the next frame start at (0,0,0) of the
    # previous frame
    for chi_idx in range(1, 4):
      if chi_angles_mask[restype][chi_idx]:
        axis_end_atom_name = chi_angles_atoms[resname][chi_idx][2]
        axis_end_atom_position = atom_positions[axis_end_atom_name]
        mat = _make_rigid_transformation_4x4(
            ex=axis_end_atom_position,
            ey=np.array([-1., 0., 0.]),
            translation=axis_end_atom_position)
        restype_rigid_group_default_frame[restype, 4 + chi_idx, :, :] = mat


_make_rigid_group_constants()


def make_atom14_dists_bounds(overlap_tolerance=1.5,
                             bond_length_tolerance_factor=15):
  """compute upper and lower bounds for bonds to assess violations."""
  restype_atom14_bond_lower_bound = np.zeros([21, 14, 14], np.float32)
  restype_atom14_bond_upper_bound = np.zeros([21, 14, 14], np.float32)
  restype_atom14_bond_stddev = np.zeros([21, 14, 14], np.float32)
  residue_bonds, residue_virtual_bonds, _ = load_stereo_chemical_props()
  for restype, restype_letter in enumerate(restypes):
    resname = restype_1to3[restype_letter]
    atom_list = restype_name_to_atom14_names[resname]

    # create lower and upper bounds for clashes
    for atom1_idx, atom1_name in enumerate(atom_list):
      if not atom1_name:
        continue
      atom1_radius = van_der_waals_radius[atom1_name[0]]
      for atom2_idx, atom2_name in enumerate(atom_list):
        if (not atom2_name) or atom1_idx == atom2_idx:
          continue
        atom2_radius = van_der_waals_radius[atom2_name[0]]
        lower = atom1_radius + atom2_radius - overlap_tolerance
        upper = 1e10
        restype_atom14_bond_lower_bound[restype, atom1_idx, atom2_idx] = lower
        restype_atom14_bond_lower_bound[restype, atom2_idx, atom1_idx] = lower
        restype_atom14_bond_upper_bound[restype, atom1_idx, atom2_idx] = upper
        restype_atom14_bond_upper_bound[restype, atom2_idx, atom1_idx] = upper

    # overwrite lower and upper bounds for bonds and angles
    for b in residue_bonds[resname] + residue_virtual_bonds[resname]:
      atom1_idx = atom_list.index(b.atom1_name)
      atom2_idx = atom_list.index(b.atom2_name)
      lower = b.length - bond_length_tolerance_factor * b.stddev
      upper = b.length + bond_length_tolerance_factor * b.stddev
      restype_atom14_bond_lower_bound[restype, atom1_idx, atom2_idx] = lower
      restype_atom14_bond_lower_bound[restype, atom2_idx, atom1_idx] = lower
      restype_atom14_bond_upper_bound[restype, atom1_idx, atom2_idx] = upper
      restype_atom14_bond_upper_bound[restype, atom2_idx, atom1_idx] = upper
      restype_atom14_bond_stddev[restype, atom1_idx, atom2_idx] = b.stddev
      restype_atom14_bond_stddev[restype, atom2_idx, atom1_idx] = b.stddev
  return {'lower_bound': restype_atom14_bond_lower_bound,  # shape (21,14,14)
          'upper_bound': restype_atom14_bond_upper_bound,  # shape (21,14,14)
          'stddev': restype_atom14_bond_stddev,  # shape (21,14,14)
         }


CCD_NAME_TO_ONE_LETTER: Mapping[str, str] = {
    '00C': 'C', '01W': 'X', '02K': 'A', '03Y': 'C', '07O': 'C', '08P': 'C',
    '0A0': 'D', '0A1': 'Y', '0A2': 'K', '0A8': 'C', '0AA': 'V', '0AB': 'V',
    '0AC': 'G', '0AD': 'G', '0AF': 'W', '0AG': 'L', '0AH': 'S', '0AK': 'D',
    '0AM': 'A', '0AP': 'C', '0AU': 'U', '0AV': 'A', '0AZ': 'P', '0BN': 'F',
    '0C': 'C', '0CS': 'A', '0DC': 'C', '0DG': 'G', '0DT': 'T', '0FL': 'A',
    '0G': 'G', '0NC': 'A', '0SP': 'A', '0U': 'U', '10C': 'C', '125': 'U',
    '126': 'U', '127': 'U', '128': 'N', '12A': 'A', '143': 'C', '193': 'X',
    '1AP': 'A', '1MA': 'A', '1MG': 'G', '1PA': 'F', '1PI': 'A', '1PR': 'N',
    '1SC': 'C', '1TQ': 'W', '1TY': 'Y', '1X6': 'S', '200': 'F', '23F': 'F',
    '23S': 'X', '26B': 'T', '2AD': 'X', '2AG': 'A', '2AO': 'X', '2AR': 'A',
    '2AS': 'X', '2AT': 'T', '2AU': 'U', '2BD': 'I', '2BT': 'T', '2BU': 'A',
    '2CO': 'C', '2DA': 'A', '2DF': 'N', '2DM': 'N', '2DO': 'X', '2DT': 'T',
    '2EG': 'G', '2FE': 'N', '2FI': 'N', '2FM': 'M', '2GT': 'T', '2HF': 'H',
    '2LU': 'L', '2MA': 'A', '2MG': 'G', '2ML': 'L', '2MR': 'R', '2MT': 'P',
    '2MU': 'U', '2NT': 'T', '2OM': 'U', '2OT': 'T', '2PI': 'X', '2PR': 'G',
    '2SA': 'N', '2SI': 'X', '2ST': 'T', '2TL': 'T', '2TY': 'Y', '2VA': 'V',
    '2XA': 'C', '32S': 'X', '32T': 'X', '3AH': 'H', '3AR': 'X', '3CF': 'F',
    '3DA': 'A', '3DR': 'N', '3GA': 'A', '3MD': 'D', '3ME': 'U', '3NF': 'Y',
    '3QN': 'K', '3TY': 'X', '3XH': 'G', '4AC': 'N', '4BF': 'Y', '4CF': 'F',
    '4CY': 'M', '4DP': 'W', '4FB': 'P', '4FW': 'W', '4HT': 'W', '4IN': 'W',
    '4MF': 'N', '4MM': 'X', '4OC': 'C', '4PC': 'C', '4PD': 'C', '4PE': 'C',
    '4PH': 'F', '4SC': 'C', '4SU': 'U', '4TA': 'N', '4U7': 'A', '56A': 'H',
    '5AA': 'A', '5AB': 'A', '5AT': 'T', '5BU': 'U', '5CG': 'G', '5CM': 'C',
    '5CS': 'C', '5FA': 'A', '5FC': 'C', '5FU': 'U', '5HP': 'E', '5HT': 'T',
    '5HU': 'U', '5IC': 'C', '5IT': 'T', '5IU': 'U', '5MC': 'C', '5MD': 'N',
    '5MU': 'U', '5NC': 'C', '5PC': 'C', '5PY': 'T', '5SE': 'U', '64T': 'T',
    '6CL': 'K', '6CT': 'T', '6CW': 'W', '6HA': 'A', '6HC': 'C', '6HG': 'G',
    '6HN': 'K', '6HT': 'T', '6IA': 'A', '6MA': 'A', '6MC': 'A', '6MI': 'N',
    '6MT': 'A', '6MZ': 'N', '6OG': 'G', '70U': 'U', '7DA': 'A', '7GU': 'G',
    '7JA': 'I', '7MG': 'G', '8AN': 'A', '8FG': 'G', '8MG': 'G', '8OG': 'G',
    '9NE': 'E', '9NF': 'F', '9NR': 'R', '9NV': 'V', 'A': 'A', 'A1P': 'N',
    'A23': 'A', 'A2L': 'A', 'A2M': 'A', 'A34': 'A', 'A35': 'A', 'A38': 'A',
    'A39': 'A', 'A3A': 'A', 'A3P': 'A', 'A40': 'A', 'A43': 'A', 'A44': 'A',
    'A47': 'A', 'A5L': 'A', 'A5M': 'C', 'A5N': 'N', 'A5O': 'A', 'A66': 'X',
    'AA3': 'A', 'AA4': 'A', 'AAR': 'R', 'AB7': 'X', 'ABA': 'A', 'ABR': 'A',
    'ABS': 'A', 'ABT': 'N', 'ACB': 'D', 'ACL': 'R', 'AD2': 'A', 'ADD': 'X',
    'ADX': 'N', 'AEA': 'X', 'AEI': 'D', 'AET': 'A', 'AFA': 'N', 'AFF': 'N',
    'AFG': 'G', 'AGM': 'R', 'AGT': 'C', 'AHB': 'N', 'AHH': 'X', 'AHO': 'A',
    'AHP': 'A', 'AHS': 'X', 'AHT': 'X', 'AIB': 'A', 'AKL': 'D', 'AKZ': 'D',
    'ALA': 'A', 'ALC': 'A', 'ALM': 'A', 'ALN': 'A', 'ALO': 'T', 'ALQ': 'X',
    'ALS': 'A', 'ALT': 'A', 'ALV': 'A', 'ALY': 'K', 'AN8': 'A', 'AP7': 'A',
    'APE': 'X', 'APH': 'A', 'API': 'K', 'APK': 'K', 'APM': 'X', 'APP': 'X',
    'AR2': 'R', 'AR4': 'E', 'AR7': 'R', 'ARG': 'R', 'ARM': 'R', 'ARO': 'R',
    'ARV': 'X', 'AS': 'A', 'AS2': 'D', 'AS9': 'X', 'ASA': 'D', 'ASB': 'D',
    'ASI': 'D', 'ASK': 'D', 'ASL': 'D', 'ASM': 'X', 'ASN': 'N', 'ASP': 'D',
    'ASQ': 'D', 'ASU': 'N', 'ASX': 'B', 'ATD': 'T', 'ATL': 'T', 'ATM': 'T',
    'AVC': 'A', 'AVN': 'X', 'AYA': 'A', 'AZK': 'K', 'AZS': 'S', 'AZY': 'Y',
    'B1F': 'F', 'B1P': 'N', 'B2A': 'A', 'B2F': 'F', 'B2I': 'I', 'B2V': 'V',
    'B3A': 'A', 'B3D': 'D', 'B3E': 'E', 'B3K': 'K', 'B3L': 'X', 'B3M': 'X',
    'B3Q': 'X', 'B3S': 'S', 'B3T': 'X', 'B3U': 'H', 'B3X': 'N', 'B3Y': 'Y',
    'BB6': 'C', 'BB7': 'C', 'BB8': 'F', 'BB9': 'C', 'BBC': 'C', 'BCS': 'C',
    'BE2': 'X', 'BFD': 'D', 'BG1': 'S', 'BGM': 'G', 'BH2': 'D', 'BHD': 'D',
    'BIF': 'F', 'BIL': 'X', 'BIU': 'I', 'BJH': 'X', 'BLE': 'L', 'BLY': 'K',
    'BMP': 'N', 'BMT': 'T', 'BNN': 'F', 'BNO': 'X', 'BOE': 'T', 'BOR': 'R',
    'BPE': 'C', 'BRU': 'U', 'BSE': 'S', 'BT5': 'N', 'BTA': 'L', 'BTC': 'C',
    'BTR': 'W', 'BUC': 'C', 'BUG': 'V', 'BVP': 'U', 'BZG': 'N', 'C': 'C',
    'C1X': 'K', 'C25': 'C', 'C2L': 'C', 'C2S': 'C', 'C31': 'C', 'C32': 'C',
    'C34': 'C', 'C36': 'C', 'C37': 'C', 'C38': 'C', 'C3Y': 'C', 'C42': 'C',
    'C43': 'C', 'C45': 'C', 'C46': 'C', 'C49': 'C', 'C4R': 'C', 'C4S': 'C',
    'C5C': 'C', 'C66': 'X', 'C6C': 'C', 'CAF': 'C', 'CAL': 'X', 'CAR': 'C',
    'CAS': 'C', 'CAV': 'X', 'CAY': 'C', 'CB2': 'C', 'CBR': 'C', 'CBV': 'C',
    'CCC': 'C', 'CCL': 'K', 'CCS': 'C', 'CDE': 'X', 'CDV': 'X', 'CDW': 'C',
    'CEA': 'C', 'CFL': 'C', 'CG1': 'G', 'CGA': 'E', 'CGU': 'E', 'CH': 'C',
    'CHF': 'X', 'CHG': 'X', 'CHP': 'G', 'CHS': 'X', 'CIR': 'R', 'CLE': 'L',
    'CLG': 'K', 'CLH': 'K', 'CM0': 'N', 'CME': 'C', 'CMH': 'C', 'CML': 'C',
    'CMR': 'C', 'CMT': 'C', 'CNU': 'U', 'CP1': 'C', 'CPC': 'X', 'CPI': 'X',
    'CR5': 'G', 'CS0': 'C', 'CS1': 'C', 'CS3': 'C', 'CS4': 'C', 'CS8': 'N',
    'CSA': 'C', 'CSB': 'C', 'CSD': 'C', 'CSE': 'C', 'CSF': 'C', 'CSI': 'G',
    'CSJ': 'C', 'CSL': 'C', 'CSO': 'C', 'CSP': 'C', 'CSR': 'C', 'CSS': 'C',
    'CSU': 'C', 'CSW': 'C', 'CSX': 'C', 'CSZ': 'C', 'CTE': 'W', 'CTG': 'T',
    'CTH': 'T', 'CUC': 'X', 'CWR': 'S', 'CXM': 'M', 'CY0': 'C', 'CY1': 'C',
    'CY3': 'C', 'CY4': 'C', 'CYA': 'C', 'CYD': 'C', 'CYF': 'C', 'CYG': 'C',
    'CYJ': 'X', 'CYM': 'C', 'CYQ': 'C', 'CYR': 'C', 'CYS': 'C', 'CZ2': 'C',
    'CZZ': 'C', 'D11': 'T', 'D1P': 'N', 'D3': 'N', 'D33': 'N', 'D3P': 'G',
    'D3T': 'T', 'D4M': 'T', 'D4P': 'X', 'DA': 'A', 'DA2': 'X', 'DAB': 'A',
    'DAH': 'F', 'DAL': 'A', 'DAR': 'R', 'DAS': 'D', 'DBB': 'T', 'DBM': 'N',
    'DBS': 'S', 'DBU': 'T', 'DBY': 'Y', 'DBZ': 'A', 'DC': 'C', 'DC2': 'C',
    'DCG': 'G', 'DCI': 'X', 'DCL': 'X', 'DCT': 'C', 'DCY': 'C', 'DDE': 'H',
    'DDG': 'G', 'DDN': 'U', 'DDX': 'N', 'DFC': 'C', 'DFG': 'G', 'DFI': 'X',
    'DFO': 'X', 'DFT': 'N', 'DG': 'G', 'DGH': 'G', 'DGI': 'G', 'DGL': 'E',
    'DGN': 'Q', 'DHA': 'S', 'DHI': 'H', 'DHL': 'X', 'DHN': 'V', 'DHP': 'X',
    'DHU': 'U', 'DHV': 'V', 'DI': 'I', 'DIL': 'I', 'DIR': 'R', 'DIV': 'V',
    'DLE': 'L', 'DLS': 'K', 'DLY': 'K', 'DM0': 'K', 'DMH': 'N', 'DMK': 'D',
    'DMT': 'X', 'DN': 'N', 'DNE': 'L', 'DNG': 'L', 'DNL': 'K', 'DNM': 'L',
    'DNP': 'A', 'DNR': 'C', 'DNS': 'K', 'DOA': 'X', 'DOC': 'C', 'DOH': 'D',
    'DON': 'L', 'DPB': 'T', 'DPH': 'F', 'DPL': 'P', 'DPP': 'A', 'DPQ': 'Y',
    'DPR': 'P', 'DPY': 'N', 'DRM': 'U', 'DRP': 'N', 'DRT': 'T', 'DRZ': 'N',
    'DSE': 'S', 'DSG': 'N', 'DSN': 'S', 'DSP': 'D', 'DT': 'T', 'DTH': 'T',
    'DTR': 'W', 'DTY': 'Y', 'DU': 'U', 'DVA': 'V', 'DXD': 'N', 'DXN': 'N',
    'DYS': 'C', 'DZM': 'A', 'E': 'A', 'E1X': 'A', 'ECC': 'Q', 'EDA': 'A',
    'EFC': 'C', 'EHP': 'F', 'EIT': 'T', 'ENP': 'N', 'ESB': 'Y', 'ESC': 'M',
    'EXB': 'X', 'EXY': 'L', 'EY5': 'N', 'EYS': 'X', 'F2F': 'F', 'FA2': 'A',
    'FA5': 'N', 'FAG': 'N', 'FAI': 'N', 'FB5': 'A', 'FB6': 'A', 'FCL': 'F',
    'FFD': 'N', 'FGA': 'E', 'FGL': 'G', 'FGP': 'S', 'FHL': 'X', 'FHO': 'K',
    'FHU': 'U', 'FLA': 'A', 'FLE': 'L', 'FLT': 'Y', 'FME': 'M', 'FMG': 'G',
    'FMU': 'N', 'FOE': 'C', 'FOX': 'G', 'FP9': 'P', 'FPA': 'F', 'FRD': 'X',
    'FT6': 'W', 'FTR': 'W', 'FTY': 'Y', 'FVA': 'V', 'FZN': 'K', 'G': 'G',
    'G25': 'G', 'G2L': 'G', 'G2S': 'G', 'G31': 'G', 'G32': 'G', 'G33': 'G',
    'G36': 'G', 'G38': 'G', 'G42': 'G', 'G46': 'G', 'G47': 'G', 'G48': 'G',
    'G49': 'G', 'G4P': 'N', 'G7M': 'G', 'GAO': 'G', 'GAU': 'E', 'GCK': 'C',
    'GCM': 'X', 'GDP': 'G', 'GDR': 'G', 'GFL': 'G', 'GGL': 'E', 'GH3': 'G',
    'GHG': 'Q', 'GHP': 'G', 'GL3': 'G', 'GLH': 'Q', 'GLJ': 'E', 'GLK': 'E',
    'GLM': 'X', 'GLN': 'Q', 'GLQ': 'E', 'GLU': 'E', 'GLX': 'Z', 'GLY': 'G',
    'GLZ': 'G', 'GMA': 'E', 'GMS': 'G', 'GMU': 'U', 'GN7': 'G', 'GND': 'X',
    'GNE': 'N', 'GOM': 'G', 'GPL': 'K', 'GS': 'G', 'GSC': 'G', 'GSR': 'G',
    'GSS': 'G', 'GSU': 'E', 'GT9': 'C', 'GTP': 'G', 'GVL': 'X', 'H2U': 'U',
    'H5M': 'P', 'HAC': 'A', 'HAR': 'R', 'HBN': 'H', 'HCS': 'X', 'HDP': 'U',
    'HEU': 'U', 'HFA': 'X', 'HGL': 'X', 'HHI': 'H', 'HIA': 'H', 'HIC': 'H',
    'HIP': 'H', 'HIQ': 'H', 'HIS': 'H', 'HL2': 'L', 'HLU': 'L', 'HMR': 'R',
    'HOL': 'N', 'HPC': 'F', 'HPE': 'F', 'HPH': 'F', 'HPQ': 'F', 'HQA': 'A',
    'HRG': 'R', 'HRP': 'W', 'HS8': 'H', 'HS9': 'H', 'HSE': 'S', 'HSL': 'S',
    'HSO': 'H', 'HTI': 'C', 'HTN': 'N', 'HTR': 'W', 'HV5': 'A', 'HVA': 'V',
    'HY3': 'P', 'HYP': 'P', 'HZP': 'P', 'I': 'I', 'I2M': 'I', 'I58': 'K',
    'I5C': 'C', 'IAM': 'A', 'IAR': 'R', 'IAS': 'D', 'IC': 'C', 'IEL': 'K',
    'IG': 'G', 'IGL': 'G', 'IGU': 'G', 'IIL': 'I', 'ILE': 'I', 'ILG': 'E',
    'ILX': 'I', 'IMC': 'C', 'IML': 'I', 'IOY': 'F', 'IPG': 'G', 'IPN': 'N',
    'IRN': 'N', 'IT1': 'K', 'IU': 'U', 'IYR': 'Y', 'IYT': 'T', 'IZO': 'M',
    'JJJ': 'C', 'JJK': 'C', 'JJL': 'C', 'JW5': 'N', 'K1R': 'C', 'KAG': 'G',
    'KCX': 'K', 'KGC': 'K', 'KNB': 'A', 'KOR': 'M', 'KPI': 'K', 'KST': 'K',
    'KYQ': 'K', 'L2A': 'X', 'LA2': 'K', 'LAA': 'D', 'LAL': 'A', 'LBY': 'K',
    'LC': 'C', 'LCA': 'A', 'LCC': 'N', 'LCG': 'G', 'LCH': 'N', 'LCK': 'K',
    'LCX': 'K', 'LDH': 'K', 'LED': 'L', 'LEF': 'L', 'LEH': 'L', 'LEI': 'V',
    'LEM': 'L', 'LEN': 'L', 'LET': 'X', 'LEU': 'L', 'LEX': 'L', 'LG': 'G',
    'LGP': 'G', 'LHC': 'X', 'LHU': 'U', 'LKC': 'N', 'LLP': 'K', 'LLY': 'K',
    'LME': 'E', 'LMF': 'K', 'LMQ': 'Q', 'LMS': 'N', 'LP6': 'K', 'LPD': 'P',
    'LPG': 'G', 'LPL': 'X', 'LPS': 'S', 'LSO': 'X', 'LTA': 'X', 'LTR': 'W',
    'LVG': 'G', 'LVN': 'V', 'LYF': 'K', 'LYK': 'K', 'LYM': 'K', 'LYN': 'K',
    'LYR': 'K', 'LYS': 'K', 'LYX': 'K', 'LYZ': 'K', 'M0H': 'C', 'M1G': 'G',
    'M2G': 'G', 'M2L': 'K', 'M2S': 'M', 'M30': 'G', 'M3L': 'K', 'M5M': 'C',
    'MA': 'A', 'MA6': 'A', 'MA7': 'A', 'MAA': 'A', 'MAD': 'A', 'MAI': 'R',
    'MBQ': 'Y', 'MBZ': 'N', 'MC1': 'S', 'MCG': 'X', 'MCL': 'K', 'MCS': 'C',
    'MCY': 'C', 'MD3': 'C', 'MD6': 'G', 'MDH': 'X', 'MDR': 'N', 'MEA': 'F',
    'MED': 'M', 'MEG': 'E', 'MEN': 'N', 'MEP': 'U', 'MEQ': 'Q', 'MET': 'M',
    'MEU': 'G', 'MF3': 'X', 'MG1': 'G', 'MGG': 'R', 'MGN': 'Q', 'MGQ': 'A',
    'MGV': 'G', 'MGY': 'G', 'MHL': 'L', 'MHO': 'M', 'MHS': 'H', 'MIA': 'A',
    'MIS': 'S', 'MK8': 'L', 'ML3': 'K', 'MLE': 'L', 'MLL': 'L', 'MLY': 'K',
    'MLZ': 'K', 'MME': 'M', 'MMO': 'R', 'MMT': 'T', 'MND': 'N', 'MNL': 'L',
    'MNU': 'U', 'MNV': 'V', 'MOD': 'X', 'MP8': 'P', 'MPH': 'X', 'MPJ': 'X',
    'MPQ': 'G', 'MRG': 'G', 'MSA': 'G', 'MSE': 'M', 'MSL': 'M', 'MSO': 'M',
    'MSP': 'X', 'MT2': 'M', 'MTR': 'T', 'MTU': 'A', 'MTY': 'Y', 'MVA': 'V',
    'N': 'N', 'N10': 'S', 'N2C': 'X', 'N5I': 'N', 'N5M': 'C', 'N6G': 'G',
    'N7P': 'P', 'NA8': 'A', 'NAL': 'A', 'NAM': 'A', 'NB8': 'N', 'NBQ': 'Y',
    'NC1': 'S', 'NCB': 'A', 'NCX': 'N', 'NCY': 'X', 'NDF': 'F', 'NDN': 'U',
    'NEM': 'H', 'NEP': 'H', 'NF2': 'N', 'NFA': 'F', 'NHL': 'E', 'NIT': 'X',
    'NIY': 'Y', 'NLE': 'L', 'NLN': 'L', 'NLO': 'L', 'NLP': 'L', 'NLQ': 'Q',
    'NMC': 'G', 'NMM': 'R', 'NMS': 'T', 'NMT': 'T', 'NNH': 'R', 'NP3': 'N',
    'NPH': 'C', 'NPI': 'A', 'NSK': 'X', 'NTY': 'Y', 'NVA': 'V', 'NYM': 'N',
    'NYS': 'C', 'NZH': 'H', 'O12': 'X', 'O2C': 'N', 'O2G': 'G', 'OAD': 'N',
    'OAS': 'S', 'OBF': 'X', 'OBS': 'X', 'OCS': 'C', 'OCY': 'C', 'ODP': 'N',
    'OHI': 'H', 'OHS': 'D', 'OIC': 'X', 'OIP': 'I', 'OLE': 'X', 'OLT': 'T',
    'OLZ': 'S', 'OMC': 'C', 'OMG': 'G', 'OMT': 'M', 'OMU': 'U', 'ONE': 'U',
    'ONH': 'A', 'ONL': 'X', 'OPR': 'R', 'ORN': 'A', 'ORQ': 'R', 'OSE': 'S',
    'OTB': 'X', 'OTH': 'T', 'OTY': 'Y', 'OXX': 'D', 'P': 'G', 'P1L': 'C',
    'P1P': 'N', 'P2T': 'T', 'P2U': 'U', 'P2Y': 'P', 'P5P': 'A', 'PAQ': 'Y',
    'PAS': 'D', 'PAT': 'W', 'PAU': 'A', 'PBB': 'C', 'PBF': 'F', 'PBT': 'N',
    'PCA': 'E', 'PCC': 'P', 'PCE': 'X', 'PCS': 'F', 'PDL': 'X', 'PDU': 'U',
    'PEC': 'C', 'PF5': 'F', 'PFF': 'F', 'PFX': 'X', 'PG1': 'S', 'PG7': 'G',
    'PG9': 'G', 'PGL': 'X', 'PGN': 'G', 'PGP': 'G', 'PGY': 'G', 'PHA': 'F',
    'PHD': 'D', 'PHE': 'F', 'PHI': 'F', 'PHL': 'F', 'PHM': 'F', 'PIV': 'X',
    'PLE': 'L', 'PM3': 'F', 'PMT': 'C', 'POM': 'P', 'PPN': 'F', 'PPU': 'A',
    'PPW': 'G', 'PQ1': 'N', 'PR3': 'C', 'PR5': 'A', 'PR9': 'P', 'PRN': 'A',
    'PRO': 'P', 'PRS': 'P', 'PSA': 'F', 'PSH': 'H', 'PST': 'T', 'PSU': 'U',
    'PSW': 'C', 'PTA': 'X', 'PTH': 'Y', 'PTM': 'Y', 'PTR': 'Y', 'PU': 'A',
    'PUY': 'N', 'PVH': 'H', 'PVL': 'X', 'PYA': 'A', 'PYO': 'U', 'PYX': 'C',
    'PYY': 'N', 'QMM': 'Q', 'QPA': 'C', 'QPH': 'F', 'QUO': 'G', 'R': 'A',
    'R1A': 'C', 'R4K': 'W', 'RE0': 'W', 'RE3': 'W', 'RIA': 'A', 'RMP': 'A',
    'RON': 'X', 'RT': 'T', 'RTP': 'N', 'S1H': 'S', 'S2C': 'C', 'S2D': 'A',
    'S2M': 'T', 'S2P': 'A', 'S4A': 'A', 'S4C': 'C', 'S4G': 'G', 'S4U': 'U',
    'S6G': 'G', 'SAC': 'S', 'SAH': 'C', 'SAR': 'G', 'SBL': 'S', 'SC': 'C',
    'SCH': 'C', 'SCS': 'C', 'SCY': 'C', 'SD2': 'X', 'SDG': 'G', 'SDP': 'S',
    'SEB': 'S', 'SEC': 'A', 'SEG': 'A', 'SEL': 'S', 'SEM': 'S', 'SEN': 'S',
    'SEP': 'S', 'SER': 'S', 'SET': 'S', 'SGB': 'S', 'SHC': 'C', 'SHP': 'G',
    'SHR': 'K', 'SIB': 'C', 'SLA': 'P', 'SLR': 'P', 'SLZ': 'K', 'SMC': 'C',
    'SME': 'M', 'SMF': 'F', 'SMP': 'A', 'SMT': 'T', 'SNC': 'C', 'SNN': 'N',
    'SOC': 'C', 'SOS': 'N', 'SOY': 'S', 'SPT': 'T', 'SRA': 'A', 'SSU': 'U',
    'STY': 'Y', 'SUB': 'X', 'SUN': 'S', 'SUR': 'U', 'SVA': 'S', 'SVV': 'S',
    'SVW': 'S', 'SVX': 'S', 'SVY': 'S', 'SVZ': 'X', 'SYS': 'C', 'T': 'T',
    'T11': 'F', 'T23': 'T', 'T2S': 'T', 'T2T': 'N', 'T31': 'U', 'T32': 'T',
    'T36': 'T', 'T37': 'T', 'T38': 'T', 'T39': 'T', 'T3P': 'T', 'T41': 'T',
    'T48': 'T', 'T49': 'T', 'T4S': 'T', 'T5O': 'U', 'T5S': 'T', 'T66': 'X',
    'T6A': 'A', 'TA3': 'T', 'TA4': 'X', 'TAF': 'T', 'TAL': 'N', 'TAV': 'D',
    'TBG': 'V', 'TBM': 'T', 'TC1': 'C', 'TCP': 'T', 'TCQ': 'Y', 'TCR': 'W',
    'TCY': 'A', 'TDD': 'L', 'TDY': 'T', 'TFE': 'T', 'TFO': 'A', 'TFQ': 'F',
    'TFT': 'T', 'TGP': 'G', 'TH6': 'T', 'THC': 'T', 'THO': 'X', 'THR': 'T',
    'THX': 'N', 'THZ': 'R', 'TIH': 'A', 'TLB': 'N', 'TLC': 'T', 'TLN': 'U',
    'TMB': 'T', 'TMD': 'T', 'TNB': 'C', 'TNR': 'S', 'TOX': 'W', 'TP1': 'T',
    'TPC': 'C', 'TPG': 'G', 'TPH': 'X', 'TPL': 'W', 'TPO': 'T', 'TPQ': 'Y',
    'TQI': 'W', 'TQQ': 'W', 'TRF': 'W', 'TRG': 'K', 'TRN': 'W', 'TRO': 'W',
    'TRP': 'W', 'TRQ': 'W', 'TRW': 'W', 'TRX': 'W', 'TS': 'N', 'TST': 'X',
    'TT': 'N', 'TTD': 'T', 'TTI': 'U', 'TTM': 'T', 'TTQ': 'W', 'TTS': 'Y',
    'TY1': 'Y', 'TY2': 'Y', 'TY3': 'Y', 'TY5': 'Y', 'TYB': 'Y', 'TYI': 'Y',
    'TYJ': 'Y', 'TYN': 'Y', 'TYO': 'Y', 'TYQ': 'Y', 'TYR': 'Y', 'TYS': 'Y',
    'TYT': 'Y', 'TYU': 'N', 'TYW': 'Y', 'TYX': 'X', 'TYY': 'Y', 'TZB': 'X',
    'TZO': 'X', 'U': 'U', 'U25': 'U', 'U2L': 'U', 'U2N': 'U', 'U2P': 'U',
    'U31': 'U', 'U33': 'U', 'U34': 'U', 'U36': 'U', 'U37': 'U', 'U8U': 'U',
    'UAR': 'U', 'UCL': 'U', 'UD5': 'U', 'UDP': 'N', 'UFP': 'N', 'UFR': 'U',
    'UFT': 'U', 'UMA': 'A', 'UMP': 'U', 'UMS': 'U', 'UN1': 'X', 'UN2': 'X',
    'UNK': 'X', 'UR3': 'U', 'URD': 'U', 'US1': 'U', 'US2': 'U', 'US3': 'T',
    'US5': 'U', 'USM': 'U', 'VAD': 'V', 'VAF': 'V', 'VAL': 'V', 'VB1': 'K',
    'VDL': 'X', 'VLL': 'X', 'VLM': 'X', 'VMS': 'X', 'VOL': 'X', 'X': 'G',
    'X2W': 'E', 'X4A': 'N', 'XAD': 'A', 'XAE': 'N', 'XAL': 'A', 'XAR': 'N',
    'XCL': 'C', 'XCN': 'C', 'XCP': 'X', 'XCR': 'C', 'XCS': 'N', 'XCT': 'C',
    'XCY': 'C', 'XGA': 'N', 'XGL': 'G', 'XGR': 'G', 'XGU': 'G', 'XPR': 'P',
    'XSN': 'N', 'XTH': 'T', 'XTL': 'T', 'XTR': 'T', 'XTS': 'G', 'XTY': 'N',
    'XUA': 'A', 'XUG': 'G', 'XX1': 'K', 'Y': 'A', 'YCM': 'C', 'YG': 'G',
    'YOF': 'Y', 'YRR': 'N', 'YYG': 'G', 'Z': 'C', 'Z01': 'A', 'ZAD': 'A',
    'ZAL': 'A', 'ZBC': 'C', 'ZBU': 'U', 'ZCL': 'F', 'ZCY': 'C', 'ZDU': 'U',
    'ZFB': 'X', 'ZGU': 'G', 'ZHP': 'N', 'ZTH': 'T', 'ZU0': 'T', 'ZZJ': 'A',
}
