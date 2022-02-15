# -*- coding: utf-8 -*-
# The turbomoleio package, a python interface to Turbomole
# for preparing inputs, parsing outputs and other related tools.
#
# Copyright (C) 2018-2021 BASF SE, Matgenix SRL.
#
# This file is part of turbomoleio.
#
# Turbomoleio is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Turbomoleio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with turbomoleio (see ~turbomoleio/COPYING). If not,
# see <https://www.gnu.org/licenses/>.

import os
import abc
import numpy as np

from collections import namedtuple
from fnmatch import fnmatch
from monty.json import MSONable
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Molecule
from pymatgen.core.structure import Structure
from pymatgen.core.units import ang_to_bohr, bohr_to_ang


def get_mol_and_indices_frozen(string, cell_string=None, periodic_string=None,
                               lattice_string=None,
                               periodic_extension=5.0):
    """
    Extracts a pymatgen Molecule and indices of the frozen atoms from the
    $coord datagroup.

    Args:
        string (str): The string containing the list of coordinates.
        cell_string (str): The string containing the cell parameters.
            For 1D periodic systems, cell consists of a single
            number: "a" lattice parameter.
            For 2D periodic systems, cell consists of 3 numbers:
            "a" and "b" lattice parameters and Gamma angle between
            "a" and "b".
            For 3D periodic systems, cell consists of 6 numbers:
            "a", "b" and "c" lattice parameters and Alpha, Beta and
            Gamma angles between lattice parameters.
            Only cell_string or lattice_string should be set for
            periodic systems (not both).
        periodic_string (str): The string describing periodicity.
        lattice_string (str): The string containing the lattice
            parameters as a matrix.
            Only cell_string or lattice_string should be set for
            periodic systems (not both).
        periodic_extension (float): For 1D and 2D periodic systems,
            periodic_extension is the amount (in Angstroms of
            artifically added vacuum in the Structure object in
            the non-periodic directions).

    Returns:
        namedtuple:
            molecule: a pymatgen molecule or structure (when cell and
            periodic are provided).
            frozen_indices: a set with the indices of the frozen atoms.
    """

    if cell_string is not None and lattice_string is not None:
        raise ValueError('Only cell_string or lattice_string should be set.')

    lines = [line.strip() for line in string.splitlines()]

    coords = []
    species = []
    frozen_indices = set()
    iatom = -1
    # here parses lines of this form:
    # 1.0 0.2 0.1 o
    # 1.0 0.2 0.1 si
    # 1.0 0.2 0.1 h f
    # with the final f meaning that the coordinate is frozen.
    for l in lines:
        s = l.split()
        if not s or s[0].startswith("#"):
            continue
        iatom += 1
        coords.append([float(i) * bohr_to_ang for i in s[0:3]])
        species.append(s[3].capitalize())

        if len(s) == 5 and s[-1].lower() == "f":
            frozen_indices.add(iatom)

    coords = np.array(coords)
    if cell_string is not None:
        periodicity = int(periodic_string.strip())
        cell = [float(sp) for sp in cell_string.strip().split()]
        if periodicity == 1:
            if len(cell) != 1:
                raise ValueError('The $cell data group should contain 1 number for 1D periodic systems')
            yrange = np.max(coords[:, 1])-np.min(coords[:, 1])
            zrange = np.max(coords[:, 2])-np.min(coords[:, 2])
            lattice = Lattice.orthorhombic(a=bohr_to_ang*cell[0],
                                           b=yrange+periodic_extension,
                                           c=zrange+periodic_extension)
        elif periodicity == 2:
            if len(cell) != 3:
                raise ValueError('The $cell data group should contain 3 numbers for 2D periodic systems')
            zrange = np.max(coords[:, 2])-np.min(coords[:, 2])
            lattice = Lattice.from_parameters(a=bohr_to_ang*cell[0],
                                              b=bohr_to_ang*cell[1],
                                              c=zrange+periodic_extension,
                                              alpha=90.0, beta=90.0, gamma=cell[2],
                                              vesta=True)
        elif periodicity == 3:
            if len(cell) != 6:
                raise ValueError('The $cell data group should contain 6 numbers for 3D periodic systems')
            lattice = Lattice.from_parameters(a=bohr_to_ang*cell[0],
                                              b=bohr_to_ang*cell[1],
                                              c=bohr_to_ang*cell[2],
                                              alpha=cell[3],
                                              beta=cell[4],
                                              gamma=cell[5],
                                              vesta=True)
        else:
            raise RuntimeError('Periodicity should be one of 1, 2 or 3.')
        molecule_or_structure = Structure(lattice=lattice, species=species, coords=coords, coords_are_cartesian=True,
                                          to_unit_cell=False)
    elif lattice_string is not None:
        lattice_string = lattice_string.strip()
        periodicity = int(periodic_string.strip())
        lattice_numbers = bohr_to_ang * np.array(
            [[float(ln) for ln in line.split()] for line in lattice_string.splitlines()]
        )
        if periodicity == 1:
            if lattice_numbers.shape != (1, 1):
                raise ValueError('The lattice_numbers should contain 1 number for 1D periodic systems')
            yrange = np.max(coords[:, 1])-np.min(coords[:, 1])
            zrange = np.max(coords[:, 2])-np.min(coords[:, 2])
            lattice = Lattice.orthorhombic(a=lattice_numbers[0][0],
                                           b=yrange+periodic_extension,
                                           c=zrange+periodic_extension)
        elif periodicity == 2:
            if lattice_numbers.shape != (2, 2):
                raise ValueError('The lattice_numbers should contain 2x2 numbers for 2D periodic systems')
            zrange = np.max(coords[:, 2])-np.min(coords[:, 2])
            lat_matrix = [[lattice_numbers[0, 0], lattice_numbers[0, 1], 0.0],
                          [lattice_numbers[1, 0], lattice_numbers[1, 1], 0.0],
                          [0.0, 0.0, zrange+periodic_extension]]
            lattice = Lattice(lat_matrix)
        elif periodicity == 3:
            if lattice_numbers.shape != (3, 3):
                raise ValueError('The lattice_numbers should contain 3x3 numbers for 3D periodic systems')
            lattice = Lattice(lattice_numbers)
        else:
            raise RuntimeError('Periodicity should be one of 1, 2 or 3.')
        molecule_or_structure = Structure(lattice=lattice, species=species, coords=coords, coords_are_cartesian=True,
                                          to_unit_cell=False)
    else:
        molecule_or_structure = Molecule(species, coords)
    CoordSystem = namedtuple("CoordSystem", ["molecule_or_structure", "frozen_indices"])

    return CoordSystem(molecule_or_structure, frozen_indices)


def get_coord_lines(molecule, frozen_indices=None):
    """
    Generates the lines contained in the $coord datagroup from a pymatgen Molecule.

    Args:
        molecule (Molecule): the molecule that should be converted.
        frozen_indices (set): indices of the molecule that should be frozen
            (an "f" appended at the end of the line).
    Returns:
        list: list of strings with the lines of $coord.
    """
    lines = []
    frozen_indices = frozen_indices or set()
    for i, (coords, specie) in enumerate(zip(molecule.cart_coords * ang_to_bohr, molecule.species)):
        l = "{:.14f} {:.14f} {:.14f} {}".format(*coords, specie.symbol.lower())
        if i in frozen_indices:
            l += " f"
        lines.append(l)

    return lines


class BaseSystem(MSONable):
    """Common base class for MoleculeSystem and PeriodicSystem objects."""

    def __init__(self, molecule_or_structure, frozen_indices=None, periodicity=0):
        """Common base for MoleculeSystem and PeriodicSystem objects.

        Args:
            molecule_or_structure: pymatgen Molecule or Structure.
            frozen_indices (set): a set of indices (0-based) indicating atoms
                that should have fixed cartesian coordinates.
            periodicity: Periodicity of the system. Default is 0 for MoleculeSystem and
                should be 1, 2 or 3 for PeriodicSystem.
        """
        if not molecule_or_structure.is_ordered:
            raise ValueError("Turbomoleio and turbomole do not handle disordered structures.")
        self._molecule_or_structure = molecule_or_structure
        self.frozen_indices = set(frozen_indices) if frozen_indices else set()
        self.periodicity = periodicity

    def coord_lines(self):
        """
        Gives the lines of the $coord datagroup based on the internal Structure.

        Returns:
            list: list of strings with the $coord datagroup.
        """

        return get_coord_lines(self._molecule_or_structure, self.frozen_indices)

    @abc.abstractmethod
    def to_coord_string(self):
        """
        Creates the string of a coord file for this system.

        Should be implemented separately for molecules and periodic systems.
        """

    def to_coord_file(self, filepath="coord"):
        """
        Writes the system to a coord file.

        Args:
            filepath (str): the path to the file
        """
        with open(filepath, 'wt') as f:
            f.write(self.to_coord_string())

    def to_file(self, filepath, fmt=None):
        """
        Writes the system to a file depending on the specified format. Could be
        a coord file or any other format supported by pymatgen Molecule or Structure.

        Args:
            filepath (str): path to the file.
            fmt (str): the format of the data. could be "coord" for Turbomole
                coord file or any format supported in pymatgen Molecule or Structure.
                If None the code would try to infer it from the file name.
        """
        fname = os.path.basename(filepath)
        if fmt == "coord" or (fmt is None and fnmatch(fname, "*coord*")):
            self.to_coord_file(filepath)
        else:
            self._molecule_or_structure.to(filename=filepath, fmt=fmt)
