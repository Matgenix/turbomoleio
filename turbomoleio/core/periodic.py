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
import re
import numpy as np

from collections import namedtuple
from fnmatch import fnmatch
from monty.json import MSONable, MontyDecoder
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Molecule
from pymatgen.core.structure import Structure
from pymatgen.core.units import ang_to_bohr, bohr_to_ang
from turbomoleio.core.datagroups import DataGroups
from turbomoleio.core.base import BaseSystem, get_mol_and_indices_frozen


class PeriodicSystem(BaseSystem):
    """
    Object providing information about the geometry of a periodic structure and its dynamics.
    The geometry is described using a pymatgen Structure and frozen coordinates are
    provided as additional attributes.

    Notably, since this can be used as an input for define to update the atomic
    coordinates, the constraints may be inconsistent with the actual cartesian
    coordinates. When parsed from a coord file obtained from define the information
    should be consistent. Methods are implemented to check such consistency.
    """

    def __init__(self, structure, frozen_indices=None, periodicity=3):
        """

        Args:
            structure (Structure): a pymatgen Structure object with the geometry of
                the system. Only supports ordered structures.
            frozen_indices (set): a set of indices (0-based) indicating atoms
                that should have fixed cartesian coordinates.
            periodicity (int): periodicity of the system when a Structure object is given.
                Can be 1, 2, or 3.
        """
        if not structure.is_ordered:
            raise ValueError("PeriodicSystem does not handle disordered structures.")

        self.frozen_indices = set(frozen_indices) if frozen_indices else set()
        if periodicity not in (1, 2, 3):
            raise ValueError("Periodicity should be 1, 2 or 3. For molecules (i.e. periodicity = 0), "
                             "use the MoleculeSystem)")
        if not isinstance(structure, Structure):
            raise ValueError("A Structure object should be provided for periodic systems.")
        super().__init__(molecule_or_structure=structure, frozen_indices=frozen_indices, periodicity=periodicity)

    @property
    def structure(self):
        return self._molecule_or_structure

    @classmethod
    def from_string(cls, string, fmt="coord", periodic_extension=5.0):
        """
        Creates an instance from a string. Could be the string of a coord file
        or any format supported by pymatgen Structure.

        Args:
            string (str): the string with the data.
            fmt (str): the format of the data. could be "coord" for Turbomole
                coord file or any format supported in pymatgen Structure.

        Returns:
            An instance of PeriodicSystem.
        """

        if fmt == "coord":
            dg = DataGroups(string=string)
            coordinates_str = dg.sdg("$coord", strict=True)
            if not coordinates_str:
                raise ValueError("The string does not contain $coord!")
            cell_str = dg.sdg("$cell", strict=True)
            lattice_str = dg.sdg("$lattice", strict=True)
            periodic_str = dg.sdg("$periodic", strict=True)
            if periodic_str is None:
                raise ValueError('The $periodic data group should be set for periodic systems.')
            if cell_str is None and lattice_str is None:
                raise ValueError('The $cell or $lattice data group should be set for periodic systems.')
            if cell_str is not None and lattice_str is not None:
                raise ValueError('Only one of $cell and $lattice data group should be set for periodic systems.')
            struct, fi = get_mol_and_indices_frozen(coordinates_str, cell_string=cell_str,
                                                    lattice_string=lattice_str,
                                                    periodic_string=periodic_str,
                                                    periodic_extension=periodic_extension)
            periodicity = int(periodic_str.strip())

            int_def_str = dg.sdg("$intdef", strict=True)
            if int_def_str is not None:
                raise ValueError('Internal definitions for periodic systems is not supported.')

            user_def_bonds_str = dg.sdg("$user-defined bonds", strict=True)
            if user_def_bonds_str is not None:
                raise ValueError('User-defined bonds for periodic systems is not supported.')

            return cls(struct, frozen_indices=fi, periodicity=periodicity)

        else:
            return cls(Structure.from_str(string, fmt))

    @classmethod
    def from_file(cls, filepath, fmt=None, periodic_extension=5.0):
        """
        Creates an instance from a file. Could be a coord file or any other
        format supported by pymatgen Structure.

        Args:
            filepath (str): path to the file.
            fmt (str): the format of the data. could be "coord" for Turbomole
                coord file or any format supported in pymatgen Structure.
                If None the code would try to infer it from the file name.

        Returns:
            An instance of PeriodicSystem.
        """
        fname = os.path.basename(filepath)
        if fmt == "coord" or (fmt is None and fnmatch(fname, "coord")):
            with open(filepath) as f:
                return cls.from_string(f.read(), fmt="coord", periodic_extension=periodic_extension)
        else:
            return cls(structure=Structure.from_file(filepath))

    def to_coord_string(self):
        """
        Creates the string of a coord file for this system.

        Returns:
            str: the string representing the system.
        """
        lines = ["$coord"]
        lines.extend(self.coord_lines())

        lines.append("$periodic {}".format(self.periodicity))
        lattice = self.structure.lattice
        lmat = lattice.matrix
        if (
                not np.allclose(lmat[0][1:], 0.0) or
                not np.isclose(lmat[1][2], 0.0)
        ):
            raise ValueError('Lattice should be oriented such that first vector '
                             'is aligned with x cartesian direction and second vector '
                             'is in the xy cartesian plane.')

        if self.periodicity == 1:
            lines.append(f"$cell\n  {ang_to_bohr * lattice.a}")
        elif self.periodicity == 2:
            lines.append(f"$cell\n  "
                         f"{ang_to_bohr * lattice.a}   "
                         f"{ang_to_bohr * lattice.b}   "
                         f"{lattice.gamma}")
        elif self.periodicity == 3:
            lines.append(f"$cell\n  "
                         f"{ang_to_bohr * lattice.a}   "
                         f"{ang_to_bohr * lattice.b}   "
                         f"{ang_to_bohr * lattice.c}   "
                         f"{lattice.alpha}   "
                         f"{lattice.beta}   "
                         f"{lattice.gamma}"
                         )
        else:
            raise ValueError('Periodicity should be 1, 2 or 3 for use in riper calculations.')

        lines.append("$end")

        return "\n".join(lines)

    def as_dict(self):
        """
        A JSON serializable dict representation of an object.
        Overwrites base method to handle sets.

        Returns:
            dict: json representation of the object.
        """
        # sort sets to keep consistent representations of dictionary
        d = {"@module": self.__class__.__module__, "@class": self.__class__.__name__,
             "structure": self.structure.as_dict(),
             "frozen_indices": sorted(self.frozen_indices), "periodicity": self.periodicity}

        return d

    @classmethod
    def from_dict(cls, d):
        """
        Generates object from JSON representation.
        Overwrites base method to handle sets.

        Args:
            d (dict): json representation of the object.

        Returns:
            An instance of PeriodicSystem.
        """
        d = d.copy()
        d.pop('@module', None)
        d.pop('@class', None)
        dec = MontyDecoder()
        d['structure'] = dec.process_decoded(d['structure'])
        d['frozen_indices'] = set(dec.process_decoded(i) for i in d['frozen_indices'])
        d['periodicity'] = d['periodicity']
        return cls(**d)
