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


class MoleculeSystem(BaseSystem):
    """
    Object providing information about the geometry of a molcule and its dynamics.
    The geometry is described using a pymatgen Molecule, while internal definitions,
    constraints and frozen coordinates are provided as additional attributes.

    Notably, since this can be used as an input for define to update the atomic
    coordinates, the constraints may be inconsistent with the actual cartesian
    coordinates. When parsed from a coord file obtained from define the information
    should be consistent. Methods are implemented to check such consistency.
    """

    def __init__(self, molecule, int_def=None, frozen_indices=None, user_defined_bonds=None):
        """

        Args:
            molecule (Molecule): a pymatgen Molecule object with the geometry of
                the system. Only supports ordered molecules.
            int_def (list): a list of InternalDefinition.
            frozen_indices (set): a set of indices (0-based) indicating atoms
                that should have fixed cartesian coordinates.
            user_defined_bonds (set): set of tuples with (index1, symbol, index2),
                where the atoms index are 0-based and the symbol can be "-" or "|".
        """
        self.int_def = int_def if int_def else []
        self.user_defined_bonds = set(user_defined_bonds) if user_defined_bonds else set()
        if not isinstance(molecule, Molecule):
            raise ValueError("A Molecule object should be provided for molecule systems.")
        super().__init__(molecule_or_structure=molecule, frozen_indices=frozen_indices,
                         periodicity=0)

    @property
    def molecule(self):
        return self._molecule_or_structure

    @classmethod
    def from_string(cls, string, fmt="coord"):
        """
        Creates an instance from a string. Could be the string of a coord file
        or any format supported by pymatgen Molecule.

        Args:
            string (str): the string with the data.
            fmt (str): the format of the data. could be "coord" for Turbomole
                coord file or any format supported in pymatgen Molecule.

        Returns:
            An instance of MoleculeSystem.
        """

        if fmt == "coord":
            dg = DataGroups(string=string)
            coordinates_str = dg.sdg("$coord", strict=True)
            if not coordinates_str:
                raise ValueError("The string does not contain $coord!")
            mol, fi = get_mol_and_indices_frozen(coordinates_str)

            int_def_str = dg.sdg("$intdef", strict=True)
            int_def = []
            if int_def_str:
                lines = []
                # remove empty lines and comments
                for l in int_def_str.splitlines():
                    lstrip = l.strip()
                    if lstrip and not lstrip.startswith("#"):
                        lines.append(l)
                int_def_str = "\n".join(lines)
                # split based on the presence of the index plus the status.
                # In a case like this:
                #    1 k  1.0000000000000 stre    4    1           val=   1.80084
                #    2 k  1.0000000000000 bend    4    3    1      val= 106.27756
                #         1.0000000000000 bend    3    2    1
                #         1.0000000000000 bend    2    4    1
                #    3 f  1.0000000000000 tors    1    2    3    4
                # will split in 3 groups based on the presence of the digit plus k, f, d or i
                # at the beginning of the line.
                r = r"^\s*\d+\s+[kfdi]\s+.*?(?=\s*\d+\s+[kfdi]\s+|\Z)"
                for group in re.findall(r, int_def_str, re.DOTALL|re.MULTILINE):
                    int_def.append(InternalDefinition.from_string(group))

            user_def_bonds_str = dg.sdg("$user-defined bonds", strict=True)
            user_def_bonds = set()
            if user_def_bonds_str:
                # parses a line of this form:
                # 1-2, 3-4, 5|6
                # splitting first on "," and then on "-" and "|"
                for l in user_def_bonds_str.splitlines():
                    l = l.strip()
                    if not l or l.startswith("#"):
                        continue
                    for bond in l.split(","):
                        for separator in ("-", "|"):
                            if separator in bond:
                                bond_indices = bond.split(separator)
                                if len(bond_indices) != 2:
                                    raise ValueError("Cannot parse user-defined bonds for line: {}".format(l))
                                index_1 = int(bond_indices[0])-1
                                index_2 = int(bond_indices[1])-1
                                user_def_bonds.add((index_1, separator, index_2))
                                break
                        else:
                            raise ValueError("Cannot parse user-defined bonds for line: {}".format(l))

            return cls(mol, int_def=int_def, frozen_indices=fi, user_defined_bonds=user_def_bonds)

        else:
            return cls(Molecule.from_str(string, fmt))

    @classmethod
    def from_file(cls, filepath, fmt=None):
        """
        Creates an instance from a file. Could be a coord file or any other
        format supported by pymatgen Molecule.

        Args:
            filepath (str): path to the file.
            fmt (str): the format of the data. could be "coord" for Turbomole
                coord file or any format supported in pymatgen Molecule.
                If None the code would try to infer it from the file name.

        Returns:
            An instance of MoleculeSystem.
        """
        fname = os.path.basename(filepath)
        if fmt == "coord" or (fmt is None and fnmatch(fname, "coord")):
            with open(filepath) as f:
                return cls.from_string(f.read(), fmt="coord")
        else:
            return cls(molecule=Molecule.from_file(filepath))

    def to_coord_string(self):
        """
        Creates the string of a coord file for this system.

        Returns:
            str: the string representing the system.
        """
        lines = ["$coord"]
        lines.extend(self.coord_lines())
        if self.int_def:
            lines.append("$intdef")
            for i, idef in enumerate(self.int_def, start=1):
                lines.extend(idef.to_string(index=i).splitlines())
        if self.user_defined_bonds:
            lines.append("$user-defined bonds")
            strings_list = []
            # sort to obtain reproducible outputs
            for udb in sorted(self.user_defined_bonds):
                strings_list.append("{}{}{}".format(udb[0]+1, udb[1], udb[2]+1))
            lines.append(", ".join(strings_list))

        lines.append("$end")

        return "\n".join(lines)

    def _check_index(self, indices_list):
        """
        Helper method to determine if one of the indices exceeds the number of
        atoms in the molecule or if it negative. Used when adding internal coordinates.

        Args:
            indices_list (list): list of integers with the indices.

        Raises:
            ValueError: if the index is larger than the number of sites.
        """
        n_sites = self.molecule.num_sites
        if any(n<0 or n >= n_sites for n in indices_list):
            raise ValueError("One of the indices representing the atoms is negative or larger then the number of sites")

    def add_distance(self, atom1, atom2, value=None, weights=None, status="f", add_user_def_bonds=True):
        """
        Adds a distance coordinate to the list of internal coordinates.
        Typically used to add frozen coordinates.

        Allows to set linear combinations of the coordinate. If a single value is given to the
        arguments it will be a simple internal coordinate. If a list of atoms and values are
        given this will produce a linear combination of those coordinates.

        Args:
            atom1 (int or list): index (0-based) or list on indices of the first atom(s).
            atom2 (int or list): index (0-based) or list on indices of the second atom(s).
            value (float): value of the distance.
            weights (float or list): weight or list of weights in case of linear
                combinations of coordinates. Should have the same length as the atoms.
                If None will be set to a list of 1.0 with the same length as atoms.
            status (str): the status of the coordinate, can be "k", "f", "d" and "i".
            add_user_def_bonds (bool): if True the pair of atoms will be added to the
                user defined bonds.
        """
        if np.ndim(atom1) == 0:
            atom1 = [atom1]
        if np.ndim(atom2) == 0:
            atom2 = [atom2]

        if len(atom1) != len(atom2):
            raise ValueError("list of atoms should have the same length")

        self._check_index(atom1)
        self._check_index(atom2)

        indices = np.transpose([atom1, atom2]).tolist()

        self.int_def.append(Distance(status=status, indices=indices, weights=weights, value=value))

        if add_user_def_bonds:
            for (i, j) in indices:
                self.user_defined_bonds.add((i, "-", j))

    def add_bond_angle(self, atom1, atom2, vertex, value=None, weights=None,
                       status="f", add_user_def_bonds=True):
        """
        Adds a bond angle coordinate to the list of internal coordinates.
        Typically used to add frozen coordinates.

        Allows to set linear combinations of the coordinate. If a single value is given to the
        arguments it will be a simple internal coordinate. If a list of atoms and values are
        given this will produce a linear combination of those coordinates.

        Args:
            atom1 (int or list): index (0-based) or list on indices of the first atom(s).
            atom2 (int or list): index (0-based) or list on indices of the second atom(s).
            vertex (int or list): index (0-based) or list on indices of the vertex.
            value (float): value of the distance.
            weights (float or list): weight or list of weights in case of linear
                combinations of coordinates. Should have the same length as the atoms.
                If None will be set to a list of 1.0 with the same length as atoms.
            status (str): the status of the coordinate, can be "k", "f", "d" and "i".
            add_user_def_bonds (bool): if True the pair of atoms will be added to the used defined bonds
        """
        if np.ndim(atom1) == 0:
            atom1 = [atom1]
        if np.ndim(atom2) == 0:
            atom2 = [atom2]
        if np.ndim(vertex) == 0:
            vertex = [vertex]

        if not len(atom1) == len(atom2) == len(vertex):
            raise ValueError("list of atoms should have the same length")

        self._check_index(atom1)
        self._check_index(atom2)
        self._check_index(vertex)

        indices = np.transpose([atom1, atom2, vertex]).tolist()

        self.int_def.append(BondAngle(status=status, indices=indices, weights=weights, value=value))

        if add_user_def_bonds:
            for (i, j, v) in indices:
                self.user_defined_bonds.add((i, "-", v))
                self.user_defined_bonds.add((j, "-", v))

    def add_dihedral(self, atom1, atom2, atom3, atom4, value=None, weights=None,
                     status="f", add_user_def_bonds=True):
        """
        Adds a dihedral angle coordinate to the list of internal coordinates.
        Typically used to add frozen coordinates.

        Allows to set linear combinations of the coordinate. If a single value is given to the
        arguments it will be a simple internal coordinate. If a list of atoms and values are
        given this will produce a linear combination of those coordinates.

        Args:
            atom1 (int or list): index (0-based) or list on indices of the first atom(s).
            atom2 (int or list): index (0-based) or list on indices of the second atom(s).
            atom3 (int or list): index (0-based) or list on indices of the third atom(s).
            atom4 (int or list): index (0-based) or list on indices of the fourth atom(s).
            value (float): value of the distance.
            weights (float or list): weight or list of weights in case of linear
                combinations of coordinates. Should have the same length as the atoms.
                If None will be set to a list of 1.0 with the same length as atoms.
            frozen (bool): True if the coordinate should be frozen.
            status (str): the status of the coordinate, can be "k", "f", "d" and "i".
            add_user_def_bonds (bool): if True the pair of atoms will be added to the used defined bonds
        """
        if np.ndim(atom1) == 0:
            atom1 = [atom1]
        if np.ndim(atom2) == 0:
            atom2 = [atom2]
        if np.ndim(atom3) == 0:
            atom3 = [atom3]
        if np.ndim(atom4) == 0:
            atom4 = [atom4]

        if not len(atom1) == len(atom2) == len(atom3) == len(atom4):
            raise ValueError("list of atoms should have the same length")

        self._check_index(atom1)
        self._check_index(atom2)
        self._check_index(atom3)
        self._check_index(atom4)

        indices = np.transpose([atom1, atom2, atom3, atom4]).tolist()

        self.int_def.append(DihedralAngle(status=status, indices=indices, weights=weights, value=value))

        if add_user_def_bonds:
            for (i, j, k, l) in indices:
                self.user_defined_bonds.add((i, "-", j))
                self.user_defined_bonds.add((j, "-", k))
                self.user_defined_bonds.add((k, "-", l))

    def get_int_def_inconsistencies(self, atol=1e-4, ltol=1e-4):
        """
        Gives a list of the internal coordinates whose values are not consistent
        with the values obtained from the molecule object.

        Args:
            atol (float): tolerance allowed on distances (in bohr).
            ltol (float): tolerance allowed on angles (in degrees).

        Returns:
            list: list of inconsistent internal coordinates.
        """
        errors = []

        for idef in self.int_def:
            if idef.coord_type == "length":
                tol = ltol
            elif idef.coord_type == "angle":
                tol = atol
            else:
                raise ValueError("unknow coord_type {}".format(idef.coord_type))
            if not idef.is_valid(self.molecule, tol):
                errors.append(idef)

        return errors

    def has_inconsistencies(self, atol=1e-4, ltol=1e-4):
        """
        Checks if any inconsistency in internal coordinates is present.

        Args:
            atol (float): tolerance allowed on distances (in bohr).
            ltol (float): tolerance allowed on angles (in degrees).

        Returns:
            bool: True if inconsistencies are present.
        """
        return len(self.get_int_def_inconsistencies(atol=atol, ltol=ltol)) != 0

    def as_dict(self):
        """
        A JSON serializable dict representation of an object.
        Overwrites base method to handle sets.

        Returns:
            dict: json representation of the object.
        """
        # sort sets to keep consistent representations of dictionary
        d = {"@module": self.__class__.__module__, "@class": self.__class__.__name__,
             "molecule": self.molecule.as_dict(), "int_def": [i.as_dict() for i in self.int_def],
             "frozen_indices": sorted(self.frozen_indices), "user_defined_bonds": sorted(self.user_defined_bonds)}

        return d

    @classmethod
    def from_dict(cls, d):
        """
        Generates object from JSON representation.
        Overwrites base method to handle sets.

        Args:
            d (dict): json representation of the object.

        Returns:
            An instance of MoleculeSystem.
        """
        d = d.copy()
        d.pop('@module', None)
        d.pop('@class', None)
        dec = MontyDecoder()
        d['molecule'] = dec.process_decoded(d['molecule'])
        d['int_def'] = [dec.process_decoded(i) for i in d['int_def']]
        d['frozen_indices'] = set(dec.process_decoded(i) for i in d['frozen_indices'])
        d['user_defined_bonds'] = set(tuple(i) for i in d['user_defined_bonds'])
        return cls(**d)


class InternalDefinition(abc.ABC, MSONable):
    """
    Abstract class for the definition of an internal coordinate according to
    the conventions established in turbomole.
    """
    # coord_str is the string used in the coord file to identify the type of
    # internal coord (e.g. stre, bend, tors)
    coord_str = NotImplemented
    # n_atoms is an int giving the number of atoms that should be given to each
    # internal coordinate (e.g. 2 for stre and 3 for bend)
    n_atoms = NotImplemented
    # coord_type is a string describing the type of coordinate. could be "length" or "angle".
    # used to identify which type of tolerance (angle or distance) to use in
    # general check for all the internal definitions.
    coord_type = NotImplemented

    def __init_subclass__(cls, **kwargs):
        """
        Override the method to check the presence of implemented subclass attributes.
        """
        super().__init_subclass__(**kwargs)

        if cls.coord_str is NotImplemented:
            raise NotImplementedError('Subclasses should implement the class attribute coord_str')

        if cls.n_atoms is NotImplemented:
            raise NotImplementedError('Subclasses should implement the class attribute n_atoms')

        if cls.coord_type is NotImplemented:
            raise NotImplementedError('Subclasses should implement the class attribute coord_type')

    def __init__(self, status, indices, weights=None, value=None):
        """
        Args:
            status (str): the status of the coordinate. Could be "k", "f", "d" or "i"
            indices (list): a lists of 0-based indices for the atoms involved in the
                internal coordinate. A list of lists if a linear combination should
                be considered. Stored internally as a list of lists.
            weights (float or list): weight or list of weights in case of linear
                combinations of coordinates. Should have the same length as indices.
                If None will be set to a list of 1.0 with the same length as indices.
            value (float): value of the internal coordinate. Can be None.
        """

        if status not in ("k", "f", "d", "i"):
            raise ValueError("Unsupported internal definition status: {}".format(status))

        self.status = status
        if np.ndim(indices) == 1:
            indices = [indices]

        for ind in indices:
            if len(ind) != self.n_atoms:
                raise ValueError("Wrong number of indices: {}".format(len(ind)))

        self.indices = indices
        if weights is None:
            weights = [1.0] * len(indices)
        elif np.ndim(weights) == 0:
            weights = [weights]
        self.weights = weights
        self.value = value

    @staticmethod
    def get_subclass_from_str(coord_str):
        """
        Helper method to get the correct subclass based on the string description

        Args:
            coord_str (str): the string to determine the subclass.

        Returns:
            The appropriate subclass of InternalDefinition.
        """
        for c in InternalDefinition.__subclasses__():
            if c.coord_str == coord_str:
                return c

        raise ValueError("Could not find a subclass with coord_str matching {}".format(coord_str))

    @classmethod
    def from_string(cls, string):
        """
        Generates a concrete subclass of InternalDefinition from a tm string.

        Args:
            string (str): a turbomole string used in the $intdef datagroup.

        Returns:
            An instance of a subclass of InternalDefinition.
        """

        lines = []
        for l in string.splitlines():
            l = l.strip()
            if l and not l.startswith("#"):
                lines.append(l)

        s = lines[0].split()

        # example line that will be parsed:
        # 1 f  1.0000000000000 bend    3    2    1      val= 104.95100
        status = s[1]
        weights = [float(s[2])]
        coord_cls = cls.get_subclass_from_str(s[3])
        indices = [[int(i)-1 for i in s[4:4+coord_cls.n_atoms]]]
        if "val" in string:
            last = s[-1]
            last = last.replace("val", "")
            last = last.replace("=", "")
            value = float(last)
        else:
            value = None

        if len(lines) > 1:
            for l in lines[1:]:
                s = l.split()
                weights.append(float(s[0]))
                indices.append([int(i)-1 for i in s[2:2+coord_cls.n_atoms]])

        return coord_cls(status, indices, weights, value)

    def to_string(self, index=1):
        """
        Generates the turbomole string used in $intdef.

        Args:
            index (int): the index that should be prepended to the generated string.

        Returns:
            A string that describes the coordinate and to be used inside $intdef.
        """
        first_line = "{} {} {} {} ".format(index, self.status, self.weights[0], self.coord_str)
        first_line += " ".join(str(i+1) for i in self.indices[0])
        if self.value is not None:
            first_line += " val={}".format(self.value)
        lines = [first_line]
        if len(self.indices) > 1:
            for w, ind in zip(self.weights[1:], self.indices[1:]):
                l = "    {} {} ".format(w, self.coord_str) + " ".join(str(i+1) for i in ind)
                lines.append(l)

        return "\n".join(lines)

    def is_valid(self, molecule, tol=1e-4):
        """
        If the value is not None, checks if it is consistent with the positions
        of the atoms in the given molecule.

        Args:
            molecule (Molecule): a pymatgen Molecule.
            tol (float): the tolerance of the allowed difference between the
                internal and the calculated value.

        Returns:
            bool: True if the values are consistent.
        """
        if self.value is None:
            return True

        computed_value = self.compute_value(molecule)

        # for angles use a 360 degrees module
        if self.coord_type == "angle":
            return np.abs(np.mod(computed_value, 360) - np.mod(self.value, 360)) < tol
        else:
            return np.abs(computed_value - self.value) < tol

    @abc.abstractmethod
    def compute_value(self, molecule):
        """
        Calculates the value of the internal coordinate based on the given molecule.

        Args:
            molecule (Molecule): a pymatgen Molecule.

        Returns:
            float: the calculated value of the internal coordinate.
        """
        raise NotImplementedError()


class Distance(InternalDefinition):
    """
    Class defining a distance internal coordinate (stre in TM).
    """

    coord_str = "stre"
    n_atoms = 2
    coord_type = "length"

    def compute_value(self, molecule):
        """
        Calculates the value of the internal coordinate based on the given molecule.

        Args:
            molecule (Molecule): a pymatgen Molecule.

        Returns:
            float: the calculated value of the internal coordinate.
        """
        tot_val = 0
        for weight, indices in zip(self.weights, self.indices):
            tot_val += weight * molecule.get_distance(indices[0], indices[1]) * ang_to_bohr

        tot_val /= np.abs(self.weights).sum()

        return tot_val

    def __str__(self):
        return "Distance between atoms with indices {}".format(self.indices)


class BondAngle(InternalDefinition):
    """
    Class defining a bond angle internal coordinate (bend in TM).
    """
    coord_str = "bend"
    n_atoms = 3
    coord_type = "angle"

    def compute_value(self, molecule):
        """
        Calculates the value of the internal coordinate based on the given molecule.

        Args:
            molecule (Molecule): a pymatgen Molecule.

        Returns:
            float: the calculated value of the internal coordinate.
        """

        tot_val = 0
        for weight, indices in zip(self.weights, self.indices):
            tot_val += weight * molecule.get_angle(indices[0], indices[2], indices[1])
        tot_val /= np.abs(self.weights).sum()

        return tot_val

    def __str__(self):
        return "Angle between atoms with indices {}".format(self.indices)


class DihedralAngle(InternalDefinition):
    """
    Class defining a dihedral angle internal coordinate (tors in TM).
    """
    coord_str = "tors"
    n_atoms = 4
    coord_type = "angle"

    def compute_value(self, molecule):
        """
        Calculates the value of the internal coordinate based on the given molecule.

        Args:
            molecule (Molecule): a pymatgen Molecule.

        Returns:
            float: the calculated value of the internal coordinate.
        """

        tot_val = 0
        for weight, indices in zip(self.weights, self.indices):
            tot_val += weight * molecule.get_dihedral(indices[0], indices[1], indices[2], indices[3])
        tot_val /= np.abs(self.weights).sum()

        return tot_val

    def __str__(self):
        return "Dihedral Angle between the planes defined by atoms with indices {}" .format(self.indices)


class InverseDistance(InternalDefinition):
    """
    Class defining a inverse distance internal coordinate (invr in TM).
    """
    coord_str = "invr"
    n_atoms = 2
    coord_type = "length"

    def compute_value(self, molecule):
        """
        Calculates the value of the internal coordinate based on the given molecule.

        Args:
            molecule (Molecule): a pymatgen Molecule.

        Returns:
            float: the calculated value of the internal coordinate.
        """

        tot_val = 0
        for weight, indices in zip(self.weights, self.indices):
            tot_val += weight / (molecule.get_distance(indices[0], indices[1]) * ang_to_bohr)

        tot_val /= np.abs(self.weights).sum()

        return tot_val

    def __str__(self):
        return "Inverse distance between atoms with indices {}".format(self.indices)


class OutOfPlaneAngle(InternalDefinition):
    """
    Class defining an out-of-plane angle internal coordinate (outp in TM).
    """
    coord_str = "outp"
    n_atoms = 4
    coord_type = "angle"

    def compute_value(self, molecule):
        """
        Calculates the value of the internal coordinate based on the given molecule.

        Args:
            molecule (Molecule): a pymatgen Molecule.

        Returns:
            float: the calculated value of the internal coordinate.
        """

        tot_val = 0
        coords = molecule.cart_coords
        # evaluation of outp extracted from the outp script in TM
        for weight, indices in zip(self.weights, self.indices):
            v1 = coords[indices[3]] - coords[indices[1]]
            v2 = coords[indices[3]] - coords[indices[2]]
            p = np.cross(v1, v2)
            p /= np.linalg.norm(p)
            a = coords[indices[0]] - coords[indices[3]]
            a /= np.linalg.norm(a)
            outp = np.pi/2 - np.arccos(np.dot(p, a))
            tot_val += weight * outp

        tot_val /= np.abs(self.weights).sum()
        tot_val = np.rad2deg(tot_val)

        return tot_val

    def __str__(self):
        return "Out of plane angle between atoms with indices {}".format(self.indices)


class CollinearBendingAngle(InternalDefinition):
    """
    Class defining an collinear bending angle internal coordinate (linc in TM).
    """
    coord_str = "linc"
    n_atoms = 4
    coord_type = "angle"

    def compute_value(self, molecule):
        """
        Calculates the value of the internal coordinate based on the given molecule.

        Args:
            molecule (Molecule): a pymatgen Molecule.

        Returns:
            float: the calculated value of the internal coordinate.
        """

        tot_val = 0
        coords = molecule.cart_coords
        # evaluation of outp extracted from the outp script in TM
        for weight, indices in zip(self.weights, self.indices):
            u = coords[indices[0]] - coords[indices[2]]
            v = coords[indices[3]] - coords[indices[2]]
            x = coords[indices[1]] - coords[indices[2]]
            u /= np.linalg.norm(u)
            v /= np.linalg.norm(v)
            x /= np.linalg.norm(x)
            linc = np.pi - np.arccos(np.dot(u, v)) - np.arccos(np.dot(x, v))
            tot_val += weight * linc

        tot_val /= np.abs(self.weights).sum()
        tot_val = np.rad2deg(tot_val)

        return tot_val

    def __str__(self):
        return "Collinear bending angle for atoms with indices {}".format(self.indices)


class PerpendicularBendingAngle(InternalDefinition):
    """
    Class defining an collinear bending angle internal coordinate (linp in TM).
    """
    coord_str = "linp"
    n_atoms = 4
    coord_type = "angle"

    def compute_value(self, molecule):
        """
        Calculates the value of the internal coordinate based on the given molecule.

        Args:
            molecule (Molecule): a pymatgen Molecule.

        Returns:
            float: the calculated value of the internal coordinate.
        """

        tot_val = 0
        coords = molecule.cart_coords
        # evaluation of outp extracted from the outp script in TM
        for weight, indices in zip(self.weights, self.indices):
            u = coords[indices[0]] - coords[indices[2]]
            v = coords[indices[3]] - coords[indices[2]]
            z = coords[indices[1]] - coords[indices[2]]
            u /= np.linalg.norm(u)
            v /= np.linalg.norm(v)
            z /= np.linalg.norm(z)
            w = np.cross(v,u)
            w /= np.linalg.norm(w)
            x = np.cross(z,v)
            x /= np.linalg.norm(x)
            linp = (np.pi - np.arccos(np.dot(u, x)) - np.arccos(np.dot(z, w))) / 2
            tot_val += weight * linp

        tot_val /= np.abs(self.weights).sum()
        tot_val = np.rad2deg(tot_val)

        return tot_val

    def __str__(self):
        return "Perpendicular bending angle for atoms with indices {}".format(self.indices)
