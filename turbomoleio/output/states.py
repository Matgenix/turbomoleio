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

import re
import subprocess
import numpy as np
from collections import defaultdict
from collections.abc import MutableSequence
from copy import copy
from fractions import Fraction
from monty.json import MSONable
from monty.os import cd
from turbomoleio.core.symmetry import irrep_size
from turbomoleio.core.control import Control


def get_mos_energies(dg_string):
    """
    Given the string of a mos datagroup or the content of the mos file
    extracts the list of energies associated with each state.
    Each state is identified by its irreducible representation
    and the index relative that that irrep.

    Args:
        dg_string (str): the string with the datagroup

    Returns:
        A list of energies, each element of the list contains
        the irrep, the index relative to the irred and the
        energy.
    """
    m = re.findall(r"^\s+(\d+)\s+([\w\'\"]+)\s+eigenvalue=([+-]?[0-9]*[.]?[0-9]+D[+-]\d{2})\s+",
                   dg_string, re.MULTILINE)

    if not m:
        raise RuntimeError("No mos energies could be extracted")

    energies = []

    for e in m:
        energies.append([e[1], int(e[0]), float(e[2].replace("D", "E"))])

    return energies


class State(MSONable):
    """
    Object a single state of a calculation.
    Each state is identified by its irreducible representation,
    the index in the irrep, its spin, its occupation and the
    calculated eigenvalue.
    Uses Fraction for occupation for compatibility with the
    Shells object.
    """

    def __init__(self, eigenvalue, irrep, irrep_index, occupation, spin=None):
        """
        Args:
            eigenvalue (float): the eigenvalue of the state
            irrep (str): the irreducible representation
            irrep_index (int): the index in the irreducible representation
            occupation (Fraction): occupation of the state
            spin (str): "a" or "b" if a uhf calculation, representing the spin,
                None otherwise.
        """
        self.eigenvalue = eigenvalue
        self.irrep = irrep
        self.irrep_index = irrep_index
        self.occupation = occupation
        self.spin = spin

    def __eq__(self, other):
        """
        Equal if all the attributes are the same

        Args:
            other: another object

        Returns:
            (bool): True if the objects are equals
        """
        if not isinstance(other, self.__class__):
            return False

        return self.eigenvalue == other.eigenvalue and self.irrep == other.irrep and \
               self.irrep_index == other.irrep_index and self.occupation == other.occupation \
               and self.spin == other.spin

    def __str__(self):
        return "Eigenvalue: {}, irrep: {}, index: {}, spin: {}, occupation: {}".format(
            self.eigenvalue, self.irrep, self.irrep_index, self.spin, self.occupation)

    @property
    def max_occupation(self):
        """
        The maximum possible occupation for this state based on the irrep and the spin.

        Returns:
            int: the maximum occupation.
        """
        max_occ = irrep_size[self.irrep[0].upper()]
        if not self.spin:
            max_occ *= 2

        return max_occ

    @property
    def has_fractional_occ(self):
        """
        True if the state is only fractionally occupied.
        """
        return 0 < self.occupation < self.max_occupation

    def as_dict(self):
        """
        A JSON serializable dict representation of State.
        """
        occupation = (self.occupation.numerator, self.occupation.denominator)
        d = {"eigenvalue": self.eigenvalue, "irrep": self.irrep, "irrep_index": self.irrep_index,
             "occupation": occupation, "spin": self.spin, "@module": self.__class__.__module__,
             "@class": self.__class__.__name__}

        return d

    @classmethod
    def from_dict(cls, d):
        occupation = Fraction(numerator=d["occupation"][0], denominator=d["occupation"][1])
        return cls(eigenvalue=d["eigenvalue"], irrep=d["irrep"], irrep_index=d["irrep_index"],
                 occupation=occupation, spin=d["spin"])


class States(MSONable, MutableSequence):
    """
    A sequence of State, sorted in increasing order of energy.
    Describes the eigenstates of a molecule.
    """

    def __init__(self, states):
        """
        Generates the instance. The states are sorted here.

        Args:
            states (list): a list of State.
        """
        self._states = sorted(states, key=lambda s: (s.eigenvalue, s.irrep, s.irrep_index, s.spin))

    def __getitem__(self, i):
        return self._states[i]

    def __setitem__(self, i, state):
        self._states[i] = state

    def __delitem__(self, i):
        self._states.__delitem__(i)

    def __len__(self):
        return len(self._states)

    def insert(self, index, state):
        self._states.insert(index, state)

    @staticmethod
    def _generate_states_lists(eigen_data, shells, spin=None):
        """
        Helper function that generates a list of states for a specifc type
        of shells.

        Args:
            eigen_data (list): a list with the eigenvalues data. irrep, index in the
                irrep and eigenvalues extracted from the  mos files.
            shells (Shell): the shell describing the occupation of the states.
            spin (str): "a" or "b" for alpha and beta spins if calculation is uhf,
                None otherwise.

        Returns:
            (list) a list of State
        """
        states = []
        for d in eigen_data:
            irrep = d[0]
            irrep_ind = d[1]
            eig = d[2]
            occ = 0
            if irrep in shells.irreps:
                for i, (irrep_shell, index_shell) in enumerate(shells.states):
                    if irrep_shell == irrep and index_shell == irrep_ind:
                        degeneracy = irrep_size[irrep[0].upper()]
                        occ = shells.occupations[i] * degeneracy
            states.append(State(eig, irrep, irrep_ind, occ, spin))

        return states

    @classmethod
    def from_file(cls, filename="control"):
        """
        Generates the instance from a file (usually control)

        Args:
            filename (str): name of the file

        Returns:
            An instance of States
        """
        c = Control.from_file(filename)
        states = []
        # files are read directly here for performance issues. The mos files can
        # be large and the DataGroups object performs useless operations
        # in this case. We assume that the mo will always be in an external file.
        if c.is_uhf:
            alpha_shells = c.get_shells("alpha")
            beta_shells = c.get_shells("beta")
            with open(c.show_subfile_fname("uhfmo_alpha")) as f:
                alpha_eigen = get_mos_energies(f.read())
            with open(c.show_subfile_fname("uhfmo_beta")) as f:
                beta_eigen = get_mos_energies(f.read())
            states.extend(cls._generate_states_lists(alpha_eigen, alpha_shells, spin="a"))
            states.extend(cls._generate_states_lists(beta_eigen, beta_shells, spin="b"))
        else:
            shells = c.get_shells("closed")
            with open(c.show_subfile_fname("scfmo")) as f:
                eigen = get_mos_energies(f.read())
            states = cls._generate_states_lists(eigen, shells, spin=None)

        return cls(states)

    @property
    def total_electrons(self):
        """
        Total numer of electrons. Can be a float if occupations are not integer.
        """
        return np.sum(self.occupations)

    @property
    def occupations(self):
        """
        List of all the occupations.
        """
        return [s.occupation for s in self._states]

    @property
    def irreps(self):
        """
        List of all the irreducible representations
        """
        return [s.irrep for s in self._states]

    @property
    def irrep_indices(self):
        """
        List of all the indices of the irreducible representations.
        """
        return [s.irrep_index for s in self._states]

    @property
    def eigenvalues(self):
        """
        List of all the eigenvalues.
        """
        return [s.eigenvalue for s in self._states]

    @property
    def spins(self):
        """
        List of all the spins.
        """
        return [s.spin for s in self._states]

    @property
    def is_uhf(self):
        """
        True if uhf calculations (i.e. spins are not None), False otherwise.
        """
        return self._states[0].spin is not None

    def __str__(self):
        lines = [""]
        for s in self._states:
            spin = s.spin if s.spin is not None else " "
            lines.append(" {}\t{}  {}\t{}\t\t{}".format(s.irrep, s.irrep_index, spin, s.occupation, s.eigenvalue))

        return "\n".join(lines)

    @property
    def homo_index(self):
        """
        Index of the HOMO state.
        """
        occupied_indices = np.where(np.array(self.occupations) > 0)[0]
        if len(occupied_indices) == 0:
            raise RuntimeError("No occupied states present")
        return occupied_indices[-1]

    @property
    def lumo_index(self):
        """
        Index of the LUMO state.
        None if no empty state is available.
        """
        unoccupied_indices = np.where(np.array(self.occupations) == 0)[0]
        if len(unoccupied_indices) == 0:
            return None
        return unoccupied_indices[0]

    @property
    def homo(self):
        """
        HOMO state.
        """
        return self[self.homo_index]

    @property
    def lumo(self):
        """
        LUMO state.
        None if no empty state is available.
        """
        lumo_index = self.lumo_index
        if lumo_index is None:
            return None
        return self[self.lumo_index]

    @property
    def gap(self):
        """
        Gap between the HOMO and the LUMO.
        Can be negative if a hole is present.
        None if no empty state is available.
        """

        lumo = self.lumo
        if lumo is None:
            return None
        return lumo.eigenvalue - self.homo.eigenvalue

    @property
    def n_states(self):
        """
        Number of states.
        """
        return len(self)

    @property
    def has_hole(self):
        """
        True if the system has a hole.
        """
        empty = False
        for occ in self.occupations:
            if occ == 0:
                empty = True
            elif empty:
                return True

        return False

    @property
    def has_fractional_occ(self):
        """
        True if any of the states has a fractional occupation.
        """
        return any(s.has_fractional_occ for s in self._states)

    def generate_lowest_filled_states(self, allow_fractional=None, only_occupied=False,
                                      reorder_irrep_index=False):
        """
        Creates a new States object filling the states from the lowest with
        the current number of electrons.
        This procedure can lead to fractional occupations even if initially
        not present because of states degeneracies. If not allowed an exception
        will be raised in this case.

        Args:
            allow_fractional (bool): whether to allow fractional occupations in the
                states. If None it will be set equal to self.has_fractional_occ.
            only_occupied (bool): if True only occupied states will be considered.
            reorder_irrep_index (bool): if True the irrep_index of the states
                will changed to be in ascending order.

        Returns:
            States: a list of the lowest lying states filled with electrons.

        Raises:
            RuntimeError: if the states will contain fractional occupations
                and allow_fractional is False.
        """
        if allow_fractional is None:
            allow_fractional = self.has_fractional_occ

        remaining_electrons = self.total_electrons
        states = []

        irrep_ind = defaultdict(lambda: defaultdict(int))
        for s in self._states:
            new_s = copy(s)
            if reorder_irrep_index:
                irrep_ind[s.spin][s.irrep] += 1
                new_s.irrep_index = irrep_ind[s.spin][s.irrep]
            max_occ = s.max_occupation
            if remaining_electrons <= 0:
                new_s.occupation = Fraction(0)
            elif max_occ > remaining_electrons:
                if not allow_fractional:
                    raise RuntimeError("Fractional occupancies of states are not allowed.")
                new_s.occupation = Fraction(remaining_electrons)
            else:
                new_s.occupation = Fraction(max_occ)

            states.append(new_s)

            remaining_electrons -= max_occ
            if only_occupied and remaining_electrons <= 0:
                break

        return States(states)

    def filter_states(self, spin=None, irrep=None):
        """
        Generates a States instance with a subset of the State objects filtered
        by spin and irrep.
        Notice that the same instances of State as in the current object wiil
        be used.
        Args:
            spin (str): string describing the spin. Values allowed: "alpha",
                "beta" and the shortcuts "a" and "b".
            irrep (str): the symbol of the irreducible representations.

        Returns:
            States: the filtered list of states.

        Raises:
            ValueError: if spin is requested for a non uhf calculation.
        """
        if spin and not self.is_uhf:
            raise ValueError("The states do not include spin.")

        states = []
        for s in self._states:
            if spin is None or spin[0] == s.spin:
                if irrep is None or irrep == s.irrep:
                    states.append(s)

        return States(states)

    def get_shells(self, spin=None):
        """
        Generates a Shells object with the occupied shells in the list of States.

        Args:
            spin (str): string describing the spin. Values allowed: "alpha",
                "beta" and the shortcuts "a" and "b".

        Returns:
            Shells: the occupied shells for the states
        """
        if (spin is not None) != self.is_uhf:
            raise ValueError("incompatible configuration and spin request")
        shell_states = []
        occupations = []
        for s in self._states:
            if s.occupation == 0:
                continue
            if spin is None or spin[0] == s.spin:
                shell_states.append((s.irrep, s.irrep_index))
                occupations.append(s.occupation / irrep_size[s.irrep[0].upper()])

        from turbomoleio.core.control import Shells
        return Shells(shell_states, occupations)

    def as_dict(self):
        """
        A JSON serializable dict representation of States.
        """
        d = {"states": [s.as_dict() for s in self._states],
             "@module": self.__class__.__module__,
             "@class": self.__class__.__name__}
        return d


class EigerOutput:
    """
    Class to read and store the output from the eiger command.
    """

    def __init__(self, eigenvalues, irreps, irrep_indices, occupations, spin, gap, nelec):
        """
        The lists should match among them and be sorted in ascending values with
        respect to the eigenvalues.
        Args:
            eigenvalues (list): the eigenvalues.
            irreps (list): the irreducible representations.
            irrep_indices (list): the inddices for each irreducible representation.
            occupations (list): the occupations.
            spin (list): the spins.
            gap (float): the value of the gap.
            nelec (float): the total number of electrons.
        """
        self.eigenvalues = eigenvalues
        self.irreps = irreps
        self.irrep_indices = irrep_indices
        self.occupations = occupations
        self.spin = spin
        self.gap = gap
        self.nelec = nelec

    @classmethod
    def from_string(cls, string):
        """
        Generates the object from the string of the output from eiger.

        Args:
            string (str): the string of the output of eiger.

        Returns:
            An instance of EigerOutput
        """
        eigenvalues = []
        irreps = []
        irrep_indices = []
        occupations = []
        spin = []
        gap = None
        nelec = 0
        parse_states = False
        r = re.compile(r"^\s*\d+.\s+([ab]*)\s+(\d+)\s+([^\s]*)\s+([\d.+\-]*)\s+([\d.+\-]+)\s+H")
        for l in string.splitlines():
            if "Gap" in l:
                gap = float(l.split()[2])
            if "Electrons=" in l:
                nelec = float(l.split()[5].replace(",", ""))
            if parse_states:
                if not l.strip():
                    break
                m = r.search(l)
                if not m:
                    raise RuntimeError("Could not match regex for line: {}".format(l))
                s = m.group(1) if m.group(1) else None
                spin.append(s)
                irrep_indices.append(int(m.group(2)))
                irreps.append(m.group(3))
                occ = float(m.group(4)) if m.group(4) else 0
                occupations.append(occ)
                eigenvalues.append(float(m.group(5)))
            if "Nr. " in l:
                parse_states = True

        return cls(eigenvalues, irreps, irrep_indices, occupations, spin, gap, nelec)

    @classmethod
    def from_file(cls, filename):
        """
        Generates the object from a file with the output from eiger.

        Args:
            filename (str): the name of the file containing the output of eiger.

        Returns:
            An instance of EigerOutput
        """
        with open(filename) as f:
            return cls.from_string(f.read())

    def compare_states(self, states, tol=1e-6):
        """
        Compares the current values with a States object and checks that
        the two are equivalent. Number are checked within the specified
        tolerance.

        Args:
            states (States): States object to be compared
            tol (float): tolerance on the floating numbers

        Returns:
            None if the data are equivalent, otherwise a string describing the
            difference between this instance and the States object.
        """
        if len(self.eigenvalues) != len(states):
            return "number of states"

        if abs(self.gap - states.gap) > tol:
            return "gap"

        if abs(self.nelec - states.total_electrons) > 0.01:
            return "number of electrons"

        for sp, irrep, ind, occ, eig in zip(self.spin, self.irreps, self.irrep_indices,
                                           self.occupations, self.eigenvalues):
            for s in states:
                if s.irrep == irrep and s.irrep_index == ind and s.spin == sp:
                    if abs(s.occupation - occ) > tol:
                        return "occupation for state {}".format(s)
                    if abs(s.eigenvalue - eig) > tol:
                        return "eigenvalue for state: {}".format(s)
                    break
            else:
                return "no match for values: {} {} {} {} {}".format(sp, irrep, ind, occ, eig)

        return None


class EigerRunner:
    """
    Class that runs the eiger executable,
    """

    def __init__(self, executable="eiger", data_path="."):
        """
        Args:
            executable (str): the eiger executable.
            data_path (str): path to where the TM data are stored.
        """
        self.executable = executable
        self.data_path = data_path
        self.out = None

    def run(self):
        """
        Runs eiger in the data_path directory

        Returns:
            None
        """
        with cd(self.data_path):
            result = subprocess.run([self.executable, "-a"], stdout=subprocess.PIPE)
            self.out = result.stdout.decode('utf-8')

    def get_eiger_output(self):
        """
        Generates an EigerOutput instance from the output of eiger

        Returns:
            EigerOutput

        Raises:
            ValueError: if eiger has not run
        """
        if not self.out:
            raise ValueError("No output")

        return EigerOutput.from_string(self.out)

    def to_file(self, filepath):
        """
        Writes the output of eiger to a file.

        Args:
            filepath (str): path to the file.

        Returns:
            None

        Raises:
            ValueError: if eiger has not run
        """
        if not self.out:
            raise ValueError("No output")

        with open(filepath, "w") as f:
            f.write(self.out)
