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

"""
This module contains the object and functions to interact directly with a
control file. Based on the DataGroups class.
"""

import re
import os
import shutil
import numpy as np
from fractions import Fraction
from collections import defaultdict

from monty.json import MSONable
from monty.os import makedirs_p, cd
from pymatgen.util.plotting import add_fig_kwargs, get_ax_fig_plt
from turbomoleio.core.datagroups import DataGroups
from turbomoleio.core.symmetry import irrep_size
from turbomoleio.core.molecule import get_mol_and_indices_frozen, MoleculeSystem


class Defaults:
    """
    Class collecting default values used by TM executables when not
    present in the control file.
    """

    METRIC = 3


class Energy(MSONable):
    """
    Represents the "energy" data group of TurboMole.
    Usually stored in the "energy" file or directly in the control file.
    """

    def __init__(self, scf=None, scfkin=None, scfpot=None, mp2=None):
        """

        Args:
            scf (list): energies from scf calculations.
            scfkin (list): kinetic contribution to the scf energy.
            scfpot (list): potential contribution to the scf energy.
            mp2 (list): mp2 contribution to the total energy (not counted in scf).
        """
        self.scf = np.array(scf)
        self.scfkin = np.array(scfkin)
        self.scfpot = np.array(scfpot)
        if mp2:
            self.mp2 = np.array(mp2)
            self.total = self.scf + self.mp2
        else:
            self.mp2 = None
            self.total = np.array(self.scf)

    @classmethod
    def from_string(cls, string):
        """
        Creates Energy object reading from a string.

        Args:
            string (str): the string of the "energy" datagroup.
        """
        energy_dg_sp = string.splitlines()
        header_sp = energy_dg_sp[0].split()
        h_ind = {"SCF": None, "SCFKIN": None, "SCFPOT": None, "MP2": None}
        for h in h_ind.keys():
            try:
                h_ind[h] = header_sp.index(h) + 1
            except:
                pass
        scf = []
        scfkin = []
        scfpot = []
        mp2 = [] if h_ind["MP2"] is not None else None
        for l in energy_dg_sp[1:]:
            lsp = l.split()
            if not lsp:
                continue
            scf.append(float(lsp[h_ind["SCF"]]))
            scfkin.append(float(lsp[h_ind["SCFKIN"]]))
            scfpot.append(float(lsp[h_ind["SCFPOT"]]))
            if mp2 is not None:
                mp2.append(float(lsp[h_ind["MP2"]]))
        return cls(scf=scf, scfkin=scfkin, scfpot=scfpot, mp2=mp2)

    @classmethod
    def from_file(cls, filename='energy'):
        """Creates Energy object reading from a given file.

        Args:
            filename (str): Name of the file from which this Energy object
                should be read, default is "energy".
        Raises:
            RuntimeError: if the energy file comes from a subfile and is missing.
        """
        dg = DataGroups.from_file(filename)
        string = dg.show_data_group("energy", show_from_subfile=True, raise_if_missing_subfile=True)
        return cls.from_string(string)

    @property
    def n_steps(self):
        """
        The number of steps present.
        """
        return len(self.scf)

    @property
    def delta_e(self):
        """
        The variation of the total energy in the last step.
        None if less than 2 steps presents.
        """
        if self.n_steps < 2:
            return None

        return abs(self.total[-1] - self.total[-2])

    @add_fig_kwargs
    def plot(self, ax=None, **kwargs):
        """
        Plot the evolution of the energies.

        Args:
            ax (Axes):  a matplotlib Axes or None if a new figure should be created.
            kwargs: arguments passed to matplotlib plot function.

        Returns:
            A matplotlib Figure
        """

        ax, fig, plt = get_ax_fig_plt(ax=ax)

        ax.plot(range(self.n_steps), self.total, **kwargs)
        ax.set_xlabel("Steps")
        ax.set_ylabel("Energy")

        return fig


class Gradient(MSONable):
    """
    Represents the "grad" data group of TurboMole.
    Usually stored in the "gradient" file or directly in the control file.
    """

    def __init__(self, gradients=None, energies=None, molecules=None):
        """
        Args:
            gradients (list): gradients for all the steps and all the atoms
                with shape nsteps*natoms*3.
            energies (list): total energies for each step.
            molecules (list): list of MoleculeSystem with coordinates for each step.
        """
        self.gradients = np.array(gradients)
        self.energies = np.array(energies)
        self.molecules = molecules

    @classmethod
    def from_string(cls, string):
        """
        Creates Gradient object reading from a given file.

        Args:
            string (str): the string of the "grad" datagroup.
        """
        # skip the first as it is before the first
        cycles = string.split("cycle =")[1:]

        gradients = []
        energies = []
        molecules = []

        for c in cycles:
            lines = c.splitlines()
            header = lines[0]
            match = re.search(r"energy\s+=\s+([+-]?[0-9]*[.]?[0-9]+)", header)
            energies.append(float(match.group(1)))
            coordinates = []
            grad = []
            for l in lines[1:]:
                l = l.replace("D", "E")
                lsp = l.split()
                if not lsp:
                    continue

                # can be with or without frozen indices
                if len(lsp) == 4 or len(lsp) == 5:
                    coordinates.append(l)
                elif len(lsp) == 3:
                    grad.append([float(g) for g in lsp])
                else:
                    raise RuntimeError("Encountered line with unexpected number of tokens: {}".format(l))

            gradients.append(grad)
            mol, fi = get_mol_and_indices_frozen("\n".join(coordinates))
            molecules.append(MoleculeSystem(molecule=mol, frozen_indices=fi))

        return cls(gradients=gradients, energies=energies, molecules=molecules)

    @classmethod
    def from_file(cls, filename='gradient'):
        """
        Creates Gradient object reading from a given file.

        Args:
            filename (str): Name of the file from which this Gradient object
                should be read, default is "gradient".
        """
        dg = DataGroups.from_file(filename)
        string = dg.show_data_group("grad", show_from_subfile=True, raise_if_missing_subfile=True)
        return cls.from_string(string)

    @property
    def norms(self):
        """
        Total norm of the gradients for each step.
        """
        return np.linalg.norm(self.gradients, axis=(1,2))

    @property
    def max_gradients(self):
        """
        Norms of the largest gradients for each step.
        """
        return np.max(np.linalg.norm(self.gradients, axis=(2)), axis=1)

    @property
    def n_steps(self):
        """
        The number of steps present.
        """
        return len(self.gradients)

    @property
    def last_grad_norm(self):
        """
        The value of the last gradient norm.
        None if no steps presents.
        """
        if self.n_steps == 0:
            return None

        return self.norms[-1]

    @property
    def last_grad_max(self):
        """
        The value of the maximum value of the gradient in the last step.
        None if no steps presents.
        """
        if self.n_steps == 0:
            return None

        return self.max_gradients[-1]

    @add_fig_kwargs
    def plot(self, **kwargs):
        """
        Plot the evolution of the maximum and the norm of the gradient.

        Args:
            kwargs: arguments passed to matplotlib plot function.

        Returns:
            A matplotlib Figure
        """

        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.subplots(2, 1, sharex=True)

        ax[0].plot(range(self.n_steps), self.norms, **kwargs)
        ax[0].set_ylabel("Gradient norm")

        ax[1].plot(range(self.n_steps), self.max_gradients, **kwargs)
        ax[1].set_xlabel("Steps")
        ax[1].set_ylabel("Gradient max")

        return fig


class Shells(MSONable):
    """
    Object describing the occupation of a type of shells.
    It can map $closed shells, $alpha shells and $beta shells
    ($open shells not supported).
    Handles integer, fractional and floating point occupations.
    Occupations are stored as Fraction objects
    Supports presence of the same irrep in multiple lines.
    """

    def __init__(self, states, occupations):
        """
        Args:
            states (list): list of tuples (irrep, irrep_index) with the occupied
                states for each.
            occupations (list): list of Fraction with the occupations for each state.
        """
        self.states = states
        self.occupations = occupations

    @classmethod
    def from_string(cls, string):
        """
        Generates an instance of Shells from a string.
        Examples of a string:
           a1      1-3                      ( 2 )
           b1      1,3                       ( 2 )
           b2      1                        ( 2 )

        Args:
            string (str): the string of the datagroup.

        Returns:
            Shells
        """
        states_per_irrep = {}
        states = []
        occupations = []

        # for the occupations matches integers, fractions (e.g.: 1 / 3)
        # and floats (e.g.: 0.833473)
        r = r"^\s*([\w\'\"]+)\s+([\d\-\s,]+)\s+\(([.\s\d/]+)\)"

        for l in string.splitlines():
            l = l.strip()
            if not l:
                continue

            match = re.search(r, l)
            if match is None:
                raise RuntimeError("Could not match occupation regex for line: {}".format(l))

            irrep = match.group(1)
            # states could be a string of the form "1-3,5-8,12" with "-" defining intervals
            # and "," separating groups.
            states_str = match.group(2)
            irrep_occupation = Fraction(match.group(3).replace(" ", ""))
            states_per_irrep = []
            for s in states_str.split(","):
                if "-" in s:
                    split = s.split("-")
                    states_per_irrep.extend(range(int(split[0]), int(split[1]) + 1))
                else:
                    states_per_irrep.append(int(s))

            for s in states_per_irrep:
                states.append((irrep, s))
                occupations.append(irrep_occupation)

        return cls(states, occupations)

    @classmethod
    def from_file(cls, shells_type, filename="control"):
        """
        Generates an instance of Shells from a file containing
        the datagroup.

        Args:
            shells_type (str): the type of shell. The supported type are
                "closed", "alpha" or "beta".
            filename (str): the name of the file.

        Returns:
            Shells

        Raises:
            ValueError: if the shell of the specified type is not present
        """

        dg = DataGroups.from_file(filename)
        string = dg.show_data_group("{} shells".format(shells_type))
        if not string:
            raise ValueError("No shells of type {} if file {}".format(shells_type, filename))

        return cls.from_string(string)

    @property
    def irreps(self):
        """
        Set of irreducible representations.
        """
        return set(s[0] for s in self.states)

    @property
    def total_electrons(self):
        """
        Total number of electrons present in all the shells.

        Returns:
            Fraction: the total number of electrons.
        """
        tot = 0
        for s, o in zip(self.states, self.occupations):
            irrep = s[0]
            degeneracy = irrep_size[irrep[0].upper()]
            tot += o * degeneracy

        return tot

    def to_datagroup(self):
        """
        Generates a string that can be used to set the datagroup in control.

        Returns:
             The string of the data group.
        """
        lines = []

        dg_dict = defaultdict(lambda: defaultdict(list))
        for state, occ in zip(self.states, self.occupations):
            dg_dict[state[0]][occ].append(state[1])

        # sort to make output reproducible and more easily readable
        for irrep in sorted(dg_dict.keys()):
            for occ, states in sorted(dg_dict[irrep].items(), key=lambda x: np.min(x[1])):

                states = sorted(states)
                # add one fictitious state to trigger the addition of the last step
                states.append(states[-1] + 2)
                # get the string representations of the groups
                groups = []
                start = states[0]
                for i in range(1, len(states)):
                    if states[i] - states[i - 1] > 1:
                        end = states[i - 1]
                        if start == end:
                            groups.append(str(start))
                        else:
                            groups.append("{}-{}".format(start, end))
                        start = states[i]
                group_str = ",".join(groups)
                if occ.denominator == 1:
                    str_occ = str(occ.numerator)
                else:
                    str_occ = "{} / {}".format(occ.numerator, occ.denominator)
                lines.append(" {0}      {1:<37} ( {2} )".format(irrep, group_str, str_occ))

        return "\n" + "\n".join(lines)

    def as_dict(self):
        """
        A JSON serializable dict representation of Shells.
        """
        occupations = [(o.numerator, o.denominator) for o in self.occupations]
        d = dict(states=self.states, occupations=occupations)

        return d

    @classmethod
    def from_dict(cls, d):
        """
        Generates an instance of Shell from a JSON serialized representation.

        Args:
            d (dict): the dictionary with the data to initialize Shells

        Returns:
             A Shells instance
        """
        occupations = [Fraction(numerator=o[0], denominator=o[1]) for o in d["occupations"]]
        return cls(states=d["states"], occupations=occupations)


class Control(DataGroups):
    """Represents the control string/file for TurboMole.

    The `Control` class is used to store the parameters of a control file.
    Parameters in a control file are stored as "Data Groups". A data group
    is a key-value pair storing data "value" for the specific "key". Value may
    be a single number, a string or a series of strings, numbers or
    combinations thereof.
    """

    def __init__(self, string=None, dg_list=None):
        """Initializes a `Control` class.

        Args:
            string: The control string (in the data group format).
                Usually read from file.
            dg_list: List of the data groups in the control string/file.
        """
        super(Control, self).__init__(string=string, dg_list=dg_list)

    @classmethod
    def from_file(cls, filename='control'):
        """Creates Control object reading from a given file.

        Args:
            filename (str): Name of the file from which this Control object
                should be read, default is "control".
        """
        return super().from_file(filename=filename)

    def to_file(self, filename='control'):
        """Writes this Control object to a file.
        If the file exists it will be overwritten.

        Args:
            filename (str): Name of the file to which this Control object
                should be written.
        """
        super(Control, self).to_file(filename=filename)

    def add_cosmo(self, epsilon=None, nppa=None, nspa=None, disex=None, rsolv=None, routf=None, cavity=None,
                  use_old_amat=None):
        """
        Adds the $cosmo datagroup to the control, using the keywords to set the cosmo options.
        If None all the options will not be written to the control file.

        Args:
            epsilon (float): permittivity used for scaling of the screening charges
            nppa (int): number of basis grid points per atom
            nspa (int): number of segments per atom
            disex (float): distance threshold for A matrix elements (Angstrom)
            rsolv (float): distance to outer solvent sphere for cavity construction (Angstrom)
            routf (float): factor for outer cavity construction in the outlying charge correction
            cavity (str): acceptable values are "open" (leave untidy seams between atoms) and
                "closed" (pave intersection seams with segments)
            use_old_amat (bool): if True adds the ``use_old_amat`` to the ``$cosmo`` datagroup,
                i.e. uses A matrix setup of TURBOMOLE 5.7.

        Returns:
            None
        """

        cosmo_lines = [""]

        if epsilon is not None:
            cosmo_lines.append("   epsilon={}".format(epsilon))
        if nppa is not None:
            cosmo_lines.append("   nppa={}".format(nppa))
        if nspa is not None:
            cosmo_lines.append("   nspa={}".format(nspa))
        if disex is not None:
            cosmo_lines.append("   disex={}".format(disex))
        if rsolv is not None:
            cosmo_lines.append("   rsolv={}".format(rsolv))
        if routf is not None:
            cosmo_lines.append("   routf={}".format(routf))
        if cavity is not None:
            cosmo_lines.append("   cavity={}".format(cavity))
        if use_old_amat is True:
            cosmo_lines.append("   use_old_amat")

        self.cdg("$cosmo", "\n".join(cosmo_lines))

        self.cdg("$cosmo_out", ' = out.cosmo')

    @classmethod
    def from_metric(cls, metric):
        """
        Generates a new Control containing the value of the metric in $redund_inp
        and an empty $coord.

        Args:
            metric (int): the value of the metric
        """
        metric_str = "$redund_inp\n    metric {}\n$coord\n$end".format(metric)
        return cls(string=metric_str)

    @property
    def energy(self):
        """
        An Energy object with the data from the $energy data group.
        None if the datagroup is absent/empty or if the file is missing.
        """
        energy_data_block = self.show_data_group('$energy', show_from_subfile=True, default="")
        if not energy_data_block.strip() or "file=" in energy_data_block:
            return None
        return Energy.from_string(string=energy_data_block)

    @property
    def gradient(self):
        """
        A Gradient object with the data from the $grad data group.
        None if the datagroup is absent/empty or if the file is missing.
        """
        grad_data_block = self.show_data_group('$grad', show_from_subfile=True, default="")
        if not grad_data_block.strip() or "file=" in grad_data_block:
            return None
        return Gradient.from_string(string=grad_data_block)

    def set_disp(self, dispersion_correction):
        """
        Sets the dispersion correction.

        Args:
            dispersion_correction (str): the name of the method used for the dispersion correction.
                Supported values are "DFT-D1", "DFT-D2" "DFT-D3" and "DFT-D3 BJ", case insensitive.
                If None all the "disp" datagroup will be eliminated.

        Returns:
            None
        """

        for dg in ("olddisp", "disp", "disp3"):
            self.kill_data_group(dg, strict=True)

        if not dispersion_correction:
            return

        dispersion_correction = dispersion_correction.lower()

        if dispersion_correction == "dft-d1":
            self.add_data_group("olddisp", "")
        elif dispersion_correction == "dft-d2":
            self.add_data_group("disp", "")
        elif dispersion_correction == "dft-d3":
            self.add_data_group("disp3", "")
        elif dispersion_correction == "dft-d3 bj":
            self.add_data_group("disp3", "bj")
        else:
            raise ValueError("Unsupported value for dispersion_correction: {}".format(dispersion_correction))

    def remove_last_energy(self, filename="control", backup_suffix=None):
        """
        Removes last energy step from the energy data group and writes it to
        the correct file, depending on the original location of the data group.
        N.B. if the energies are in control and not in an external file, this
        will call the to_file of the current instance of Control. If other
        modifications have been done to the dg_list will be written to the
        file as well.

        Args:
            filename (str): name of the file to which the current instance of
                Control will be written in case the energies are control.
                If energies are in external file it will be ignored.
            backup_suffix (str): if not None a file named filename+backup_suffix
                will be created with the current content of Control.

        Raises:
            RuntimeError: if there is no energy step to be removed or if energy
                file is missing.
        """
        energy_dg =  self.show_data_group("energy", default="", show_from_subfile=True)

        energy_dg_sp = energy_dg.split("\n")
        if len(energy_dg_sp) < 3:
            raise RuntimeError('No energy to remove in the "energy" data group.')
        energy_dg_sp.pop(-2)
        new_energy_dg_sp = '\n'.join(energy_dg_sp)

        block_data = self.show_data_group("energy", default="", show_from_subfile=False)

        dg_filename = self._get_subfile_fname(block_data, raise_if_regular_and_subfile=False)

        if dg_filename:
            dg = DataGroups.from_file(dg_filename)
            dg.change_data_group("energy", new_energy_dg_sp)
            if backup_suffix:
                shutil.copy2(dg_filename, dg_filename+backup_suffix)
            dg.to_file(dg_filename)
        elif not filename:
            raise ValueError("The datagroup is in the current file. An output filename"
                             "should be specified.")
        else:
            self.change_data_group("energy", new_energy_dg_sp)
            if backup_suffix:
                shutil.copy2(filename, filename+backup_suffix)
            self.to_file(filename)

    def remove_last_gradient(self, filename="control", backup_suffix=None):
        """
        Removes last gradient step from the grad data group and writes it to
        the correct file, depending on the original location of the data group.
        N.B. if the gradients are in control and not in an external file, this
        will call the to_file of the current instance of Control. If other
        modifications have been done to the dg_list will be written to the
        file as well.

        Args:
            filename (str): name of the file to which the current instance of
                Control will be written in case the gradients are control.
                If gradients are in external file it will be ignored.
            backup_suffix (str): if not None a file named filename+backup_suffix
                will be created with the current content of Control.

        Raises:
            RuntimeError: if there is no gradient step to be removed or if
                gradient file is missing.
        """
        grad_dg =  self.show_data_group("grad", default="", show_from_subfile=True)

        if "cycle" not in grad_dg:
            raise RuntimeError('No gradient to remove in the "grad" data group.')

        # take the string until before the last "cycle"
        m = re.match(r"(.+)\n\s*cycle.+$", grad_dg, re.DOTALL)
        new_grad_dg = m.group(1)

        block_data = self.show_data_group("grad", default="", show_from_subfile=False)

        dg_filename = self._get_subfile_fname(block_data, raise_if_regular_and_subfile=False)

        if dg_filename:
            dg = DataGroups.from_file(dg_filename)
            dg.change_data_group("grad", new_grad_dg)
            if backup_suffix:
                shutil.copy2(dg_filename, dg_filename+backup_suffix)
            dg.to_file(dg_filename)
        elif not filename:
            raise ValueError("The datagroup is in the current file. An output filename"
                             "should be specified.")
        else:
            self.change_data_group("grad", new_grad_dg)
            if backup_suffix:
                shutil.copy2(filename, filename+backup_suffix)
            self.to_file(filename)

    def get_shells(self, shells_type):
        """
        Extract a Shell object for a specific type of shell.
        Args:
            shells_type: the type of shell. Can be "closed", "alpha" or "beta".

        Returns:
            An instance of Shells
        """
        string = self.show_data_group("{} shells".format(shells_type))
        if not string:
            raise ValueError("No shells of type {}".format(shells_type))

        return Shells.from_string(string)

    @property
    def is_uhf(self):
        """
        True if the $uhf datagroup is present in control.
        """
        return self.show_data_group("uhf") is not None

    def get_charge(self):
        """
        Extract the charge from the "$charge from xxx" datagroup.

        Returns:
            float: the charge in the control file. None if not present.
        """
        dg = self.show_data_group("charge from")
        if not dg:
            return None

        return float(dg.splitlines()[1].split()[0])

    def get_subfiles_list(self):
        """
        Extracts the list of files referenced in the datagroups with
        "file=xxx" and "file xxx".

        Returns:
            list: list of files referenced in the datagroups.
        """
        files = set()
        for dg in self.dg_list:
            for l in dg.splitlines():
                m = re.search(" file[ =](.*)", l)
                if m:
                    files.add(m.group(1).strip())

        return list(files)

    def cpc(self, dest_dir, force_overwrite=False):
        """
        Copies the control file and all files referenced here.
        Creates the destination folder if not already existing.

        Args:
            dest_dir (str): path to the destination folder.
            force_overwrite (bool): if True files already present in the
                destination folder will be overwritten.
        """
        dest_dir = os.path.abspath(dest_dir)
        makedirs_p(dest_dir)
        dest_control_path = os.path.join(dest_dir, "control")
        if not os.path.isfile(dest_control_path) or force_overwrite:
            self.to_file(dest_control_path)

        for fn in self.get_subfiles_list():
            dest_file_path = os.path.join(dest_dir, fn)
            # copy only if the file exists
            if os.path.isfile(fn) and (not os.path.isfile(dest_file_path) or force_overwrite):
                shutil.copy2(fn, dest_file_path)


def kdg(data_group, filepath="control", strict=True, backup_file=None):
    """
    Function that removes `data_group` from the control file specified by the path.
    Comments and empty lines will be removed. The file will be overwritten.

    Args:
        data_group (str): Data group to be removed.
        filepath (str): path to the control file
        strict (bool): If True `data_group` should be an exact match.
            If False, any Data group starting with `data_group` will be removed.
        backup_file (str): filename to which the original file will be backed up.
            If None no backup will be created.

    Returns:
        None
    """

    if backup_file:
        shutil.copy2(filepath, backup_file)
    c = Control.from_file(filepath)
    c.kdg(data_group, strict)
    c.to_file(filepath)


def sdg(data_group, filepath="control", default=None,
        show_from_subfile=True,
        raise_if_multiple_subfiles=False,
        raise_if_missing_subfile=False,
        raise_if_regular_and_subfile=False):
    """
    Function that shows the `data_group` from the control file specified by the path.
    Comments and empty lines will be removed.
    If the filepath is a full path (i.e. not just a file name), it will change
    directory to the specified one, so that the subfiles are correctly read.

    Args:
        data_group (str): Data group to be removed.
        filepath (str): path to the control file
        default (str): the default value that will be returned if the
                data_group is not present in the list. Default is None.
        show_from_subfile (bool): If True, will show `data_group` from within the
            "subfile" if the data block contains a "file=FILENAME". If
            False the "file=FILENAME" block is simply returned. This
            supposes that file exists in the current directory.
        raise_if_multiple_subfiles (bool): Whether to raise an error if
            multiple "file=" directives are present in this data group. If
            False, just returns the standard control data block.
        raise_if_missing_subfile (bool): Whether to raise an error if the
                subfile does not exist in the current directory or exists but
                is empty. If False, just returns the standard control data block.
        raise_if_regular_and_subfile (bool): Whether to raise an error if
            data group contains both a reference to a file with a
            "file=FILENAME" and regular data block options. If
            False, just returns the standard control data block.
    Returns:
        str: the value of the selected datagroup. None if `no data_group` is found.
    """

    # cd to the directory so that the referenced subfiles are working
    dirname = os.path.dirname(filepath) or os.getcwd()
    with cd(dirname):
        c = Control.from_file(filepath)
        str_dg = c.sdg(data_group, default=default,
                       show_from_subfile=show_from_subfile,
                       raise_if_multiple_subfiles=raise_if_multiple_subfiles,
                       raise_if_missing_subfile=raise_if_missing_subfile,
                       raise_if_regular_and_subfile=raise_if_regular_and_subfile)

    return str_dg


def adg(data_group, data_block, filepath="control", backup_file=None):
    """
    Function that adds `data_group` with `data_block` from the control file specified by the path.
    Comments and empty lines will be removed. The file will be overwritten.

    Args:
        data_group (str): Data group (key) to be added. The dollar ("$")
            sign will be automatically added if not present.
        data_block (str): Data block corresponding to the data group.
        filepath (str): path to the control file
        backup_file (str): filename to which the original file will be backed up.
            If None no backup will be created.

    Returns:
        None

    Raises:
            RuntimeError: if data_group already exists.
    """

    if backup_file:
        shutil.copy2(filepath, backup_file)
    c = Control.from_file(filepath)
    c.adg(data_group, data_block)
    c.to_file(filepath)


def cdg(data_group, data_block, filepath="control", backup_file=None):
    """
    Function that replaces an existing `data_group` with `data_block` from the control file specified
    by the path. Comments and empty lines will be removed. The file will be overwritten.

    If `data_block` is None it is equivalent to kdg. To set a data group with no
    explicit value use and empty string for `data_block` (e.g. to add $uhf the input
    should be data_group="uhf" and data_block="").

    Args:
        data_group (str): Data group (key) to be added. The dollar ("$")
            sign will be automatically added if not present.
        data_block (str): Data block corresponding to the data group.
            If None it will simply kill_data_group for the specified data_group value.
        filepath (str): path to the control file
        backup_file (str): filename to which the original file will be backed up.
            If None no backup will be created.
    Returns:
        None

    """

    if backup_file:
        shutil.copy2(filepath, backup_file)
    c = Control.from_file(filepath)
    c.cdg(data_group, data_block)
    c.to_file(filepath)


def mdgo(data_group, options, filepath="control", backup_file=None):
    """
    Function that, given a data group that allows several options on separate
    line (e.g. $dft), updates the values of the options according to the
    dictionary provided.
    The option dictionary should have the form
    {"option_name1": "option_name1 option_value",
     "option_name2": "option_name2=option_value"}.
    The key will be used to identify the line to be modified and that line
    would be entirely replaced by the value.
    Since different options may be defined in different ways, no
    attempt is made here to identify the suitable format for the option.
    It is responsibility of the caller to specify the line for the option
    in the correct format.

    If the datagroup is not present will be created with the specified options.

    If the entire data group should be modified it would be safer to use
    change_data_group.

    Args:
        data_group (str): Data group (key) to be changed. The dollar ("$")
            sign will be automatically added if not present.
        options (dict): The options that should be added or updated. If the
            value of one key is None the option will be removed if present.
        filepath (str): path to the control file
        backup_file (str): filename to which the original file will be backed up.
            If None no backup will be created.
    Returns:
        None

    """

    if backup_file:
        shutil.copy2(filepath, backup_file)
    c = Control.from_file(filepath)
    c.mdgo(data_group, options)
    c.to_file(filepath)


def sdgo(data_group, option, filepath="control", default=None):
    """
    Given a data group that allows several options on separate line
    (e.g. $dft), returns the value of the option provided.

    Since each value may be defined in a different way, the returned
    value will include anything after the option specified. It is up
    to the caller determine if there are symbols like "=" that should
    be removed, depending on the data group and option that is queried.

    Args:
        data_group (str): Data group (key) to be changed. The dollar ("$")
            sign will be automatically added if not present.
        option (str): The options that should be returned.
        filepath (str): path to the control file.
        default (str): the default value that will be returned if the
            data_group or the option are not present in the list.

    Returns:
        str: The value of the option.
    """

    c = Control.from_file(filepath)
    return c.sdgo(data_group, option=option, default=default)


def cpc(dest_dir, force_overwrite=False, control_dir=None):
    """
    Copies a control file from the specified folder and all files
    referenced there to the destination folder.
    Creates the destination folder if not already existing.

    Args:
        dest_dir (str): path to the destination folder.
        force_overwrite (bool): if True files already present in the
            destination folder will be overwritten.
        control_dir (str): path to the directory containing the control file
            and the other files to be copied. If None the current dir will
            be used.
    """
    if control_dir is None:
        control_dir = os.getcwd()

    dest_dir = os.path.abspath(dest_dir)

    with cd(control_dir):
        c = Control.from_file("control")
        c.cpc(dest_dir=dest_dir, force_overwrite=force_overwrite)
