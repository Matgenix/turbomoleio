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
Module with basic objects extracted from the stdout of Turbomole executables.
Rely on the Parser object to extract the data.
"""
import abc
import pprint
import numpy as np
import pandas as pd
from turbomoleio.output.parser import Parser
from pymatgen.core.structure import Molecule
from pymatgen.core.units import bohr_to_ang
from pymatgen.util.plotting import add_fig_kwargs, get_ax_fig_plt
from monty.json import MSONable


class BaseData(MSONable, abc.ABC):
    """
    Base class for the data extracted from the parser of the output files.
    """

    @classmethod
    def from_file(cls, filepath):
        """
        Generates an instance of the class from a file containing the stdout
        of a Turbomole executable.

        Args:
            filepath (str): path to the file.

        Returns:
            An instance of the class.
        """
        try:
            with open(filepath, 'r') as f:
                return cls.from_string(f.read())
        except UnicodeDecodeError:
            with open(filepath, 'r', errors='ignore') as f:
                return cls.from_string(f.read())

    @classmethod
    def from_string(cls, string):
        """
        Generates an instance of the class from a string containing the stdout
        of a Turbomole executable.

        Args:
            string (str): the string with output.

        Returns:
            An instance of the class.
        """
        p = Parser(string)
        return cls.from_parser(p)

    @classmethod
    @abc.abstractmethod
    def from_parser(cls, parser):
        """
        Generates an instance of the class from a parser based on the stdout
        of a Turbomole executable. Should return None if no data could be parsed.

        Subclasses should use the parser to extract the relevant information.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            An instance of the class.
        """
        pass

    def __str__(self):
        """
        String representation of the object with print of the as_dict()
        """
        string = "Content of the {} object".format(self.__class__.__name__)
        string += pprint.pformat(self.as_dict())

        return string


class CosmoData(BaseData):
    """
    Information about the cosmo run. Inputs and outputs.
    """

    def __init__(self, info=None, parameters=None, screening_charge=None,
                 energies=None, element_radius=None):
        """
        Args:
            info (dict): initial info provided in the output by COSMO. Contains the
                keys area and volume.
            parameters (dict): COSMO parameters. Contains the keys nppa, nspa, nsph,
                npspher, disex, disex2, rsolv, routf, phsran, ampran, cavity, epsilon,
                refind, fepsi.
            screening_charge (dict): values of the screening charge. Contains
                the keys cosmo, correction, total.
            energies (dict): total and dielectric energies. Contains the keys total_energy,
                total_energy_oc_corr, dielectric_energy, dielectric_energy_oc_corr.
            element_radius (dict): to an Element object correspond values of the
                radius and a list of sites of the corresponding element in the
                molecule.
        """

        self.info = info
        self.parameters = parameters
        self.screening_charge = screening_charge
        self.energies = energies
        self.element_radius = element_radius

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of CosmoData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.
        Can be used for scf and escf/egrad.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            CosmoData.
        """
        info = parser.cosmo_header

        results = parser.cosmo_results
        parameters = screening_charge = energies = element_radius = None

        if results:
            parameters = results["parameters"]
            screening_charge = results["screening_charge"]
            energies = results["energies"]
            if results["element_radius"] is not None:
                element_radius = {}
                for el, el_data in results["element_radius"].items():
                    element_radius[el.capitalize()] = el_data

        if all(i is None for i in [info, parameters, screening_charge, energies, element_radius]):
            return None

        return cls(info=info, parameters=parameters, screening_charge=screening_charge,
                   energies=energies, element_radius=element_radius)


class TurbomoleData(BaseData):
    """
    Information about the Turbomole version and executable used.
    Can be used for all Turbomole executables.
    """

    def __init__(self, version=None, build=None, executable=None):
        """
        Args:
            version (str): number of the Turbomole version.
            build (str): Turbomole build.
            executable (str): name of the executable in the header of the output
        """
        self.version = version
        self.build = build
        self.executable = executable

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of TurbomoleData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            TurbomoleData.
        """

        header = parser.header

        if not header:
            return None

        return cls(version=header["tm_version"], build=header["tm_build"], executable=header["executable"])


class RunData(BaseData):
    """
    Information about where the calculation was executed and the timings.
    Can be used for all Turbomole executables.
    """

    def __init__(self, host=None, start_time=None, end_time=None, cpu_time=None, wall_time=None):
        """
        Args:
            host (str): name of the host.
            start_time (datetime): initial time of the run.
            end_time (datetime): end time of the run.
            cpu_time (float): cpu time in seconds.
            wall_time (float): wall time in seconds.
        """
        self.host = host
        self.start_time = start_time
        self.end_time = end_time
        self.cpu_time = cpu_time
        self.wall_time = wall_time


    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of RunData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            RunData.
        """
        kwargs = {}

        header = parser.header
        time_data = parser.timings
        if not header and not time_data:
            return None

        kwargs = {"start_time": None,
                  "host": None,
                  "end_time": None,
                  "cpu_time": None,
                  "wall_time": None}

        if header:
            kwargs["start_time"] = header["start_time"]
            kwargs["host"] = header["host"]

        time_data = parser.timings
        if time_data:
            kwargs["end_time"] = time_data["end_time"]
            kwargs["cpu_time"] = time_data["cpu_time"]
            kwargs["wall_time"] = time_data["wall_time"]

        return cls(**kwargs)


class BasisData(BaseData):
    """
    Information about the basis used for the calculation.
    Can be used for most of the Turbomole executables (including the scf, escf, grad).
    """
    def __init__(self, basis_per_specie=None, aux_basis_per_specie=None,
                 number_scf_basis_func=None, number_scf_aux_basis_func=None):
        """
        Args:
            basis_per_specie (dict): dict with species as keys and name of the basis as values.
            aux_basis_per_specie (dict): dict with species as keys and name of the auxiliary
                basis as values.
            number_scf_basis_func (int): number of scf basis functions.
            number_scf_aux_basis_func (int): number of auxialiary scf basis functions.
        """
        self.basis_per_specie = basis_per_specie
        self.aux_basis_per_specie = aux_basis_per_specie
        self.number_scf_basis_func = number_scf_basis_func
        self.number_scf_aux_basis_func = number_scf_aux_basis_func

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of BasisData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            BasisData.
        """
        basis_data = parser.basis

        if not basis_data:
            return None

        return cls(**basis_data)


class SymmetryData(BaseData):
    """
    Information on the symmetry of the molecule.
    Can be used for all the TM executables (in some case only a part of the
    information might be available).
    """
    def __init__(self, symbol=None, n_reps=None, reps=None):
        """
        Args:
            symbol (str): symbol of the symmetry of the molecule.
            n_reps (int): number of representations
            reps (list): the symbols of the representations.
        """
        self.symbol = symbol
        self.n_reps = n_reps
        self.reps = reps

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of SymmetryData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            SymmetryData.
        """
        sym_data = parser.symmetry

        if not sym_data:
            return None

        return cls(**sym_data)


class FunctionalData(BaseData):
    """
    Information about the exchange-correlation functional.
    Can be used for for scf, gradient and relax executables.
    """

    def __init__(self, msg=None, name=None, func_type=None, xcfun=None):
        """
        Args:
            msg (str): the full message in the output concerning the XC functional.
                In case the parsing failed to identify the correct name it will provide
                a way to figure out the information.
            name (str): name of the XC funtional.
            func_type (str): type of the XC functional (e.g. GGA, LDA, ...).
            xcfun (str): version of xcfun used, if present.
        """
        self.msg = msg
        self.name = name
        self.func_type = func_type
        self.xcfun = xcfun

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of FunctionalData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            FunctionalData.
        """
        func_data = parser.density_functional_data
        if not func_data:
            return None

        return cls(msg=func_data["functional_msg"], name=func_data["functional_name"],
                   func_type=func_data["functional_type"], xcfun=func_data["xcfun"])


class RiData(BaseData):
    """
    Information about RI calculations.
    Can be used for ridft and escf/egrad
    """

    def __init__(self, ricore=None, marij=None, rij_memory=None, rik=None):
        """
        Args:
            ricore (int): memory used in ricore.
            marij (bool): True if marij approximation used.
            rij_memory (int): memory for rij.
            rik (bool): True if rik approximation used.
        """
        self.ricore = ricore
        self.marij = marij
        self.rij_memory = rij_memory
        self.rik = rik

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of RiData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            RiData.
        """
        ri_data = parser.rij_info

        if not ri_data:
            return None

        return cls(**ri_data)


class DispersionCorrectionData(BaseData):
    """
    Information about the dispersion correction used in the calculation.
    Can be used for scf executables.
    """

    def __init__(self, correction=None, en_corr=None):
        """
        Args:
            correction (str): the name of the correction (e.g. D1, D2, ...).
            en_corr (float): correction on the total energy.
        """
        self.correction = correction
        self.en_corr = en_corr

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of DispersionCorrectionData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            DispersionCorrectionData.
        """
        d = parser.dftd
        if not d:
            return None
        return cls(**d)


class DFTData(BaseData):
    """
    The information about a dft calculation.
    Can be used for scf, gradients and escf executables.
    """

    def __init__(self, functional=None, ri=None, spherical_gridsize=None,
                 gridpoints=None, dispersion_correction=None):
        """

        Args:
            functional (FunctionalData): information about the exchange-correlation functional.
            ri (RiData): information RI calculations.
            spherical_gridsize (int): size of the grid for spherical integration.
            gridpoints (int): number of points for spherical integration.
            dispersion_correction (DispersionCorrectionData): information about
                dispersion correction.
        """
        self.functional = functional
        self.ri = ri
        self.spherical_gridsize = spherical_gridsize
        self.gridpoints = gridpoints
        self.dispersion_correction = dispersion_correction

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of DFTData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            DFTData.
        """
        data = parser.density_functional_data
        if not data:
            spherical_gridsize = None
            gridpoints = None
        else:
            spherical_gridsize = data["spherical_gridsize"]
            gridpoints = data["gridpoints"]

        func = FunctionalData.from_parser(parser)
        ri = RiData.from_parser(parser)
        disp = DispersionCorrectionData.from_parser(parser)

        if all(i is None for i in [spherical_gridsize, gridpoints, func, ri, disp]):
            return None

        return cls(functional=func, ri=ri, spherical_gridsize=spherical_gridsize,
                   gridpoints=gridpoints, dispersion_correction=disp)


class ScfIterationData(BaseData):
    """
    Details about the iteration in a scf calculation.
    It contains the value for each step as lists, but only keeps the initial
    and final indices to avoid wasting space. The index may start from
    values higher than 1 if it is a restart of a previous scf calculation.
    """

    def __init__(self, energies=None, first_index=None, n_steps=None, dampings=None,
                 converged=None):
        """
        Args:
            energies (list): values of the energies for each step of the scf loop.
            first_index (int): first index of the loop.
            n_steps (int): number of scf steps.
            dampings (list): values of the dampings for each step of the scf loop.
            converged (bool): True if the calculation converged.
        """
        self.energies = energies
        self.first_index = first_index
        self.n_steps = n_steps
        self.dampings = dampings
        self.converged = converged

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of ScfIterationData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            ScfIterationData.
        """
        iter_data = parser.scf_iterations
        if not iter_data:
            return None

        return cls(energies=iter_data["energies"], first_index=iter_data["first_index"],
                   n_steps=iter_data["n_steps"], dampings=iter_data["dampings"], converged=iter_data["converged"])

    @add_fig_kwargs
    def plot_energies(self, ax=None, **kwargs):
        """
        Plots the evolution of the energies in the scf loop.

        Args:
            ax (Axes):  a matplotlib Axes or None if a new figure should be created.
            kwargs: arguments passed to matplotlib plot function.

        Returns:
            A matplotlib Figure.
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        ax.plot(range(len(self.energies)), self.energies, **kwargs)
        ax.set_xlabel("Steps")
        ax.set_ylabel("Energy")

        return fig


class ScfData(BaseData):
    """
    Information about options and operations in an scf calculation.
    Can be used for scf calculations.
    """

    def __init__(self, iterations=None, diis=None, diis_error_vect=None, conv_tot_en=None,
                 conv_one_e_en=None, virtual_orbital_shift_on=None, virtual_orbital_shift_limit=None,
                 orbital_characterization=None, restart_file=None, n_occupied_orbitals=None):
        """
        Args:
            iterations (ScfIterationData):  details about the iteration in the scf loop.
            diis (bool): True if DIIS is switched on.
            diis_error_vect (str): type of DIIS error vector.
            conv_tot_en (float): criterion for scf convergence on increment of total energy.
            conv_one_e_en (float): criterion for scf convergence on increment of one-electron energy.
            virtual_orbital_shift_on (bool): True if automatic virtual orbital shift is switched on.
            virtual_orbital_shift_limit (float): automatic virtual orbital shift switched when
                e(lumo)-e(homo) lower than this value.
            orbital_characterization (str): type of orbital characterization.
            restart_file (str): file to which restart information will be dumped.
            n_occupied_orbitals (int): number of occupied orbitals.
        """
        self.iterations = iterations
        self.diis = diis
        self.diis_error_vect = diis_error_vect
        self.conv_tot_en = conv_tot_en
        self.conv_one_e_en = conv_one_e_en
        self.virtual_orbital_shift_on = virtual_orbital_shift_on
        self.virtual_orbital_shift_limit = virtual_orbital_shift_limit
        self.orbital_characterization = orbital_characterization
        self.restart_file = restart_file
        self.n_occupied_orbitals = n_occupied_orbitals

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of ScfData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            ScfData.
        """
        iterations = ScfIterationData.from_parser(parser)

        scf_data = parser.pre_scf_run

        if not iterations and not scf_data:
            return None

        kwargs = {}

        if iterations:
            kwargs["iterations"] = iterations

        if scf_data:
            kwargs.update(scf_data)

        return cls(**kwargs)


class ScfEnergiesData(BaseData):
    """
    Final energies and different contributions obtained from an scf calculation.
    Can be used for scf calculations.
    """

    def __init__(self, total_energy=None, kinetic_energy=None, potential_energy=None,
                 virial_theorem=None, wavefunction_norm=None, coulomb_energy=None,
                 xc_energy=None, ts_energy=None, free_energy=None, sigma0_energy=None):
        """
        Args:
            total_energy (float):
            kinetic_energy (float):
            potential_energy (float):
            virial_theorem (float):
            wavefunction_norm (float):
            coulomb_energy (float):
            xc_energy (float):
            ts_energy (float):
            free_energy (float):
            sigma0_energy (float):
        """
        self.total_energy = total_energy
        self.kinetic_energy = kinetic_energy
        self.potential_energy = potential_energy
        self.virial_theorem = virial_theorem
        self.wavefunction_norm = wavefunction_norm
        self.coulomb_energy = coulomb_energy
        self.xc_energy = xc_energy
        self.ts_energy = ts_energy
        self.free_energy = free_energy
        self.sigma0_energy = sigma0_energy

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of ScfEnergiesData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            ScfEnergiesData.
        """

        en_data = parser.scf_energies
        riper_en_data = parser.riper_scf_energies

        if en_data:
            if riper_en_data:
                raise RuntimeError('Found scf energy data from dscf/ridft as well as from riper.')
            return cls(**en_data)
        elif riper_en_data:
            return cls(**riper_en_data)

        return None


class ElectrostaticMomentsData(BaseData):
    """
    The data of the electrostatic moments (charge, dipole and quadrupole).
    Can be used for scf executables.
    """

    def __init__(self, charge=None, unrestricted_electrons=None, dipole_vector=None,
                 dipole_norm=None, quadrupole_tensor=None, quadrupole_trace=None,
                 quadrupole_anisotropy=None):
        """
        Args:
            charge (float): total charge.
            unrestricted_electrons (float): number of unrestricted electrons.
            dipole_vector (list): the 3 components of the dipole moment.
            dipole_norm (float): norm of the dipole moment.
            quadrupole_tensor (list): 3x3 matrix with the quadrupole tensor.
            quadrupole_trace (float): trace of the quadrupole tensor.
            quadrupole_anisotropy (float): anisotropy of the quadrupole tensor.
        """
        self.charge = charge
        self.unrestricted_electrons = unrestricted_electrons
        self.dipole_vector = dipole_vector
        self.dipole_norm = dipole_norm
        self.quadrupole_tensor = quadrupole_tensor
        self.quadrupole_trace = quadrupole_trace
        self.quadrupole_anisotropy = quadrupole_anisotropy

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of ElectrostaticMomentsData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            ElectrostaticMomentsData.
        """
        mom_data = parser.electrostatic_moments

        if not mom_data:
            return None

        dip = mom_data["dipole"]
        quad = mom_data["quadrupole"]
        return cls(charge=mom_data["charge"], unrestricted_electrons=mom_data["unrestricted_electrons"],
                   dipole_vector=dip["moment"], dipole_norm=dip["norm"], quadrupole_tensor=quad["moment"],
                   quadrupole_trace=quad["trace"], quadrupole_anisotropy=quad["anisotropy"])


class GeometryData(BaseData):
    """
    Data of the geometry of the system: molecule and centers.
    Can be used for most of the Turbomole executables (including the scf, escf, grad).
    """

    def __init__(self, center_of_mass=None, center_of_charge=None, molecule=None):
        """
        Args:
            center_of_mass (list): 3D vector with the position of the center of mass.
            center_of_charge: 3D vector with the position of the center of charge.
            molecule (Molecule): the molecule in the system.
        """
        self.center_of_mass = center_of_mass
        self.center_of_charge = center_of_charge
        self.molecule = molecule

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of GeometryData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            GeometryData.
        """
        centers_data = parser.centers

        coords_data = parser.coordinates

        if not centers_data and not coords_data:
            return None

        kwargs = dict(center_of_mass=None, center_of_charge=None, molecule=None)

        if centers_data:
            kwargs.update(centers_data)

        if coords_data:
            species = [s.capitalize() for s in coords_data["species"]]
            coords = np.array(coords_data["coords"]) * bohr_to_ang
            mol = Molecule(species, coords)
            kwargs["molecule"] = mol

        return cls(**kwargs)


class SpinData(BaseData):
    """
    Information about the spin in the calculation.
    Can be used for scf, gradient and escf executables.
    """

    def __init__(self, unrestricted=None, s2=None):
        """
        Args:
            unrestricted (bool): True if is an uhf calculation.
            s2 (float): S^2 value of the spin.
        """
        self.unrestricted = unrestricted
        self.s2 = s2

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of SpinData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            SpinData.
        """
        return cls(unrestricted=parser.is_uhf, s2=parser.s2)


class SmearingData(BaseData):
    """
    Information about the smearing ($fermi datagroup).
    Can be used for scf executables.
    """

    def __init__(self, initial_elec_temp=None, final_elec_temp=None, annealing_factor=None,
                 annealing_homo_lumo_gap_limit=None, smearing_de_limit=None):
        """
        Args:
            initial_elec_temp (float): initial electron temperature.
            final_elec_temp (float): final electron temperature.
            annealing_factor (float): annealing factor
            annealing_homo_lumo_gap_limit (float): annealing if HOMO-LUMO gap lower than this value.
            smearing_de_limit (float): smearing switched off if DE lower than this value.
        """
        self.initial_elec_temp = initial_elec_temp
        self.final_elec_temp = final_elec_temp
        self.annealing_factor = annealing_factor
        self.annealing_homo_lumo_gap_limit = annealing_homo_lumo_gap_limit
        self.smearing_de_limit = smearing_de_limit

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of SmearingData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            SmearingData.
        """
        fermi_data = parser.fermi

        if not fermi_data:
            return None

        return cls(**fermi_data)


class IntegralData(BaseData):
    """
    Data about the thresholds for integrals.
    Can be used for for scf, gradien and escf executables.
    """

    def __init__(self, integral_neglect_threshold=None, thize=None, thime=None):
        """
        Args:
            integral_neglect_threshold (float): integral neglect threshold.
            thize (float): integral storage threshold THIZE.
            thime (float): integral storage threshold THIME.
        """
        self.integral_neglect_threshold = integral_neglect_threshold
        self.thize = thize
        self.thime = thime

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of IntegralData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            IntegralData.
        """
        integral_data = parser.integral

        if not integral_data:
            return None

        return cls(**integral_data)


class EscfIterationData(BaseData):
    """
    Details about the iteration in an escf calculation.
    It contains the data for a list of steps. For each one a sublist of the convergence
    information converning each of the irreps treated as excited states.
    """

    def __init__(self, steps=None, converged=None):
        """
        Args:
            steps (list): list of lists of lists (3 dimensions). One element of the list for
                each escf step. One element of the sublist for each irrep dealt with. The inner list
                contains the convergence information:
                [name of the irrep, number of converged roots, euclidean residual norm]
            converged (bool): True if converged.
        """
        self.steps = steps
        self.converged = converged

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of EscfIterationData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            EscfIterationData.
        """
        iter_data = parser.escf_iterations
        if not iter_data:
            return None

        return cls(**iter_data)


class SingleExcitation(MSONable):
    """
    Auxiliary object to store the information about an escf/egrad calculation.
    It cannot be extracted directly from the string, being present several
    times in the output of an escf/egrad calculation.
    """

    def __init__(self, tot_en=None, osc_stre=None, rot_stre=None, dominant_contributions=None,
                 moments_columns=None):
        """
        Args:
            tot_en (float): total energy.
            osc_stre (float): oscillator strength, length representation.
            rot_stre (float): rotatory strength, length representation.
            dominant_contributions (list): list of all the dominant contributions.
                Each element is a dictionary with "occ_orb" (a dictionary with index of
                the occupied orbital, irrep, energy and spin), "virt_orb" (same as "occ_orb")
                and "coeff" (the coefficient of the contribution).
            moments_columns (list): a list of dictionaries containing the electric an magnetic
                moments for each column.
                Each dictionary should contain the following keys:
                    electric_dipole (list): 3D vector with the electronic dipole.
                    magnetic_dipole (list): 3D vector with the magnetic dipole.
                    electric_quadrupole (dict): description of the electronic quadrupole with keys
                        moment (the 3x3 matrix of the quadrupole moment), "trace" and "anisotropy".
        """
        self.tot_en = tot_en
        self.osc_stre = osc_stre
        self.rot_stre = rot_stre
        self.dominant_contributions = dominant_contributions
        self.moments_columns = moments_columns


class EscfData(BaseData):
    """
    Output of an escf calculation.
    Can be used for escf and egrad.
    """

    def __init__(self, calc_type=None, iterations=None, residuum_convergence_criterium=None,
                 n_occupied_orbitals=None, orbital_characterization=None, max_davidson_iter=None,
                 machine_precision=None, max_core_mem=None, max_cao_basis_vectors=None,
                 max_treated_vectors=None, irrep_data=None, gs_tot_en=None, excitations=None):
        """
        Args:
            calc_type (str): string describing the calculation type (e.g. 'RPA SINGLET-EXCITATION').
            iterations (EscfIterationData): the data about the escf iterations.
            residuum_convergence_criterium (float): residuum convergence criterium.
            n_occupied_orbitals (int): number of occupied orbitals.
            orbital_characterization (str): description of how the orbital were converged in scf.
            max_davidson_iter (int): maximum number of Davidson iterations.
            machine_precision (float): machine precision.
            max_core_mem (int): maximum core memory in MB.
            max_cao_basis_vectors (int): number of maximum CAO vectors in memory.
            max_treated_vectors (int): maximum number of simultaneously treated vectors,
                including degeneracy.
            irrep_data (dict): keys are strings with the name of the irrep, values are tuples
                with [tensor space dimension, number of roots].
            gs_tot_en (float): total energy of the ground state. Turbomole extracts it from the
                output of the scf calculation. If that is missing Turbomole sets this value to 0
                in the output of escf.
            excitations (dict): keys are the name of the irreps and values are lists of SingleExcitation.
        """
        self.calc_type = calc_type
        self.iterations = iterations
        self.residuum_convergence_criterium = residuum_convergence_criterium
        self.n_occupied_orbitals = n_occupied_orbitals
        self.orbital_characterization = orbital_characterization
        self.max_davidson_iter = max_davidson_iter
        self.machine_precision = machine_precision
        self.max_core_mem = max_core_mem
        self.max_cao_basis_vectors = max_cao_basis_vectors
        self.max_treated_vectors = max_treated_vectors
        self.irrep_data = irrep_data
        self.gs_tot_en = gs_tot_en
        self.excitations = excitations

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of EscfIterationData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            EscfIterationData.
        """
        iterations = EscfIterationData.from_parser(parser)

        excitations = parser.escf_excitations

        escf_data = parser.pre_escf_run

        if not iterations and not escf_data and not excitations:
            return None

        kwargs = dict(
            calc_type=None,
            residuum_convergence_criterium=None,
            n_occupied_orbitals=None,
            orbital_characterization=None,
            max_davidson_iter=None,
            machine_precision=None,
            max_core_mem=None,
            max_cao_basis_vectors=None,
            max_treated_vectors=None,
            irrep_data=None
        )

        kwargs["iterations"] = iterations
        converted_excitations = {}
        for irrep, l in excitations.items():
            converted_excitations[irrep] = [SingleExcitation(**d) for d in l]
        kwargs["excitations"] = converted_excitations
        kwargs["gs_tot_en"] = parser.escf_gs_total_en

        if escf_data:
            kwargs.update(escf_data)

        return cls(**kwargs)


class StatptData(BaseData):
    """
    Initial information provided in statpt.
    Can be used only for statpt.
    """

    def __init__(self, max_trust_radius=None, min_trust_radius=None, init_trust_radius=None,
                 min_grad_norm_for_gdiis=None, prev_steps_for_gdiis=None, hessian_update_method=None,
                 thr_energy_change=None, thr_max_displ=None, thr_max_grad=None, thr_rms_displ=None,
                 thr_rms_grad=None, use_defaults=None, final_step_radius=None):
        """
        Args:
            max_trust_radius (float): maximum allowed trust radius.
            min_trust_radius (float): minimum allowed trust radius
            init_trust_radius (float): initial trust radius.
            min_grad_norm_for_gdiis (float): GDIIS used if gradient norm is lower
                than this value.
            prev_steps_for_gdiis (int): number of previous steps for GDIIS.
            hessian_update_method (str): hessian update method.
            thr_energy_change (float): threshold for energy change.
            thr_max_displ (float): threshold for max displacement element.
            thr_max_grad (float): threshold for max gradient element.
            thr_rms_displ (float): threshold for RMS of displacement.
            thr_rms_grad (float): threshold for RMS of gradient.
            use_defaults (bool): True if $statpt was not defined and using default options.
            final_step_radius (float): final step radius.
        """
        self.max_trust_radius = max_trust_radius
        self.min_trust_radius = min_trust_radius
        self.init_trust_radius = init_trust_radius
        self.min_grad_norm_for_gdiis = min_grad_norm_for_gdiis
        self.prev_steps_for_gdiis = prev_steps_for_gdiis
        self.hessian_update_method = hessian_update_method
        self.thr_energy_change = thr_energy_change
        self.thr_max_displ = thr_max_displ
        self.thr_max_grad = thr_max_grad
        self.thr_rms_displ = thr_rms_displ
        self.thr_rms_grad = thr_rms_grad
        self.use_defaults = use_defaults
        self.final_step_radius = final_step_radius

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of StatptData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            StatptData.
        """
        statpt_data = parser.statpt_info
        if not statpt_data:
            return None

        return cls(**statpt_data)


class RelaxData(BaseData):
    """
    Initial information provided in relax.
    Can be used only for relax.
    """

    def __init__(self, optimizations=None, thr_int_coord=None):
        """
        Args:
            optimizations (list): list of strings describing with respect to what the
                optimization has been performed.
            thr_int_coord (float): convergence criterion for internal coordinates.
        """
        self.optimizations = optimizations
        self.thr_int_coord = thr_int_coord

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of RelaxData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            RelaxData.
        """
        relax_data = parser.relax_info
        if not relax_data:
            return None

        return cls(**relax_data)


class RelaxGradientsData(BaseData):
    """
    Gradient values extracted from the relax/stapt output that take into account
    the frozen coordinates.
    Can be used for relax and statpt.
    """

    def __init__(self, norm_cartesian=None, norm_internal=None, max_internal=None):
        """
        Args:
            norm_cartesian (float): norm of the cartesian gradient.
            norm_internal (float): norm of the internal gradient.
            max_internal (float): maximum norm of the internal gradient (available
                only in relax).
        """
        self.norm_cartesian = norm_cartesian
        self.norm_internal = norm_internal
        self.max_internal = max_internal

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of RelaxGradientsData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            RelaxGradientsData.
        """
        grad_data = parser.relax_gradient_values
        if not grad_data:
            return None

        return cls(**grad_data)


class RelaxConvergenceData(BaseData):
    """
    Final information about convergence according to what is defined in the
    control file.
    Can be used for relax and statpt.
    """

    def __init__(self, energy_change=None, rms_displ=None, max_displ=None, rms_grad=None, max_grad=None):
        """
        Args:
            energy_change (dict): information about the convergence of the energy change.
                Contains "value" with the actual value, "thr" with the specified threshold,
                and "conv" a bool specifying wether the convergence for this threshold has been
                achieved or not.
            rms_displ (dict): information about the convergence of the rms of the displacements.
                Same structure as "energy_change".
            max_displ (dict): information about the convergence of the maximum of the displacements.
                Same structure as "energy_change".
            rms_grad (dict): information about the convergence of the rms of the gradient.
                Same structure as "energy_change".
            max_grad (dict): information about the convergence of the maximum of the gradient.
                Same structure as "energy_change".
        """
        self.energy_change = energy_change
        self.rms_displ = rms_displ
        self.max_displ = max_displ
        self.rms_grad = rms_grad
        self.max_grad = max_grad

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of RelaxConvergenceData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            RelaxConvergenceData.
        """
        conv_data = parser.relax_conv_info
        if not conv_data:
            return None

        return cls(**conv_data)


class AoforceNumericalIntegrationData(BaseData):
    """
    Information about the numerical integration in aoforce.
    Will be present only if a proper aoforce is run. Absent if running after numforce.
    """

    def __init__(self, core_memory_dft=None, memory_per_atom=None, atoms_per_loop=None,
                 construction_timings=None):
        """
        Args:
            core_memory_dft (int): remaining core memory for DFT in MB.
            memory_per_atom (int): memory needed per atom in KB.
            atoms_per_loop (int): atoms per loop corresponding to the memory.
            construction_timings (list): a list of lists with construction timings.
                For each element:
                [str with description of the timing, cpu time in seconds, wall time in seconds]
        """
        self.core_memory_dft = core_memory_dft
        self.memory_per_atom = memory_per_atom
        self.atoms_per_loop = atoms_per_loop
        self.construction_timings = construction_timings

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of AoforceNumericalIntegrationData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            AoforceNumericalIntegrationData.
        """
        data = parser.aoforce_numerical_integration
        if not data:
            return None

        kwargs = {"construction_timings": data["construction_timings"]}
        kwargs.update(data["memory"])

        return cls(**kwargs)


class AoforceRotationalData(BaseData):
    """
    Analysis of rotational states in aoforce.
    """

    def __init__(self, b=None, intensities=None, m=None, dipole_moment=None):
        """
        Args:
            b (float): 3D vector with rotational constants b in cm^-1.
            intensities (list): 3D vector with optical intensities in a.u.
            m (list): 3x3 matrix.
            dipole_moment (list): 3D vector with dipole moment in principle axis system in a.u.
        """
        self.b = b
        self.intensities = intensities
        self.m = m
        self.dipole_moment = dipole_moment

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of AoforceRotationalData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            AoforceRotationalData.
        """
        data = parser.aoforce_analysis
        if not data or "rotational" not in data:
            return None

        return cls(**data["rotational"])


class AoforceVibrationalData(BaseData):
    """
    Analysis of vibrational states in aoforce.
    Results are stored as lists, each value corresponds to the frequency
    in the "frequencies" list with the same index.
    """

    def __init__(self, frequencies=None, symmetries=None, ir=None, dDIP_dQ=None,
                 intensities=None, intensities_perc=None, raman=None, eigenvectors=None,
                 reduced_masses=None, energies=None):
        """
        Args:
            frequencies (list): vibrational frequencies in cm^-1.
            symmetries (list): string with the name of the symmetry of the modes.
                None if not specified, like for the first 6 modes.
            ir (list): booleans stating whether if the mode is ir or not.
                None if "-" in the output file.
            dDIP_dQ (list): normal mode derivatives of the dipole moment in a.u.
            intensities (list): intensity in km/mol.
            intensities_perc (list): intensity percentage.
            raman (list):booleans stating whether if the mode is raman or not.
                None if "-" in the output file.
            eigenvectors (list): matrix with shape (num frequencies, number of atoms, 3)
                with the eigenvectors for each atom.
            reduced_masses (list): reduced masses in g/mol.
            energies (dict): the values of the vibrational energies. A dictionary with "zpve"
                "scf" and "total" as keys.
        """
        self.frequencies = frequencies
        self.symmetries = symmetries
        self.ir = ir
        self.dDIP_dQ = dDIP_dQ
        self.intensities = intensities
        self.intensities_perc = intensities_perc
        self.raman = raman
        self.eigenvectors = eigenvectors
        self.reduced_masses = reduced_masses
        self.energies = energies

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of AoforceVibrationalData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            AoforceVibrationalData.
        """
        data = parser.aoforce_analysis
        if not data or "vibrational" not in data:
            return None

        return cls(**data["vibrational"])

    def get_freqs_df(self):
        """
        Generates a pandas DataFrame with the frequencies and their properties

        Returns:
            DataFrame
        """
        cols = np.column_stack([self.frequencies, self.symmetries, self.intensities, self.intensities_perc, self.ir,
                                self.raman, self.dDIP_dQ, self.reduced_masses])
        labels = ["frequency", "symmetry", "intensity", "intensity %", "IR", "Raman", "|dDIP/dQ|", "red. mass"]
        return pd.DataFrame(cols, columns=labels)

    def n_negative_freqs(self, tol=0.1):
        """
        The number of frequencies that satisfy the condition
        f < -abs(tol)

        Args:
            tol (float): tolerance for considering a frequency as negative.

        Returns:
            int: the number of negative frequencies.
        """

        return np.count_nonzero(np.array(self.frequencies) < -np.abs(tol))

    def n_positive_freqs(self, tol=0.1):
        """
        The number of frequencies that satisfy the condition
        f > abs(tol)

        Args:
            tol (float): tolerance for considering a frequency as positive.

        Returns:
            int: the number of positive frequencies.
        """

        return np.count_nonzero(np.array(self.frequencies) > np.abs(tol))

    def n_zero_freqs(self, tol=0.1):
        """
        The number of frequencies that satisfy the condition
        -abs(tol) > f > abs(tol)

        Args:
            tol (float): tolerance for considering a frequency as zero.

        Returns:
            int: the number of zero frequencies.
        """
        f_arr = np.array(self.frequencies)
        return np.count_nonzero(np.logical_and(-np.abs(tol) <= f_arr, f_arr <= np.abs(tol)))


class MP2Data(BaseData):
    """
    MP2 data object containing parameters used to run the MP2 calculation.
    """

    def __init__(self, energy_only=None):
        """
        Args:
            energy_only (bool): whether this is an energy-only MP2 calculation.
        """
        self.energy_only = energy_only

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of MP2Data from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            MP2Data.
        """
        data = parser.mp2_data
        if not data:
            return None

        return cls(**data)


class MP2Results(BaseData):
    """
    MP2 results object.
    """

    def __init__(self, energy=None):
        """
        Args:
            energy (float): MP2 energy.
        """
        self.energy = energy

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of MP2Results from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            MP2Results.
        """
        data = parser.mp2_results
        if not data or "energy" not in data:
            return None

        return cls(**data)


class PeriodicityData(BaseData):
    """
    Information about periodicity
    """

    def __init__(self, periodicity=None, lattice_params=None,
                 shortest_interatomic_distance=None,
                 direct_space_vectors=None, reciprocal_space_vectors=None):
        """
        Args:
            periodicity (int): Periodicity of the system (1, 2 or 3).
            lattice_params: Lattice parameters of the system (in Angstroms).
                For 1D systems, a single number: [a]
                For 2D systems, three numbers: [a, b, gamma], i.e. the two lattice
                    parameters and the angle between them.
                For 3D systems, three numbers: [a, b, c, alpha, beta, gamma], i.e.
                    three two lattice parameters and the angles between them.
            shortest_interatomic_distance: Shortest interatomic distance in the
                system (in Angstroms).
            direct_space_vectors: Lattice vectors of the system (in Angstroms).
                The first vector will always be aligned with the x cartesian direction.
                For 2D and 3D, the second vector is always in the xy cartesian plane.
            reciprocal_space_vectors: Reciprocal lattice vectors of the system (in Angstroms^-1).
                The physics definition of the reciprocal lattice vectors is used here, i.e.
                    b1 = 2*pi/V * (a2 x a3)
                    b2 = 2*pi/V * (a3 x a1)
                    b3 = 2*pi/V * (a1 x a2)
                The crystallographer's definition is easily recovered as:
                    b'1 = b1 / (2*pi)
                    b'2 = b2 / (2*pi)
                    b'3 = b3 / (2*pi)

        Returns:
            PeriodicityData.

        """
        self.periodicity = periodicity
        self.lattice_params = lattice_params
        self.shortest_interatomic_distance = shortest_interatomic_distance
        self.direct_space_vectors = direct_space_vectors
        self.reciprocal_space_vectors = reciprocal_space_vectors

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of PeriodicityData from a parser based on the stdout
        of a Turbomole executable. Returns None if no data could be parsed.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            PeriodicityData.
        """
        data = parser.periodicity_data
        if not data:
            return None

        periodicity = data["periodicity"]
        if periodicity == 1:
            # a
            lattice_params = [bohr_to_ang * data["tm_lattice_params"][0]]
        elif periodicity == 2:
            # a, b, gamma
            lattice_params = [bohr_to_ang * data["tm_lattice_params"][0],
                              bohr_to_ang * data["tm_lattice_params"][1],
                              data["tm_lattice_params"][2]]
        elif periodicity == 3:
            # a, b, c, alpha, beta, gamma
            lattice_params = [bohr_to_ang * data["tm_lattice_params"][0],
                              bohr_to_ang * data["tm_lattice_params"][1],
                              bohr_to_ang * data["tm_lattice_params"][2],
                              data["tm_lattice_params"][3],
                              data["tm_lattice_params"][4],
                              data["tm_lattice_params"][5]]
        else:
            raise RuntimeError('"periodicity" should be 1, 2 or 3.')

        direct_space_vectors = [[bohr_to_ang * xx for xx in vect]
                                for vect in data["direct_space_vectors"]]
        reciprocal_space_vectors = [[xx / bohr_to_ang for xx in vect]
                                    for vect in data["reciprocal_space_vectors"]]

        return cls(periodicity=periodicity,
                   lattice_params=lattice_params,
                   shortest_interatomic_distance=bohr_to_ang*data["shortest_interatomic_distance"],
                   direct_space_vectors=direct_space_vectors,
                   reciprocal_space_vectors=reciprocal_space_vectors
                   )
