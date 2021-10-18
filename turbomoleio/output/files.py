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

from turbomoleio.output.data import BaseData, CosmoData, TurbomoleData, RunData, BasisData, IntegralData
from turbomoleio.output.data import DFTData, ScfData, ScfEnergiesData, ElectrostaticMomentsData, GeometryData, SpinData
from turbomoleio.output.data import SmearingData, EscfData, AoforceNumericalIntegrationData, AoforceRotationalData
from turbomoleio.output.data import AoforceVibrationalData, RelaxConvergenceData, RelaxData, RelaxGradientsData
from turbomoleio.output.data import StatptData, SymmetryData
from turbomoleio.output.data import MP2Data, MP2Results


class ScfOutput(BaseData):
    """
    Object containing the data of the output of an scf calculation, i.e. dscf or ridft.
    """

    def __init__(self, dft, scf, energies, electrostatic, geometry, basis, run, tm,
                 cosmo, spin, integral, smearing, symmetry):
        """
        Args:
            dft (DFTData): information about dft calculation.
            scf (ScfData): information about scf loop.
            energies (ScfEnergiesData): final energy of scf calculation.
            electrostatic (ElectrostaticMomentsData): information about electrostatic moments.
            geometry (GeometryData): the geometry of the system.
            basis (BasisData): the basis used for the calculation.
            run (RunData): information about calculation running (e.g. timings, ...)
            tm (TurbomoleData): information about the turbomole used for the calculation.
            cosmo (CosmoData): cosmo approximations and results.
            spin (SpinData): information about the spin of the system.
            integral (IntegralData): information about the thresholds for integrals.
            smearing (SmearingData): information about the smearing ($fermi datagroup).
            symmetry (SymmetryData): information about the symmetry of the molecule
        """
        self.dft = dft
        self.scf = scf
        self.energies = energies
        self.electrostatic = electrostatic
        self.geometry = geometry
        self.basis = basis
        self.run = run
        self.tm = tm
        self.cosmo = cosmo
        self.spin = spin
        self.integral = integral
        self.smearing = smearing
        self.symmetry = symmetry

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of ScfOutput from a parser based on the stdout
        of a Turbomole executable.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            ScfOutput.
        """
        dft = DFTData.from_parser(parser)
        scf = ScfData.from_parser(parser)
        energies = ScfEnergiesData.from_parser(parser)
        electrostatic = ElectrostaticMomentsData.from_parser(parser)
        geometry = GeometryData.from_parser(parser)
        basis = BasisData.from_parser(parser)
        run = RunData.from_parser(parser)
        tm = TurbomoleData.from_parser(parser)
        cosmo = CosmoData.from_parser(parser)
        spin = SpinData.from_parser(parser)
        integral = IntegralData.from_parser(parser)
        smearing = SmearingData.from_parser(parser)
        symmetry = SymmetryData.from_parser(parser)

        return cls(dft=dft, scf=scf, energies=energies, electrostatic=electrostatic, geometry=geometry,
                   basis=basis, run=run, tm=tm, cosmo=cosmo, spin=spin, integral=integral,
                   smearing=smearing, symmetry=symmetry)


class EscfOutput(BaseData):
    """
    Object containing the data of the output of an escf calculation.
    Parses also the information relative to the dft part, which is common to the
    output of the scf calculation that preceeds the escf.
    """

    def __init__(self, dft, escf, geometry, basis, run, tm, cosmo, integral, symmetry):
        """
        Args:
            dft (DFTData): information about dft calculation.
            escf (EscfData): information about the excited state calculation.
            geometry (GeometryData): the geometry of the system.
            basis (BasisData): the basis used for the calculation.
            run (RunData): information about calculation running (e.g. timings, ...)
            tm (TurbomoleData): information about the turbomole used for the calculation.
            cosmo (CosmoData): cosmo approximations and results.
            integral (IntegralData): information about the thresholds for integrals.
            symmetry (SymmetryData): information about the symmetry of the molecule
        """
        self.dft = dft
        self.escf = escf
        self.geometry = geometry
        self.basis = basis
        self.run = run
        self.tm = tm
        self.cosmo = cosmo
        self.integral = integral
        self.symmetry = symmetry

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of EscfOutput from a parser based on the stdout
        of a Turbomole executable.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            EscfOutput.
        """
        dft = DFTData.from_parser(parser)
        escf = EscfData.from_parser(parser)
        run = RunData.from_parser(parser)
        tm = TurbomoleData.from_parser(parser)
        geometry = GeometryData.from_parser(parser)
        basis = BasisData.from_parser(parser)
        cosmo = CosmoData.from_parser(parser)
        integral = IntegralData.from_parser(parser)
        symmetry = SymmetryData.from_parser(parser)

        return cls(dft=dft, escf=escf, geometry=geometry, basis=basis, run=run, tm=tm,
                   cosmo=cosmo, integral=integral, symmetry=symmetry)


class EscfOnlyOutput(BaseData):
    """
    Minimal version of the object containing the data of the output of an escf calculation.
    Parses only the information relative to escf calculation, ignoring the part that is common
    to the output of the scf calculation that preceeds the escf.
    """
    def __init__(self, escf, run, tm):
        """
        Args:
            escf (EscfData): information about the excited state calculation.
            run (RunData): information about calculation running (e.g. timings, ...)
            tm (TurbomoleData): information about the turbomole used for the calculation.
        """
        self.escf = escf
        self.run = run
        self.tm = tm

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of EscfOnlyOutput from a parser based on the stdout
        of a Turbomole executable.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            EscfOnlyOutput.
        """
        escf = EscfData.from_parser(parser)
        run = RunData.from_parser(parser)
        tm = TurbomoleData.from_parser(parser)

        return cls(escf=escf, run=run, tm=tm)


class GradOutput(BaseData):
    """
    Object containing the data of the output of a simple gradient calculation,
    i.e. grad and rdgrad.
    """
    def __init__(self, gradients, dielectric, run, tm, memory):
        """
        Args:
            gradients (list): matrix of shape (natoms, 3) containing the values of
                the cartesian gradients.
            dielectric (list): 3x3 matrix containing the values of the dielectric tensor.
            run (RunData): information about calculation running (e.g. timings, ...)
            tm (TurbomoleData): information about the turbomole used for the calculation.
            memory (int): memory in MB used. Available only for rdgrad.
        """
        self.gradients = gradients
        self.dielectric = dielectric
        self.run = run
        self.tm = tm
        self.memory = memory

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of GradOutput from a parser based on the stdout
        of a Turbomole executable.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            GradOutput.
        """
        grad_data = parser.gradient
        if not grad_data:
            return None

        run = RunData.from_parser(parser)
        tm = TurbomoleData.from_parser(parser)
        memory = parser.rdgrad_memory

        return cls(gradients=grad_data["gradients"], dielectric=grad_data["dielectric"],
                   run=run, tm=tm, memory=memory)


class EgradOutput(BaseData):
    """
    Object containing the data of the output of an egrad calculation.
    It combines the outputs of standard gradients calculations with the
    output of escf.
    """
    def __init__(self, gradients, dielectric, dft, escf, geometry, basis, run, tm, cosmo,
                 integral, ex_state):
        """
        Args:
            gradients (list): matrix of shape (natoms, 3) containing the values of
                the cartesian gradients.
            dielectric (list): 3x3 matrix containing the values of the dielectric tensor.
            dft (DFTData): information about dft calculation.
            escf (EscfData): information about the excited state calculation.
            geometry (GeometryData): the geometry of the system.
            basis (BasisData): the basis used for the calculation.
            run (RunData): information about calculation running (e.g. timings, ...)
            tm (TurbomoleData): information about the turbomole used for the calculation.
            cosmo (CosmoData): cosmo approximations and results.
            integral (IntegralData): information about the thresholds for integrals.
            ex_state (int): the (1-based) index of the excited state used for the gradients.
        """
        self.gradients = gradients
        self.dielectric = dielectric
        self.dft = dft
        self.escf = escf
        self.geometry = geometry
        self.basis = basis
        self.run = run
        self.tm = tm
        self.cosmo = cosmo
        self.integral = integral
        self.ex_state = ex_state

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of EgradOutput from a parser based on the stdout
        of a Turbomole executable.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            EgradOutput.
        """
        grad_data = parser.gradient
        if not grad_data:
            return None

        dft = DFTData.from_parser(parser)
        escf = EscfData.from_parser(parser)
        run = RunData.from_parser(parser)
        tm = TurbomoleData.from_parser(parser)
        geometry = GeometryData.from_parser(parser)
        basis = BasisData.from_parser(parser)
        cosmo = CosmoData.from_parser(parser)
        integral = IntegralData.from_parser(parser)
        egrad_ex = parser.egrad_excited_state
        ex_state = egrad_ex["index"] if egrad_ex is not None else None

        return cls(gradients=grad_data["gradients"], dielectric=grad_data["dielectric"],
                   dft=dft, escf=escf, geometry=geometry, basis=basis, run=run, tm=tm,
                   cosmo=cosmo, integral=integral, ex_state=ex_state)


class RelaxOutput(BaseData):
    """
    Object containing the data of the output of a relax calculation (not statpt).
    """
    def __init__(self, info, gradients, convergence, run, tm):
        """
        Args:
            info (RelaxData): initial information provided in relax.
            gradients (RelaxGradientsData): gradient values extracted from the relax
                output that take into account the frozen coordinates.
            convergence (RelaxConvergenceData): final information about convergence
                according to what is defined in the control file.
            run (RunData): information about calculation running (e.g. timings, ...)
            tm (TurbomoleData): information about the turbomole used for the calculation.
        """
        self.info = info
        self.gradients = gradients
        self.convergence = convergence
        self.run = run
        self.tm = tm

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of RelaxOutput from a parser based on the stdout
        of a Turbomole executable.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            RelaxOutput.
        """
        info = RelaxData.from_parser(parser)
        gradients = RelaxGradientsData.from_parser(parser)
        convergence = RelaxConvergenceData.from_parser(parser)
        run = RunData.from_parser(parser)
        tm = TurbomoleData.from_parser(parser)

        return cls(info=info, gradients=gradients, convergence=convergence, run=run, tm=tm)


class StatptOutput(BaseData):
    """
    Object containing the data of the output of a statpt calculation (not relax).
    """
    def __init__(self, info, gradients, convergence, run, tm):
        """
        Args:
            info (StatptData): initial information provided in statpt.
            gradients (RelaxGradientsData): gradient values extracted from the statpt
                output that take into account the frozen coordinates.
            convergence (RelaxConvergenceData): final information about convergence
                according to what is defined in the control file.
            run (RunData): information about calculation running (e.g. timings, ...)
            tm (TurbomoleData): information about the turbomole used for the calculation.
        """
        self.info = info
        self.gradients = gradients
        self.convergence = convergence
        self.run = run
        self.tm = tm

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of StatptOutput from a parser based on the stdout
        of a Turbomole executable.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            StatptOutput.
        """
        info = StatptData.from_parser(parser)
        gradients = RelaxGradientsData.from_parser(parser)
        convergence = RelaxConvergenceData.from_parser(parser)
        run = RunData.from_parser(parser)
        tm = TurbomoleData.from_parser(parser)

        return cls(info=info, gradients=gradients, convergence=convergence, run=run, tm=tm)


class AoforceOutput(BaseData):
    """
    Object containing the data of the output of an aforce calculation.
    Works native aoforce calculations and after a numforce.
    """
    def __init__(self, numerical_integration, rotational, vibrational, run, tm):
        """
        Args:
            numerical_integration (AoforceNumericalIntegrationData): information about the
                numerical integration in aoforce. Missing if from a numforce calculation.
            rotational (AoforceRotationalData): analysis of rotational states in aoforce.
            vibrational (AoforceVibrationalData): analysis of vibrational states in aoforce.
            run (RunData): information about calculation running (e.g. timings, ...)
            tm (TurbomoleData): information about the turbomole used for the calculation.
        """
        self.numerical_integration = numerical_integration
        self.rotational = rotational
        self.vibrational = vibrational
        self.run = run
        self.tm = tm

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of AoforceOutput from a parser based on the stdout
        of a Turbomole executable.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            AoforceOutput.
        """
        ni = AoforceNumericalIntegrationData.from_parser(parser)
        rot = AoforceRotationalData.from_parser(parser)
        vib = AoforceVibrationalData.from_parser(parser)
        run = RunData.from_parser(parser)
        tm = TurbomoleData.from_parser(parser)

        return cls(numerical_integration=ni, rotational=rot, vibrational=vib, run=run, tm=tm)


class Ricc2Output(BaseData):
    """
    Object containing the data of the output of a ricc2 calculation.
    """
    def __init__(self, mp2, run, tm):
        """
        Args:
            mp2 (MP2Results): MP2 results.
            run (RunData): information about calculation running (e.g. timings, ...)
            tm (TurbomoleData): information about the turbomole used for the calculation.
        """
        self.mp2 = mp2
        self.run = run
        self.tm = tm

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of Ricc2Output from a parser based on the stdout
        of a Turbomole executable.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            Ricc2Output.
        """
        mp2 = MP2Results.from_parser(parser)
        run = RunData.from_parser(parser)
        tm = TurbomoleData.from_parser(parser)

        return cls(mp2=mp2, run=run, tm=tm)


class MP2Output(BaseData):
    """
    Object containing the data of the output of an MP2 calculation, i.e. mpgrad, ricc2 (with proper MP2 options),
    or pnoccsd (with proper MP2 options).

    Note: Parsing of PNO-based MP2 calculations (i.e. performed with the pnoccsd program) is not (yet) supported.
    """
    def __init__(self, info, results, geometry, basis, symmetry, cosmo, run, tm):
        """
        Args:
            info (MP2Data): Information about MP2 calculation.
            results (MP2Results): Results of an MP2 calculation.
            geometry (GeometryData): the geometry of the system.
            basis (BasisData): the basis used for the calculation.
            symmetry (SymmetryData): information about the symmetry of the molecule.
            cosmo (CosmoData): cosmo approximations and results.
            run (RunData): information about calculation running (e.g. timings, ...)
            tm (TurbomoleData): information about the turbomole used for the calculation.
        """
        self.info = info
        self.results = results
        self.geometry = geometry
        self.basis = basis
        self.symmetry = symmetry
        self.cosmo = cosmo
        self.run = run
        self.tm = tm

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of MP2Output from a parser based on the stdout
        of a Turbomole executable.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            MP2Output.
        """
        info = MP2Data.from_parser(parser)
        results = MP2Results.from_parser(parser)
        run = RunData.from_parser(parser)
        tm = TurbomoleData.from_parser(parser)
        geometry = GeometryData.from_parser(parser)
        basis = BasisData.from_parser(parser)
        symmetry = SymmetryData.from_parser(parser)
        cosmo = CosmoData.from_parser(parser)
        # TODO: fix parsing of Cosmo data in non-scf executables or add a separate parsing

        return cls(info=info, results=results, geometry=geometry,
                   basis=basis, symmetry=symmetry, cosmo=cosmo, run=run, tm=tm)


class JobexOutput(BaseData):
    """
    Object containing the data of the output of the last step of a jobex calculation.
    Namely the outputs of the last energy, gradient and relax calculations stored in
    the "job.last" file.
    """

    def __init__(self, energy, gradient, relax):
        """
        Args:
            energy (ScfOutput): the output data of the energy calculation.
            gradient: the output data of the energy calculation. The type depends on
                the type of calculation. Can be GradOutput, EgradOutput or coming from
                other executables (e.g. Ricc2Output or MP2Output).
            relax (RelaxOutput or StatptOutput): the output data of the relax calculation.
        """
        self.energy = energy
        self.gradient = gradient
        self.relax = relax

    @classmethod
    def from_parser(cls, parser):
        """
        Generates an instance of JobexOutput from a parser based on the output
        contained in the "job.last" file generated by jobex. Splits the output
        in the different types of executables and generates suitable objects
        for them.

        Args:
            parser (Parser): the parser to be used to extract the data.

        Returns:
            JobexOutput.
        """
        exec_en, p_en, exec_grad, p_grad, exec_relax, p_relax = parser.get_split_jobex_parsers()

        energy = gradient = relax = None
        if exec_en in exec_to_out_obj:
            energy = exec_to_out_obj[exec_en].from_parser(p_en)

        if exec_grad in exec_to_out_obj:
            gradient = exec_to_out_obj[exec_grad].from_parser(p_grad)

        if exec_relax in exec_to_out_obj:
            relax = exec_to_out_obj[exec_relax].from_parser(p_relax)

        return cls(energy=energy, gradient=gradient, relax=relax)


# correspondence between name of the tm executable/name of the log file
# and the object used to parse the output.
# Does not include jobex since the output of jobex is in the "job.last" file.
exec_to_out_obj = {
    "dscf": ScfOutput,
    "ridft": ScfOutput,
    "escf": EscfOutput,
    "grad": GradOutput,
    "rdgrad": GradOutput,
    "egrad": EgradOutput,
    "relax": RelaxOutput,
    "statpt": StatptOutput,
    "aoforce": AoforceOutput,
    "force": AoforceOutput,  # because when executed aoforce has the name "force" in the header.
    "mpgrad": MP2Output,
    "ricc2": Ricc2Output,
    "rimp2": MP2Output,
}
