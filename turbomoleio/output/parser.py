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
Module with the main parsing utilities for the stdout of Turbomole executables.
"""

import re
import datetime
import numpy as np
from collections import namedtuple
from monty.functools import lazy_property


# common pattern for regex
date_format = "%Y-%m-%d %H:%M:%S.%f"
float_number_re = r"[+-]?[0-9]*[.]?[0-9]+"
float_number_d_re = r"[+-]?[0-9]*[.]?[0-9]+D[+-]\d{2}"
float_number_e_re = r"[+-]?[0-9]*[.]?[0-9]+E[+-]\d{2}"
float_number_all_re = r"[+-]?[0-9]*[.]?[0-9]+(?:[ED][+-]\d{2})?"
irrep_re_group = r"\w'\""


# pre compiled regex
asterisks_re = re.compile(r"\*+")


# strings present in the output for all the types of standard xc functionals.
# NB This should be updated if other standard xc functionals are added to TM.
functional_strings = {
    "s-vwn": r"Slater.+?Dirac.+?exchange.+?with.+?VWN.+?corr\..+?functional",
    "s-vwn_Gaussian": r"Slater.+?Dirac.+?exchange.+?with.+?VWN.+?corr\.functional.+?\(fit.+?III.+?as.+?in.+?Gaussian\)",
    "pwlda": r"Slater.+?Dirac.+?exchange.+?with.+?PW.+?LDA.+?corr\..+?functional",
    "b-lyp": r"B\-LYP.+?functional.+?exchange\:.+?LDA.+?\+.+?Becke.+?\(B88\).+?correlation\:.+?Lee\-Yang\-Parr.+?\(LYP\)",
    "b-vwn": r"B\-VWN.+?functional.+?exchange\:.+?LDA.+?\+.+?Becke.+?\(B88\).+?correlation\:.+?LDA.+?\(VWN\)",
    "b-p": r"B\-P86.+?functional.+?exchange\:.+?LDA.+?\+.+?Becke.+?\(B88\).+?correlation\:.+?LDA.+?\(VWN\).+?\+.+?Perdew.+?\(P86\)",
    "pbe": r"PBE.+?functional.+?exchange\:.+?LDA.+?\+.+?PBE.+?correlation\:.+?LDA.+?\(PW\).+?\+.+?PBE",
    "tpss": r"TPSS.+?meta\-GGA.+?functional.+?exchange\:.+?LDA.+?\+.+?TPSS.+?correlation\:.+?LDA.+?\(PW\).+?\+.+?TPSS",
    "bh-lyp": r"Becke\-Half\-and\-Half\-LYP.+?hybrid.+?functional\:.+?BH\-LYP.+?exchange\:.+?1\/2.+?\(LDA.+?\+.+?Becke.+?\(B88\)\).+?\+.+?1\/2.+?HF.+?correlation\:.+?Lee\-Yang\-Parr.+?\(LYP\)",
    "b3-lyp": r"Becke\-3\-Parameter.+?hybrid.+?functional\:.+?B3\-LYP.+?exchange\:.+?0\.8\*LDA.+?\+.+?0\.72\*B88.+?\+.+?0\.2\*HF.+?correlation\:.+?0\.19\*LDA\(VWN\).+?\+.+?0\.81\*LYP",
    "b3-lyp_Gaussian": r"Becke\-3\-Parameter.+?hybrid.+?functional\:.+?B3\-LYP.+?exchange\:.+?0\.8\*LDA.+?\+.+?0\.72\*B88.+?\+.+?0\.2\*HF.+?correlation\:.+?0\.19\*LDA\(VWNIII\).+?\+.+?0\.81\*LYP.+?\(VWNIII.+?fit.+?as.+?in.+?Gaussian\)",
    "pbe0": r"PBE0.+?hybrid.+?functional.+?exchange\:.+?3\/4.+?\(LDA\+PBE\).+?\+.+?1\/4.+?HF.+?correlation\:.+?LDA.+?\(PW\).+?\+.+?PBE",
    "tpssh": r"TPSS.+?global.+?hybrid.+?functional.+?exchange\:.+?9\/10.+?\(LDA\+TPSS\).+?\+.+?1\/10.+?HF.+?correlation\:.+?LDA.+?\(PW\).+?\+.+?TPSS",
    "pw6b95": r"PW6B95.+?global.+?meta.+?hybrid.+?functional.+?code.+?by.+?Stefan.+?Grimme\,.+?University.+?of.+?Muenster.+?Zhao.+?and.+?Truhlar\,.+?J.+?Phys.+?Chem.+?A\,.+?109\,.+?25\,.+?2005\,.+?5656\.",
    "m06": r"M06.+?meta\-GGA.+?functional.+?Truhlar.+?functional.+?with.+?27\%.+?HF.+?exchange.+?USING.+?XCfun.+?library\,.+?see.+?documentation.+?XCFun.+?library.+?is.+?being.+?used\,.+?version\:.+?1\.99000000000000.+?XCFun.+?DFT.+?library.+?Copyright.+?2009\-2011.+?Ulf.+?Ekstrom.+?and.+?contributors\..+?See.+?http\:\/\/admol\.org\/xcfun.+?for.+?more.+?information\..+?This.+?is.+?free.+?software\;.+?see.+?the.+?source.+?code.+?for.+?copying.+?conditions\..+?There.+?is.+?ABSOLUTELY.+?NO.+?WARRANTY\;.+?not.+?even.+?for.+?MERCHANTABILITY.+?or.+?FITNESS.+?FOR.+?A.+?PARTICULAR.+?PURPOSE\..+?For.+?details.+?see.+?the.+?documentation\..+?Scientific.+?users.+?of.+?this.+?library.+?should.+?cite.+?U\..+?Ekstrom\,.+?L\..+?Visscher\,.+?R\..+?Bast\,.+?A\..+?J\..+?Thorvaldsen.+?and.+?K\..+?Ruud\;.+?J\.Chem\.Theor\.Comp\..+?2010\,.+?DOI\:.+?10\.1021\/ct100117s.+?XCFun.+?uses.+?functional\:.+?m06",
    "m06-l": r"M06\-L.+?meta\-GGA.+?functional.+?Truhlar.+?functional.+?without.+?HF.+?exchange.+?USING.+?XCfun.+?library\,.+?see.+?documentation.+?XCFun.+?library.+?is.+?being.+?used\,.+?version\:.+?1\.99000000000000.+?XCFun.+?DFT.+?library.+?Copyright.+?2009\-2011.+?Ulf.+?Ekstrom.+?and.+?contributors\..+?See.+?http\:\/\/admol\.org\/xcfun.+?for.+?more.+?information\..+?This.+?is.+?free.+?software\;.+?see.+?the.+?source.+?code.+?for.+?copying.+?conditions\..+?There.+?is.+?ABSOLUTELY.+?NO.+?WARRANTY\;.+?not.+?even.+?for.+?MERCHANTABILITY.+?or.+?FITNESS.+?FOR.+?A.+?PARTICULAR.+?PURPOSE\..+?For.+?details.+?see.+?the.+?documentation\..+?Scientific.+?users.+?of.+?this.+?library.+?should.+?cite.+?U\..+?Ekstrom\,.+?L\..+?Visscher\,.+?R\..+?Bast\,.+?A\..+?J\..+?Thorvaldsen.+?and.+?K\..+?Ruud\;.+?J\.Chem\.Theor\.Comp\..+?2010\,.+?DOI\:.+?10\.1021\/ct100117s.+?XCFun.+?uses.+?functional\:.+?m06l",
    "m06-2x": r"M06\-2X.+?meta\-GGA.+?functional.+?Truhlar.+?functional.+?with.+?54\%.+?HF.+?exchange.+?USING.+?XCfun.+?library\,.+?see.+?documentation.+?XCFun.+?library.+?is.+?being.+?used\,.+?version\:.+?1\.99000000000000.+?XCFun.+?DFT.+?library.+?Copyright.+?2009\-2011.+?Ulf.+?Ekstrom.+?and.+?contributors\..+?See.+?http\:\/\/admol\.org\/xcfun.+?for.+?more.+?information\..+?This.+?is.+?free.+?software\;.+?see.+?the.+?source.+?code.+?for.+?copying.+?conditions\..+?There.+?is.+?ABSOLUTELY.+?NO.+?WARRANTY\;.+?not.+?even.+?for.+?MERCHANTABILITY.+?or.+?FITNESS.+?FOR.+?A.+?PARTICULAR.+?PURPOSE\..+?For.+?details.+?see.+?the.+?documentation\..+?Scientific.+?users.+?of.+?this.+?library.+?should.+?cite.+?U\..+?Ekstrom\,.+?L\..+?Visscher\,.+?R\..+?Bast\,.+?A\..+?J\..+?Thorvaldsen.+?and.+?K\..+?Ruud\;.+?J\.Chem\.Theor\.Comp\..+?2010\,.+?DOI\:.+?10\.1021\/ct100117s.+?XCFun.+?uses.+?functional\:.+?m062x",
    "lhf": r"Localized.+?Hartree\-Fock.+?Methods\:.+?F\..+?Della.+?Sala.+?and.+?A\..+?Goerling\,.+?J\..+?Chem\..+?Phys\..+?115\,.+?5718.+?\(2001\).+?F\..+?Della.+?Sala.+?and.+?A\..+?Goerling\,.+?J\..+?Chem\..+?Phys\..+?116\,.+?5374.+?\(2002\)",
    "oep": r"Exact\-Exchange.+?Optimized.+?Effective.+?Potential.+?Method\:.+?Hesselmann\,.+?A\.\,.+?Goetz\,.+?A\.\,.+?Della.+?Sala\,.+?F\.\,.+?Goerling\,.+?A\.\,.+?J\..+?Chem\..+?Phys\.\,.+?127.+?\(2007\)\,.+?054102",
    "b97-d": r"exchange\:.+?B97GGA.+?vdW.+?refit.+?correlation\:.+?\".+?\".+?\".+?S\..+?Grimme\,.+?J\.Comput\..+?Chem\..+?27\,.+?\(2006\)\,.+?1787\-1799",
    "pbeh-3c": r"PBE0.+?modified.+?by.+?S\..+?Grimme.+?for.+?D3.+?and.+?gCP.+?exchange\:.+?PBE.+?\(kappa\=1\.0245\,.+?mu\=0\.12345679\).+?correlation\:.+?LDA.+?\(PW\).+?\+.+?PBE",
    "b97-3c": r"exchange\:.+?B97GGA.+?vdW.+?refit.+?correlation\:.+?\".+?\".+?\".+?S\..+?Grimme\,.+?modifications.+?B97\-3c.+?\(2016\)",
    "lh07t-svwn": r"Lh07t\-SVWN.+?local.+?hybrid.+?functional.+?Local.+?hybrid.+?semi\-numerical.+?integral.+?thresholds\:.+?S\-junctions.+?\:.+?0\.10E\-05.+?P\-junctions.+?\:.+?0\.10E\-05",
    "lh07s-svwn": r"Lh07s\-SVWN.+?local.+?hybrid.+?functional.+?Local.+?hybrid.+?semi\-numerical.+?integral.+?thresholds\:.+?S\-junctions.+?\:.+?0\.10E\-05.+?P\-junctions.+?\:.+?0\.10E\-05",
    "lh12ct-ssirpw92": r"Lh12ct\-SsirPW92.+?local.+?hybrid.+?functional.+?Local.+?hybrid.+?semi\-numerical.+?integral.+?thresholds\:.+?S\-junctions.+?\:.+?0\.10E\-05.+?P\-junctions.+?\:.+?0\.10E\-05",
    "lh12ct-ssifpw92": r"Lh12ct\-SsifPW92.+?local.+?hybrid.+?functional.+?Local.+?hybrid.+?semi\-numerical.+?integral.+?thresholds\:.+?S\-junctions.+?\:.+?0\.10E\-05.+?P\-junctions.+?\:.+?0\.10E\-05",
    "lh14t-calpbe": r"Lh14t\-calPBE.+?local.+?hybrid.+?functional.+?Local.+?hybrid.+?semi\-numerical.+?integral.+?thresholds\:.+?S\-junctions.+?\:.+?0\.10E\-05.+?P\-junctions.+?\:.+?0\.10E\-05",
    "b2-plyp": r"Hybrid.+?part.+?of.+?B2\-PLYP.+?double.+?hybrid.+?functional.+?exchange\:.+?0\.47\(LDA.+?\+.+?Becke.+?\(B88\)\).+?\+.+?0\.53.+?HF.+?correlation\:.+?0\.73.+?LYP.+?\+.+?0\.27.+?PT2.+?\(MP2.+?program\).+?S\..+?Grimme\,.+?JCP.+?124\,.+?\(2006\)\,.+?034108\-16"
}


# type of xc functionals
# NB This should be updated if other standard xc functionals are added to TM.
functional_types = {
 's-vwn': 'LDA',
 's-vwn_Gaussian': 'LDA',
 'pwlda': 'LDA',
 'b-lyp': 'GGA',
 'b-vwn': 'GGA',
 'b-p': 'GGA',
 'pbe': 'GGA',
 'tpss': 'MGGA',
 'bh-lyp': 'HYB',
 'b3-lyp': 'HYB',
 'b3-lyp_Gaussian': 'HYB',
 'pbe0': 'HYB',
 'tpssh': 'MHYB',
 'pw6b95': 'MHYB',
 'm06': 'MHYB',
 'm06-l': 'MGGA',
 'm06-2x': 'MHYB',
 'lhf': 'ODFT',
 'oep': 'ODFT',
 'b97-d': 'GGA',
 'pbeh-3c': 'HYB',
 'b97-3c': 'GGA',
 'lh07t-svwn': 'LHYB',
 'lh07s-svwn': 'LHYB',
 'lh12ct-ssirpw92': 'LHYB',
 'lh12ct-ssifpw92': 'LHYB',
 'lh14t-calpbe': 'LHYB',
 'b2-plyp': 'DHYB'
}


def convert_float(f):
    """
    Helper function to convert a float from the output.
    Handles simple float and exponential notation with D and E
    (e.g 11.23 1.123E+01 1.123D+01).
    Returns None if the string is composed of only asterisks,
    as it might be in case of formatting problems in fortran.

    Args:
        f (str): the string to be converted to a float.

    Returns:
        float.
    """
    if asterisks_re.match(f.strip()):
        return None

    return float(f.replace("D", "E"))


def convert_int(i):
    """
    Helper function to convert an int from the output.
    Returns None if the string is composed of only asterisks,
    as it might be in case of formatting problems in fortran.

    Args:
        i (str): the string to be converted to an int.

    Returns:
        int.
    """
    if asterisks_re.match(i.strip()):
        return None

    return int(i)


def convert_time_string(l):
    """
    Helper function to convert the string of the cpu and wall time given in the
    Turbomole outputs to seconds.

    Example: 2 days 1 hours 12 minutes and 55 seconds

    Args:
        l (str): the line with the timings.

    Returns:
        float: the value of the time in seconds.
    """
    split = l.split()

    time = 0
    for i, s in enumerate(split):
        if s == "days":
            time += float(split[i-1]) * 86400
        if s == "hours":
            time += float(split[i-1]) * 3600
        if s == "minutes":
            time += float(split[i-1]) * 60
        if s == "seconds":
            time += float(split[i-1])

    return time


class Parser:
    """
    Main Parser object for Turbomole output files.

    This object allows to parse the different outputs (the stdout) produced by
    the Turbomole executables. Each method parses a specific portion of the output
    aiming at being as reusable as possible across the different types of outputs.
    The quantities are accessible as lazy properties, to avoid parsing the same
    sections multiple times.
    The parser expects the output of a completed calculation. Parsing of crashed
    calculations may lead to wrong output.
    The methods rely heavily on regex to select portions of the output and extract
    the quantities.
    """

    def __init__(self, string):
        """

        Args:
            string (str): the string of the output from Turbomole.
        """

        self._string = string

    @property
    def string(self):
        """
        The string given as an input.
        """
        return self._string

    @classmethod
    def from_file(cls, filepath, check_all_done=True):
        """
        Generates an instance from a file path.

        Args:
            filepath (str): path to the output file to read.
            check_all_done (bool): if True it will be checked that the job
                was "all done". If not an exception is raised.

        Returns:
            an instance of Parser.
        """

        try:
            with open(filepath, 'r') as f:
                string = f.read()
        except UnicodeDecodeError:
            with open(filepath, 'r', errors='ignore') as f:
                string = f.read()

        if check_all_done and string.rfind("all done") < 0:
            raise ValueError("The string does not contain data for a completed calculation")

        return cls(string=string)

    @lazy_property
    def all_done(self):
        """
        True if "all done" is present in the string.
        """
        return self.string.rfind("all done") >= 0

    @lazy_property
    def header(self):
        """
        The information contained in the header of the output, like the TM version and timings.
        Valid for all the TM executables.

        Returns:
            dict with "executable" name, "host" of execution, "tm_version", "tm_build", "start_time".
        """
        # Example string parsed here:
        #  escf (node001) : TURBOMOLE V7.3 ( 22118 ) 1 Jul 2018 at 20:38:15
        #  Copyright (C) 2018 TURBOMOLE GmbH, Karlsruhe
        #    2018-11-13 09:30:15.283
        r = r"^.*?([\w]+)\s+\((.*?)\).*?TURBOMOLE(.*?)(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}.\d{3})"
        match = re.search(r, self.string, re.DOTALL)
        if not match:
            return None

        # the section of the version can miss some of the required information.
        # In this case turbomole version and build will be set to None.
        # For example it may be:
        # escf (node001) : TURBOMOLE rev. compiled 1 Jul 2018 at 20:38:15
        r_version = r'V([\d\.]+\d)\s+.*'
        r_build = r'V[\d\.]+\d\s+\((.*)\)'
        match_version = re.search(r_version, match.group(3).strip())
        match_build = re.search(r_build, match.group(3).strip())
        if match_version:
            tm_version = match_version.group(1).strip()
        else:
            tm_version = None
        if match_build:
            tm_build = match_build.group(1).strip() or None  # When the build info inside the parenthesis is empty
        else:
            tm_build = None

        d = dict(executable=match.group(1),
                 host=match.group(2).strip(),
                 tm_version=tm_version,
                 tm_build=tm_build,
                 start_time=datetime.datetime.strptime(match.group(4), date_format))

        return d

    @lazy_property
    def centers(self):
        """
        Center of mass and center of charge values.
        Valid for most of the TM executables (including the scf, escf, grad).

        Returns:
            dict with "center_of_mass" and "center_of_charge".
        """
        d = dict(center_of_mass=None, center_of_charge=None)

        r = r"center of nuclear mass[\s:]*" + (r"\s*("+float_number_re+r")")*3

        match = re.search(r, self.string)

        if match:
            d["center_of_mass"] = [convert_float(match.group(i)) for i in range(1,4)]

        r = r"center of nuclear charge[\s:]*" + (r"\s*("+float_number_re+r")")*3

        match = re.search(r, self.string)

        if match:
            d["center_of_charge"] = [convert_float(match.group(i)) for i in range(1,4)]

        return d

    @lazy_property
    def coordinates(self):
        """
        Coordinates, species and charges of the atoms.
        Valid for most of the TM executables (including the scf, escf, grad).

        Returns:
            dict with "coord", "species", "charges", and "isotopes" as lists.
        """
        # Example string parsed here:
        #               +--------------------------------------------------+
        #               | Atomic coordinate, charge and isotop information |
        #               +--------------------------------------------------+
        #
        #                     atomic coordinates            atom    charge  isotop
        #           0.00000000    0.00000000    1.00000000    h      1.000     0
        #           0.00000000    0.00000000   -1.00000000    h      1.000     0
        r = r"atomic coordinates\s+atom\s+charge\s+isotop\s+(.+?)\s+center of nuclear mass"

        match = re.search(r, self.string, re.DOTALL)
        if not match:
            return None

        coords = []
        species = []
        charges = []
        isotopes = []

        for l in match.group(1).splitlines():
            l = l.strip()
            if not l:
                continue

            split = l.split()
            coords.append([convert_float(f) for f in split[:3]])
            species.append(split[3])
            charges.append(convert_float(split[4]))
            isotopes.append(convert_int(split[5]))

        return dict(coords=coords, species=species, charges=charges, isotopes=isotopes)

    @lazy_property
    def basis(self):
        """
        Data of the basis used for the calculation. Also the auxiliary base in
        case of RI.
        Valid for most of the TM executables (including the scf, escf, grad).

        Returns:
            dict with "basis_per_specie" a dict with species as keys and name of the
            basis as values, "aux_basis_per_specie" the same for auxiliary basis,
            "number_scf_basis_func" and "number_scf_aux_basis_func".
        """
        # See a test file as an example. Information parsed from the section starting with:
        #      +--------------------------------------------------+
        #      |               basis set information              |
        #      +--------------------------------------------------+
        r = r"basis set information.*?symmetry group of the molecule"

        match = re.search(r, self.string, re.DOTALL)

        if match is None:
            return None

        basis_string = match.group()

        # takes the section:
        #    type   atoms  prim   cont   basis
        #    ---------------------------------------------------------------------------
        #     h        2     10      5   def-SV(P)   [2s1p|4s2p]
        #    ---------------------------------------------------------------------------
        r_atoms = r"type   atoms  prim   cont   basis\s*-{5,}(.*?)-{5,}"

        match_atoms = re.search(r_atoms, basis_string, re.DOTALL)

        basis_per_specie = {}
        if match_atoms:
            for l in match_atoms.group(1).splitlines():
                l = l.strip()
                if not l:
                    continue

                l_split = l.split()
                basis_per_specie[l_split[0]] = l_split[4]

        # For ridft can be present the section starting with "RI-J AUXILIARY BASIS SET information".
        r_atoms_aux = r"RI-J[K]? AUXILIARY BASIS SET.*?" + r_atoms

        match_atoms_aux = re.search(r_atoms_aux, basis_string, re.DOTALL)

        aux_basis_per_specie = {}
        if match_atoms_aux:
            for l in match_atoms_aux.group(1).splitlines():
                l = l.strip()
                if not l:
                    continue
                l_split = l.split()
                aux_basis_per_specie[l_split[0]] = l_split[4]

        # search the total number of basis function. This should be present at maximum two
        # one for the standard basis and one for the auxiliary.
        list_n_basis_func = re.findall(r"total number of SCF-basis functions[\s:]+?(\d+)", basis_string)

        len_n_basis = len(list_n_basis_func)
        if len_n_basis > 2:
            raise RuntimeError("Found {} instances of 'total number of SCF-basis functions', while "
                               "expecting at most 2.".format(len_n_basis))

        if len_n_basis > 0:
            number_scf_basis_func = convert_int(list_n_basis_func[0])
        else:
            number_scf_basis_func = None

        if len_n_basis == 2:
            number_scf_aux_basis_func = convert_int(list_n_basis_func[1])
        else:
            number_scf_aux_basis_func = None

        d = dict(basis_per_specie=basis_per_specie,
                 aux_basis_per_specie=aux_basis_per_specie,
                 number_scf_basis_func=number_scf_basis_func,
                 number_scf_aux_basis_func=number_scf_aux_basis_func)
        return d

    @lazy_property
    def symmetry(self):
        """
        Information about the symmetry of the molecule and the irreducible representations.
        Valid for all the TM executables (in some case only a part of the information might
        be available).

        Returns:
            dict with symmetry "symbol", "n_reps" and list of representation symbols "reps".
        """

        match_symm = re.search(r"symmetry group of the molecule[\s:]+?([\w]+)", self.string)
        if match_symm:
            mol_sym_group = match_symm.group(1)
        else:
            mol_sym_group = None

        # Example: there are 4 real representations :   a1   a2   b1   b2
        match_reps = re.search(r"there are\s+(\d+)\s+real representations", self.string)
        if match_reps:
            n_reps = convert_int(match_reps.group(1))
            match_reps = re.search(r"there are\s+\d+\s+real representations\s*:(\s+[\w\'\"]+){"+str(n_reps)+"}",
                                   self.string, re.MULTILINE)
            if not match_reps:
                raise RuntimeError("Could not parse correctly the list of representations")
            reps = match_reps.group().split()[-n_reps:]
        else:
            n_reps = None
            reps = None

        if all(x is None for x in (mol_sym_group, n_reps, reps)):
            return None

        d = dict(symbol=mol_sym_group,
                 n_reps=n_reps,
                 reps=reps)

        return d

    @lazy_property
    def cosmo_header(self):
        """
        Information from the header of the cosmo section (area and volume)
        when cosmo is activated. Valid for scf and escf/egrad.

        Returns:
            dict with "area" and "volume".
        """
        # See a test file as an example. Selects the section of the cosmo header starting with:
        # ==============================================================================
        #                       COSMO switched on
        # ==============================================================================
        r = r"COSMO switched on.*?(---)+"

        match = re.search(r, self.string, re.DOTALL)

        if match is None:
            return None

        d = dict(area=None,
                 volume=None)

        for l in match.group().splitlines():
            if "area" in l:
                d["area"] = convert_float(l.split()[-1])
            elif "volume" in l:
                d["volume"] = convert_float(l.split()[-1])

        return d

    @lazy_property
    def density_functional_data(self):
        """
        Information relative to DFT calculations (xf functional and grids).
        Valid for scf, gradient and relax executables.

        Returns:
            dict with the full string describing the functional "functional_msg",
            "functional_name", "functional_type", the version of "xcfun",
            "spherical_gridsize, number of "gridpoints" for spherical integration.
        """
        # See a test file as an example. Selects the section of dft information starting with:
        #           ------------------
        #           density functional
        #           ------------------
        r = r"density functional[\s]*?-{5,}.*?(biggest|-{5,})"

        match = re.search(r, self.string, re.DOTALL)

        if match is None:
            return None

        dft_string = match.group()

        d = dict(functional_msg=None,
                 functional_name=None,
                 functional_type=None,
                 xcfun=None,
                 spherical_gridsize=None,
                 gridpoints=None)

        # capture the part with the message describing the xc functional.
        # "iterations will be done with small grid" is present only for modified grids (e.g. m3)
        r_func = r"--\s+(.*?)\s*(iterations will be done with small grid|spherical integration)"

        if "USE AT YOUR OWN RISK" not in dft_string:
            func_match = re.search(r_func, dft_string, re.DOTALL)

            # Matches the type of functional based on the message produced by TM.
            # The messages have been extracted for all the standard types available
            # in define. Notice that if further xc functional need be parsed, the
            # functional_strings dictionary should be updated accordingly.
            # Also arbitrary mix will be ignored (when "USE AT YOUR OWN RISK" is present).
            if func_match:
                func_msg = func_match.group(1).strip()
                d["functional_msg"] = func_msg

                for func_name, desc in functional_strings.items():
                    if re.search(desc, func_msg, re.DOTALL):
                        d["functional_name"] = func_name
                        d["functional_type"] = functional_types[func_name]
                        break

        for l in dft_string.splitlines():
            if "XCFun" in l and "version" in l:
                d["xcfun"] = l.split()[-1]

            if "spherical gridsize" in l:
                d["spherical_gridsize"] = convert_int(l.split()[-1])

            if "i.e. gridpoints" in l:
                d["gridpoints"] = convert_int(l.split()[-1])

        return d

    @lazy_property
    def rij_info(self):
        """
        Information about RI. Valid for ridft and escf/egrad.

        Returns:
            dict with "marij" to True if marij calculation, "rij_memory", "rik" to True if $rik
            and "ricore" memory in MB.
        """
        # Parses the section that can start with RI, RI-J or RI-JK
        #            ------------------------
        #                RI-JK - INFORMATION
        #            ------------------------
        r = r"RI[\-JK]* - INFORMATION.*?Memory allocated for.*?$"

        match = re.search(r, self.string, re.DOTALL|re.MULTILINE)

        if match is None:
            return None

        match_str = match.group()

        d = dict(marij=False,
                 rij_memory=None,
                 rik=False,
                 ricore=None)

        for l in match_str.splitlines():
            if "Multipole Accelerated RI-J" in l:
                d["marij"] = True

            # stopping at "RI" should match the line for both ridft and escf
            if "Memory allocated for RI" in l:
                d["rij_memory"] = convert_int(l.split()[-2])

        match_rik_1 = re.search("BLOCKING OF .*?MOS FOR RI-K", self.string, re.DOTALL) is not None # ridft
        match_rik_2 = "FOUND RI-K FLAG" in self.string # escf
        d["rik"] = match_rik_1 or match_rik_2

        # ricore in ridft (is in another section)
        r_ricore = r"Allocatable memory for RI due to \$ricore \(MB\):\s+(\d+)"
        match = re.search(r_ricore, self.string)
        if match:
            d["ricore"] = convert_int(match.group(1))
        else:
            # ricore in escf
            r_ricore = r"Core memory available \(ricore\)\s+(\d+)"
            match = re.search(r_ricore, match_str)
            if match:
                d["ricore"] = convert_int(match.group(1))

        return d

    @lazy_property
    def dftd(self):
        """
        Information about dispersion correction in dft.
        Valid for scf executables.

        Returns:
            dict with the type of "correction" (string) and the value of correction "en_corr".
        """

        # assume that the DFT-D part is between these two sections
        r = r"nuclear repulsion energy(.*?)-S,T\+V- integrals"

        match = re.search(r, self.string, re.DOTALL)

        if match is None:
            return None

        match_str = match.group(1)

        if "DFT-D" not in match_str:
            return None

        d = dict(correction=None,
                 en_corr=None)

        # The output changes depending on the type of correction. Is not always mentioned.
        # For D1 and D2 should be deduced from the value of the damping.
        if "DFT-D V3(BJ)" in match_str:
            d["correction"] = "D3-BJ"
        elif "DFT-D V3" in match_str:
            d["correction"] = "D3"
        else:
            damp_l = re.search(r"exponent of damping function \(d\)\s+(" + float_number_re + r")", match_str)
            if damp_l:
                damp_val = convert_float(damp_l.group(1))
                if damp_val == 20:
                    d["correction"] = "D2"
                elif damp_val == 23:
                    d["correction"] = "D1"

        corr_l = re.search(r"empirical dispersive energy correction[\s=]+(" + float_number_re + r")", match_str)

        if corr_l:
            d["en_corr"] = convert_float(corr_l.group(1))

        return d

    @lazy_property
    def pre_scf_run(self):
        """
        Information extracted about the scf calculation printed before running
        the scf loop. Valid only for scf executables.

        Returns:
            dict with "diis" True if activated, diis_error_vect", "conv_tot_en", "conv_one_e_en",
            "virtual_orbital_shift_on", "virtual_orbital_shift_limit", "orbital_characterization",
            "restart_file", "n_occupied_orbitals".
        """

        # large portion to include several sections
        r = r"basis set information(.*?)ITERATION  ENERGY"

        match = re.search(r, self.string, re.DOTALL)

        if match is None:
            return None

        match_str = match.group(1)

        d = dict(diis=False,
                 diis_error_vect=None,
                 conv_tot_en=None,
                 conv_one_e_en=None,
                 virtual_orbital_shift_on=False,
                 virtual_orbital_shift_limit=None,
                 orbital_characterization=None,
                 restart_file=None,
                 n_occupied_orbitals=None)

        # analyze the string line by line to find those corresponding to the values to parse.
        for l in match_str.splitlines():
            if "DIIS switched on" in l:
                d["diis"] = True
                d["diis_error_vect"] = l.split()[-1]

            if "increment of total energy" in l:
                d["conv_tot_en"] = convert_float(l.split()[-1])

            if "increment of one-electron energy" in l:
                d["conv_one_e_en"] = convert_float(l.split()[-1])

            if "automatic virtual orbital shift switched on" in l:
                d["virtual_orbital_shift_on"] = True

            if "shift if e(lumo)-e(homo)" in l:
                d["virtual_orbital_shift_limit"] = convert_float(l.split()[-1])

            if "orbital characterization" in l:
                d["orbital_characterization"] = l.split()[-1]

            if "DSCF restart information will be dumped onto file" in l:
                d["restart_file"] = l.split()[-1]

            if "number of occupied orbitals" in l:
                d["n_occupied_orbitals"] = convert_int(l.split()[-1])

        return d

    @lazy_property
    def scf_iterations(self):
        """
        Data for each iteration of the scf loop.
        Valid only for scf executables.

        Returns:
            dict with the list of values extracted for each scf step "energies", "dampings".
            The value of the "first_index" and "n_steps". "converged" True if the calculation is converged.
        """

        # parsing an iteration of a block like:
        #                                               current damping :  0.700
        #  ITERATION  ENERGY          1e-ENERGY        2e-ENERGY     RMS[dD(SAO)]   TOL
        #    1  -76.143808666911    -124.43866883     38.987705388    0.000D+00 0.175D-08

        std_float = r"\s+(" + float_number_re + ")"
        d_float = r"\s+(" + float_number_d_re + ")"
        r = r"^\s*(\d+)" + std_float*3 + d_float*2+r"\s*$"

        match_l = re.findall(r, self.string, re.MULTILINE)

        if not match_l:
            return None

        iter_indices = []
        energies = []
        for l in match_l:
            iter_indices.append(convert_int(l[0]))
            energies.append(convert_float(l[1]))

        # dscf has a ":" while ridft as a "="
        dampings = re.findall(r"current damping[\s:=]+"+std_float, self.string)
        if dampings:
            dampings = [convert_float(d) for d in dampings]

        converged = "convergence criteria satisfied after " in self.string

        d = dict(energies=energies,
                 first_index=iter_indices[0],
                 n_steps=len(iter_indices),
                 dampings=dampings,
                 converged=converged)

        return d

    @lazy_property
    def scf_energies(self):
        """
        Final energies for scf calculation.
        Valid only for scf executables.

        Returns:
            dict with the energy and its contributions: "total_energy", "kinetic_energy",
            "potential_energy", "virial_theorem", "wavefunction_norm".
        """
        # Parsing the section:
        #              ------------------------------------------
        #              |  total energy      =    -76.34276215784  |
        #              ------------------------------------------
        #              :  kinetic energy    =     75.98641545802  :
        #              :  potential energy  =   -152.32917761586  :
        #              :  virial theorem    =      1.99533227919  :
        #              :  wavefunction norm =      1.00000000000  :
        #              ..........................................
        r = r"\|\s+total energy(.*?)\.{3,}"

        match = re.search(r, self.string, re.DOTALL)

        if match is None:
            return None

        d = dict(total_energy=None,
                 kinetic_energy=None,
                 potential_energy=None,
                 virial_theorem=None,
                 wavefunction_norm=None)

        for l in match.group().splitlines():
            for k in d.keys():
                if " ".join(k.split("_")) in l:
                    d[k] = convert_float(l.split()[-2])
                    break

        return d

    @lazy_property
    def cosmo_results(self):
        """
        Results of cosmo.
        Valid only for scf executables.

        Returns:
            dict with several subdictionaries: "parameters" ("nppa", "nspa", "nsph", "npspher", "disex", "disex2",
            "rsolv", "routf", "phsran", "ampran", "cavity", "epsilon", "refind", "fepsi"), "screening_charge"
            ("cosmo", "correction", "total"), "energies" ("total_energy", "total_energy_oc_corr",
            "dielectric_energy", "dielectric_energy_oc_corr"), "element_radius" ("radius, "sites").
        """
        # See test file for an example. Parses the section starting with:
        #  ==============================================================================
        #                                   COSMO RESULTS
        #  ==============================================================================
        # and splits in block for the section "CAVITY VOLUME/AREA", "SCREENING CHARGE", "ENERGIES", "ELEMENT RADIUS"
        r = r"COSMO RESULTS\s*(===)+.*?PARAMETER:(.*?)CAVITY VOLUME/AREA(.*?)SCREENING CHARGE:(.*?)ENERGIES(.*?)ELEMENT RADIUS \[A\]: ATOM LIST(.*?)(===)+"

        match = re.search(r, self.string, re.DOTALL)

        if match is None:
            return None

        # In this section a list of keyword_name: value are parsed, one for each line.
        # this dictionary allows to loop over all the lines and cast the value to the
        # correct type.
        key_types_parameters = dict(nppa=int,
                                    nspa=int,
                                    nsph=int,
                                    nps=int,
                                    npspher=int,
                                    disex=float,
                                    disex2=float,
                                    rsolv=float,
                                    routf=float,
                                    phsran=float,
                                    ampran=float,
                                    cavity=str,
                                    epsilon=str,
                                    refind=float,
                                    fepsi=float)

        d_parameters = {k: None for k in key_types_parameters.keys()}

        for l in match.group(2).splitlines():
            for k, t in key_types_parameters.items():
                # NB here the regex is tailored to explicitly capture the string that
                # need to be parsed. For example:
                # disex2: 3538.50
                # rsolv [A]: 1.3000
                # if more are added the regex may need an update.
                if re.search(k+r"[\s\[A\]]*:", l):
                    d_parameters[k] = t(l.split()[-1])
                    break

        d_screening = dict(cosmo=None,
                           correction=None,
                           total=None)

        for l in match.group(4).splitlines():
            if "cosmo" in l:
                d_screening["cosmo"] = convert_float(l.split()[-1])
            if "correction" in l:
                d_screening["correction"] = convert_float(l.split()[-1])
            if "total" in l:
                d_screening["total"] = convert_float(l.split()[-1])

        d_energies = dict(total_energy=None,
                          total_energy_oc_corr=None,
                          dielectric_energy=None,
                          dielectric_energy_oc_corr=None)

        for l in match.group(5).splitlines():
            if re.search(r"Total energy\s*=", l):
                d_energies["total_energy"] = convert_float(l.split()[-1])
            if re.search(r"Total energy \+ OC corr\.\s*=",l):
                d_energies["total_energy_oc_corr"] = convert_float(l.split()[-1])
            if re.search(r"Dielectric energy\s*=", l):
                d_energies["dielectric_energy"] = convert_float(l.split()[-1])
            if re.search(r"Diel\. energy \+ OC corr\.\s+=", l):
                d_energies["dielectric_energy_oc_corr"] = convert_float(l.split()[-1])

        d_element_radius = {}

        for l in match.group(6).splitlines():
            if not l.strip():
                continue

            split = l.split()
            # list of sites are separated by "," and can be defined as intervals with "-"
            # example: 1,3-5,7
            sites = []
            for s in split[2].split(","):
                if "-" in s:
                    sites_split = s.split("-")
                    sites.extend(range(convert_int(sites_split[0]), convert_int(sites_split[1]) + 1))
                else:
                    sites.append(convert_int(s))
            d_element_radius[split[0]] = {"radius": convert_float(split[1].replace(":", "")),
                                          "sites": sites}

        d = dict(parameters=d_parameters, screening_charge=d_screening,
                 energies=d_energies, element_radius=d_element_radius)

        return d

    @lazy_property
    def electrostatic_moments(self):
        """
        Information about the electrostatic moment (dipole and quadrupole).
        Valid only for scf executables.

        Returns:
            dict with single "dipole" and "quadrupole" moments.
        """
        # See test files for an example. Parses the section starting with:
        #  ==============================================================================
        #                            electrostatic moments
        #  ==============================================================================
        r = r"electrostatic moments\s*(==)+.*?charge\s*-+(.*?)-+\s*dipole moment(.*?)" \
            r"quadrupole moment(.*?anisotropy=\s*" + float_number_all_re + r")"
        # does not stop at "============" since there might be the PABOON section before that.

        match = re.search(r, self.string, re.DOTALL)

        if match is None:
            return None

        d = {"unrestricted_electrons": None}

        for l in match.group(2).splitlines():
            if not l.strip():
                continue
            elif "a-b" in l:
                d["unrestricted_electrons"] = convert_float(l.split()[-1])
            else:
                d["charge"] = convert_float(l.split()[-1])

        dipole = [0, 0, 0]
        d_dipole = {}
        for l in match.group(3).splitlines():
            if " x " in l:
                dipole[0] = convert_float(l.split()[-1])
            elif " y " in l:
                dipole[1] = convert_float(l.split()[-1])
            elif " z " in l:
                dipole[2] = convert_float(l.split()[-1])
            elif "| dipole moment |" in l:
                d_dipole["norm"] = convert_float(l.split()[5])

        d_dipole["moment"] = dipole

        quadrupole = np.zeros((3,3))
        d_quadrupole = {}
        for l in match.group(4).splitlines():
            if " xx " in l:
                quadrupole[0, 0] = convert_float(l.split()[-1])
            elif " yy " in l:
                quadrupole[1, 1] = convert_float(l.split()[-1])
            elif " zz " in l:
                quadrupole[2, 2] = convert_float(l.split()[-1])
            elif " xy " in l:
                quadrupole[0, 1] = quadrupole[1, 0] = convert_float(l.split()[-1])
            elif " xz " in l:
                quadrupole[0, 2] = quadrupole[2, 0] = convert_float(l.split()[-1])
            elif " yz " in l:
                quadrupole[1, 2] = quadrupole[2, 1] = convert_float(l.split()[-1])
            elif "trace=" in l:
                d_quadrupole["trace"] = 3 * convert_float(l.split()[-1])
            elif "anisotropy" in l:
                d_quadrupole["anisotropy"] = convert_float(l.split()[-1])

        d_quadrupole["moment"] = quadrupole.tolist()

        d["dipole"] = d_dipole
        d["quadrupole"] = d_quadrupole

        return d

    @lazy_property
    def timings(self):
        """
        Cpu and wall time for the calculation.
        Valid for all the Turbomole executables.

        Returns:
            dict with "cpu_time" and "wall_time" in seconds.
        """
        #          total  cpu-time :   0.74 seconds
        #          total wall-time :   0.79 seconds
        r = r"(total\s+cpu-time.*?)(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}.\d{3})"

        match = re.search(r, self.string, re.DOTALL)
        if not match:
            return None

        d = dict(cpu_time=None,
                 wall_time=None,
                 end_time=datetime.datetime.strptime(match.group(2), date_format))

        for l in match.group(1).splitlines():
            if "cpu-time" in l:
                d["cpu_time"] = convert_time_string(l.split(":")[1])
            if "wall-time" in l:
                d["wall_time"] = convert_time_string(l.split(":")[1])

        return d

    @lazy_property
    def s2(self):
        """
        Value S^2 of the spin.
        Valid only for scf executables.

        Returns:
            dict containing "s2".
        """
        # Parsing:  <S*S>     0.75205422
        r = r"<S\*S>\s+?(" + float_number_re + r")"

        match = re.search(r, self.string)

        if not match:
            return None

        return {"s2": convert_float(match.group(1))}

    @lazy_property
    def is_uhf(self):
        """
        True if UHF calculation.
        Valid for scf, gradient and escf executables.

        Returns:
            bool.
        """

        # "mode" appears in dscf and escf, while "modus" in ridft and grad.
        is_uhf = re.search("UHF (mode|modus) switched on", self.string) is not None
        return is_uhf

    @lazy_property
    def fermi(self):
        """
        Information about smearing of occupations.
        Valid only for scf executables.

        Returns:
            dict with "initial_elec_temp", "final_elec_temp", "annealing_homo_lumo_gap_limit",
            "smearing_de_limit".
        """

        # parses the section:
        #     Fermi smearing switched on
        #       Initial electron temperature:   500.00
        #         Final electron temperature:    50.00
        #                   Annealing factor:    0.900
        #        Annealing if HOMO-LUMO gap <  0.2E+00
        #       Smearing switched off if DE <  0.1E-02

        # this shows up before the start of the scf iterations. Since the beginning of that
        # may change (at least depending if it's a new calculation or a restart)
        # the end of the regex is on the last line and the value.
        r = r"Fermi smearing switched on.*?Smearing switched off if DE\s+<\s+" + float_number_e_re

        match = re.search(r, self.string, re.DOTALL)

        if not match:
            return None

        d = dict(initial_elec_temp=None,
                 final_elec_temp=None,
                 annealing_homo_lumo_gap_limit=None,
                 smearing_de_limit=None)

        for l in match.group().splitlines():
            if "Initial electron temperature" in l:
                d["initial_elec_temp"] = convert_float(l.split()[-1])
            if "Final electron temperature" in l:
                d["final_elec_temp"] = convert_float(l.split()[-1])
            if "Annealing factor" in l:
                d["annealing_factor"] = convert_float(l.split()[-1])
            if "Annealing if HOMO-LUMO gap" in l:
                d["annealing_homo_lumo_gap_limit"] = convert_float(l.split()[-1])
            if "Smearing switched off if DE" in l:
                d["smearing_de_limit"] = convert_float(l.split()[-1])

        return d

    @lazy_property
    def integral(self):
        """
        Thresholds for integrals.
        Valid for scf, gradient and escf executables.

        Returns:
            dict with "integral_neglect_threshold", "thize", "thime".
        """

        # parses the section
        #  integral neglect threshold       :  0.18E-08
        #  integral storage threshold THIZE :  0.10E-04
        #  integral storage threshold THIME :         5

        d = dict(integral_neglect_threshold=None,
                 thize=None,
                 thime=None)

        match = re.search(r"integral neglect threshold.*?$", self.string, re.MULTILINE)
        if match is not None:
            d["integral_neglect_threshold"] = convert_float(match.group().split()[-1])

        match = re.search(r"integral storage threshold THIZE.*?$", self.string, re.MULTILINE)
        if match is not None:
            d["thize"] = convert_float(match.group().split()[-1])

        match = re.search(r"integral storage threshold THIME.*?$", self.string, re.MULTILINE)
        if match is not None:
            d["thime"] = convert_int(match.group().split()[-1])

        if all(v is None for v in d.values()):
            return None

        return d

    @lazy_property
    def pre_escf_run(self):
        """
        Information extracted about the escf calculation printed before running
        the escf iterations. Valid for escf and egrad executables.

        Returns:
            dict with "calc_type", "residuum_convergence_criterium", "n_occupied_orbitals",
            "orbital_characterization", "max_davidson_iter", "machine_precision", "max_core_mem",
            "max_cao_basis_vectors", "max_treated_vectors", "irrep_data".
        """
        # See test files for examples. The values spread over a relatively large number of lines
        # and the lines are matched individually.

        # large portion of the code considered to accomodate the various pieces of information.
        r = r"basis set information.+Iteration\s+IRREP\s+Converged"

        match = re.search(r, self.string, re.DOTALL)

        if not match:
            return None

        match_str = match.group()

        d = dict(calc_type=None,
                 residuum_convergence_criterium=None,
                 n_occupied_orbitals=None,
                 orbital_characterization=None,
                 max_davidson_iter=None,
                 machine_precision=None,
                 max_core_mem=None,
                 max_cao_basis_vectors=None,
                 max_treated_vectors=None,
                 irrep_data=None)

        match = re.search(r"^(.*)[\-\s]CALCULATION\s+", match_str, re.MULTILINE)
        if match is not None:
            d["calc_type"] =match.group(1).strip()

        match = re.search(r"residuum convergence criterium.*?$", match_str, re.MULTILINE)
        if match is not None:
            d["residuum_convergence_criterium"] = convert_float(match.group().split()[-1])

        match = re.search(r"number of occupied orbitals.*?$", match_str, re.MULTILINE)
        if match is not None:
            d["n_occupied_orbitals"] = convert_int(match.group().split()[-1])

        match = re.search(r"orbital characterization.*?$", match_str, re.MULTILINE)
        if match is not None:
            d["orbital_characterization"] = match.group().split()[-1]

        match = re.search(r"maximum number of Davidson iterations set to.*?$", match_str, re.MULTILINE)
        if match is not None:
            d["max_davidson_iter"] = convert_int(match.group().split()[-1])

        match = re.search(r"machine precision.*?$", match_str,re.MULTILINE)
        if match is not None:
            d["machine_precision"] = convert_float(match.group().split()[-1])

        match = re.search(r"maximum core memory set to.*?$", match_str, re.MULTILINE)
        if match is not None:
            d["max_core_mem"] = convert_float(match.group().split()[-2])

        match = re.search(r"corresponding to\s+(\d+)\s+vectors in CAO basis", match_str)
        if match is not None:
            d["max_cao_basis_vectors"] = convert_int(match.group(1))

        match = re.search(r"maximum number of simultaneously treated vectors.*?$", match_str, re.MULTILINE)
        if match is not None:
            d["max_treated_vectors"] = convert_int(match.group().split()[-1])

        r = r"IRREP\s+tensor\s+space\s+dimension\s+number\s+of\s+roots(.+)maximum number of Davidson"
        match = re.search(r, match_str, re.DOTALL)
        if match is not None:
            irrep_data = {}
            for l in match.group(1).splitlines():
                l = l.strip()
                if not l:
                    continue
                split = l.split()
                irrep_data[split[0]] = {'tensor_space_dim': convert_int(split[1]), 'n_roots': convert_int(split[2])}

            d["irrep_data"] = irrep_data

        return d

    @lazy_property
    def escf_iterations(self):
        """
        Data for each iteration of the escf loop.
        Valid for escf and egrad executables.

        Returns:
            dict "converged" to True if the calculation is converged and a list of "steps". One element of the list for
                each escf step. One element of the sublist for each irrep. The inner list contains the convergence information:
                [name of the irrep, number of converged roots, euclidean residual norm].
        """
        # splits the different iterations, parsing the section starting with
        #  Iteration IRREP Converged      Max. Euclidean
        #                  roots          residual norm
        #
        #     1       a1       0        3.367862225107080D-01
        #             a2       9        6.722471267583007D-14

        r = r"Iteration\s+IRREP\s+Converged\s+Max\.\s+Euclidean\s+roots\s+residual\s+norm\s+(.*?)(converged|Warning)"
        match = re.search(r, self.string, re.DOTALL)

        if not match:
            return None

        iter_str = match.group(1)

        def convert_line(ls):
            ls[1] = convert_int(ls[1])
            ls[2] = convert_float(ls[2])
            return ls

        iter_steps_list = []
        iter_step = None
        for l in iter_str.splitlines():
            l = l.strip()
            if not l:
                continue

            split = l.split()
            if len(split) == 4:
                if iter_step is not None:
                    iter_steps_list.append(iter_step)
                iter_step = [convert_line(split[1:])]
            elif len(split) == 3:
                iter_step.append(convert_line(split))
            else:
                raise RuntimeError("wrong line {} in escf iterations".format(l))
        if iter_step:
            iter_steps_list.append(iter_step)

        converged = "converged" in match.group(2)

        d = dict(steps=iter_steps_list,
                 converged=converged)

        return d

    @lazy_property
    def escf_gs_total_en(self):
        """
        Value of the final total energy of escf.
        Valid for escf and egrad executables.

        Returns:
            float.
        """
        # Parses:
        #                             Ground state
        #  Total energy:                           -76.34301407423000
        r = r"Ground\s+state\s+Total\s+energy:\s+("+float_number_re+r")"

        match = re.search(r, self.string, re.DOTALL)

        if not match:
            return None
        else:
            return convert_float(match.group(1))

    @lazy_property
    def escf_excitations(self):
        """
        All the values for each calculated excitation.
        Valid for escf and egrad executables.

        Returns:
            dict with name of irreps as keywords and list of excitations as values. For each excitation
            a dict with: "tot_en", "osc_stre", "rot_stre", "dominant_contributions" (list of contributions),
            "moments_columns" (list of dicts value of electric and magnetic moment for each column).
        """
        # See test files for example. Splits the sections for each irrep delimited by
        #  ==============================================================================
        #
        #                               I R R E P   a1
        #
        #  ==============================================================================
        r = r"={3,}\s+I R R E P\s+["+ irrep_re_group +r"]+\s+={3,}.*?(?=={3,}|SUMMARY|-{5,})"

        match_irreps = re.findall(r, self.string, re.DOTALL)

        if not match_irreps:
            return None

        # Compile the regex beforehand since they will be executed several times.

        # Splits the excitations for each irrep. For example a group starting with:
        #                          1 singlet a1 excitation
        exc_head_group = r"^\s+\d+\s+(?:singlet\s+|triplet\s+|)[" + irrep_re_group + r"]+?\s+excitation"
        r_excited = exc_head_group + r".*?(?=" + exc_head_group + r"|\Z)"
        regex_excited = re.compile(r_excited, re.DOTALL | re.MULTILINE)

        # Finds the value of the total energy
        r_tot_en = r"Total energy:\s+(" + float_number_all_re + r")"
        regex_tot_en = re.compile(r_tot_en)

        # Parses:
        #  Oscillator strength:
        #     velocity representation:             0.1689055727340166
        #     length representation:               0.7845154315911747E-01
        #     mixed representation:                0.1151125659046751
        r_osc = r"Oscillator strength.*?length representation:\s+(" + float_number_all_re + r")"
        regex_osc = re.compile(r_osc, re.DOTALL)

        # Parses:
        #  Rotatory strength:
        #     velocity representation:              0.000000000000000
        #     velocity rep. / 10^(-40)erg*cm^3:     0.000000000000000
        #     length representation:                0.000000000000000
        #     length rep. / 10^(-40)erg*cm^3:       0.000000000000000
        r_rot = r"Rotatory strength.*?length representation:\s+(" + float_number_all_re + r")"
        regex_rot = re.compile(r_rot, re.DOTALL)

        # Parses the dominant contributions. Example:
        #  Dominant contributions:
        #       occ. orbital   energy / eV   virt. orbital     energy / eV   |coeff.|^2*100
        #         3 a1             -8.44           4 a1              0.92       99.2
        #
        # Different points to stop depending on singlet, triplet or uhf.
        r_dom_contrib = r"Dominant contributions:\s+.*?coeff\.\|\^2\*100\s+(.*?)(?=Change of electron number|<|S\^2|1st)"
        regex_dom_contrib = re.compile(r_dom_contrib, re.DOTALL | re.MULTILINE)

        # The following parse the section with the electric and magnetic moments
        r_elec_dip = r"Electric transition dipole moment \(length rep\.\):(.*?)Magnetic transition"
        regec_elec_dip = re.compile(r_elec_dip, re.DOTALL)

        r_mag_dip = r"Magnetic transition dipole moment / i:(.*?)Electric quadrupole"
        regec_mag_dip = re.compile(r_mag_dip, re.DOTALL)

        r_elec_quad = r"Electric quadrupole transition moment:(.*)"
        regec_elec_quad = re.compile(r_elec_quad, re.DOTALL)

        excited_data = {}

        for irrep_group in match_irreps:
            r = r"I R R E P\s+([" + irrep_re_group + r"]+)"
            irrep = re.search(r, irrep_group).group(1)

            match_excitations = regex_excited.findall(irrep_group)
            if not match_excitations:
                raise RuntimeError("Could not extract the single excitation contributions")

            single_excitations_list = []
            # loop over the excitations
            for exc_group in match_excitations:
                exc_data = dict(tot_en=None,
                                osc_stre=None,
                                rot_stre=None,
                                dominant_contributions=None,
                                moments_columns=None)

                match_tot_en = regex_tot_en.search(exc_group)
                if match_tot_en:
                    exc_data["tot_en"] = convert_float(match_tot_en.group(1))

                match_osc_stre = regex_osc.search(exc_group)
                if match_osc_stre:
                    exc_data["osc_stre"] = convert_float(match_osc_stre.group(1))

                match_rot_stre = regex_rot.search(exc_group)
                if match_rot_stre:
                    exc_data["rot_stre"] = convert_float(match_rot_stre.group(1))

                match_dominant = regex_dom_contrib.search(exc_group)
                if match_dominant:
                    dominant_contributions = []
                    # loop over dominant contributions to the excitation
                    for l in match_dominant.group(1).splitlines():
                        l = l.strip()
                        if not l:
                            continue
                        s = l.split()
                        # there could be 7 or 9 columns depending if it is a rhf or uhf, respectively.
                        if len(s) == 7:
                            occ_orb = dict(index=convert_int(s[0]), irrep=s[1], energy=convert_float(s[2]), spin=None)
                            virt_orb = dict(index=convert_int(s[3]), irrep=s[4], energy=convert_float(s[5]), spin=None)
                            dominant_contributions.append(dict(occ_orb=occ_orb, virt_orb=virt_orb,
                                                               coeff=convert_float(s[6])))
                        elif len(s) == 9:
                            occ_orb = dict(index=convert_int(s[0]), irrep=s[1], energy=convert_float(s[3]), spin=s[2])
                            virt_orb = dict(index=convert_int(s[4]), irrep=s[5], energy=convert_float(s[7]), spin=s[6])
                            dominant_contributions.append(dict(occ_orb=occ_orb, virt_orb=virt_orb,
                                                               coeff=convert_float(s[8])))
                        else:
                            raise RuntimeError("error while parsing line of escf dominant contribution: {}".format(l))

                    exc_data["dominant_contributions"] = dominant_contributions

                # Depending on the irrep there might be more "columns" with the dipole
                # and quadrupole values. Generate a list of dictionaries with these data.
                moments_columns = []
                if " column" in exc_group:
                    columns_to_parse = exc_group.split(" column")[1:]
                else:
                    columns_to_parse = [exc_group]

                # loop over the columns for the electric and magnetic moments.
                for str_c in columns_to_parse:
                    moments_data = {}
                    match_elec_dip = regec_elec_dip.search(str_c)
                    if match_elec_dip:
                        el_dip = []
                        for l in match_elec_dip.group(1).splitlines():
                            l = l.strip()
                            if l and l[0] in ("x", "y", "z"):
                                el_dip.append(convert_float(l.split()[1]))

                        moments_data["electric_dipole"] = el_dip
    
                    match_mag_dip = regec_mag_dip.search(str_c)
                    if match_mag_dip:
                        mag_dip = []
                        for l in match_mag_dip.group(1).splitlines():
                            l = l.strip()
                            if l and l[0] in ("x", "y", "z"):
                                mag_dip.append(convert_float(l.split()[1]))

                        moments_data["magnetic_dipole"] = mag_dip
    
                    match_elec_quad = regec_elec_quad.search(str_c)
                    if match_elec_quad:
                        quadrupole = np.zeros((3, 3))
                        d_quadrupole = {}
                        for l in match_elec_quad.group(1).splitlines():
                            s = l.split()
                            if " xx " in l:
                                quadrupole[0, 0] = convert_float(s[1])
                            elif " yy " in l:
                                quadrupole[1, 1] = convert_float(s[1])
                                d_quadrupole["trace"] = 3 * convert_float(s[-1])
                            elif " zz " in l:
                                quadrupole[2, 2] = convert_float(s[1])
                            elif " xy " in l:
                                quadrupole[0, 1] = quadrupole[1, 0] = convert_float(s[1])
                            elif " xz " in l:
                                quadrupole[0, 2] = quadrupole[2, 0] = convert_float(s[1])
                                d_quadrupole["anisotropy"] = convert_float(s[-1])
                            elif " yz " in l:
                                quadrupole[1, 2] = quadrupole[2, 1] = convert_float(s[1])
    
                        d_quadrupole["moment"] = quadrupole.tolist()
                        moments_data["electric_quadrupole"] = d_quadrupole

                    if moments_data:
                        moments_columns.append(moments_data)
                    
                if moments_columns:
                    exc_data["moments_columns"] = moments_columns

                single_excitations_list.append(exc_data)

            excited_data[irrep] = single_excitations_list

        return excited_data

    @lazy_property
    def rdgrad_memory(self):
        """
        Memory allocated for rdgrad.
        Valid for rdgrad.

        Returns:
            int.
        """
        # Parses: Memory allocated for RDGRAD: 1 MiB
        r = r"Memory allocated for RDGRAD.*$"

        match = re.search(r, self.string, re.MULTILINE)

        if not match:
            return None
        else:
            return convert_int(match.group().split()[-2])

    @lazy_property
    def gradient(self):
        """
        Values of the gradients and dielectric data.
        Valid for gradient executables.

        Returns:
            dict with "gradients", a matrix of shape (natoms, 3) containing the values of
                the cartesian gradients, and "dielectric".
        """
        # Parse the section starting with:
        #           ------------------------------------------------
        #            cartesian gradient of the energy (hartree/bohr)
        #           ------------------------------------------------
        r = "cartesian gradient of the energy.*cartesian gradients written onto"
        match = re.search(r, self.string, re.DOTALL)

        if not match:
            return None

        grad_str = match.group()

        grad_data = {}

        # Parse lines of the form:
        # dE/dx  0.0000000D+00 -0.2705802D-01  0.2705802D-01
        # dE/dy  0.0000000D+00  0.0000000D+00  0.0000000D+00
        # dE/dz  0.4144103D-01 -0.2070390D-01 -0.2070390D-01
        #
        # Notice that they can span multiple sets of lines.

        gradients = [[], [], []]
        for i, d in enumerate(["x", "y", "z"]):
            match = re.findall(r"dE/d"+d+r"(.*?)$", grad_str, re.MULTILINE)
            for g_vals in match:
                g = [convert_float(f) for f in g_vals.split()]
                gradients[i].extend(g)

        grad_data["gradients"] = np.transpose(gradients).tolist()

        # Parses the lines:
        #  exx =      -0.076690 eyy =       0.000000 ezz =      -0.045091
        #  eyz =       0.000000 exz =      -0.000000 exy =       0.000000

        dielectric = np.zeros((3,3))

        for l in grad_str.splitlines():
            if "exx" in l:
                s = l.split()
                dielectric[0,0] = convert_float(s[2])
                dielectric[1,1] = convert_float(s[5])
                dielectric[2,2] = convert_float(s[8])
            elif "eyz" in l:
                s = l.split()
                dielectric[1,2] = dielectric[2,1] = convert_float(s[2])
                dielectric[0,2] = dielectric[2,0] = convert_float(s[5])
                dielectric[0,1] = dielectric[1,0] = convert_float(s[8])

        grad_data["dielectric"] = dielectric.tolist()

        return grad_data

    @lazy_property
    def egrad_excited_state(self):
        """
        Information about the excited state chosen for optimization in egrad.
        Valid only for egrad.

        Returns:
            dict with "index" indicating the (1-based) index of the state used for
            optimization and "default" to True if the default values used for $exopt.
        """
        # Parses:  Excited state no.    3 chosen for optimization
        r = r"Excited state no.\s+(\d+)\s+chosen for optimization"
        match = re.search(r, self.string)

        if not match:
            return None

        data = dict(index=convert_int(match.group(1)))
        data["default"] = "Data group $exopt missing or empty" in self.string

        return data

    @lazy_property
    def statpt_info(self):
        """
        Initial information provided in statpt.
        Valid only for statpt.

        Returns:
            dict with "max_trust_radius", "min_trust_radius", "init_trust_radius", "min_grad_norm_for_gdiis",
            "prev_steps_for_gdiis", "hessian_update_method", "thr_energy_change", "thr_max_displ", "thr_max_grad",
            "thr_rms_displ", "thr_rms_grad", "use_defaults", "final_step_radius".
        """
        # See test files for examples. Parses the block starting with:
        # ***************  Stationary point options ******************
        r = r"Stationary point options[\s\*]{60,}(.+?)\*{60,}"
        match = re.search(r, self.string, re.DOTALL)

        data = dict(max_trust_radius=None,
                    min_trust_radius=None,
                    init_trust_radius=None,
                    min_grad_norm_for_gdiis=None,
                    prev_steps_for_gdiis=None,
                    hessian_update_method=None,
                    thr_energy_change=None,
                    thr_max_displ=None,
                    thr_max_grad=None,
                    thr_rms_displ=None,
                    thr_rms_grad=None,
                    use_defaults=False,
                    final_step_radius=None)

        if match:
            for l in match.group(1).splitlines():
                if "Maximum allowed trust radius" in l:
                    data["max_trust_radius"] = convert_float(l.split()[-1])
                if "Minimum allowed trust radius" in l:
                    data["min_trust_radius"] = convert_float(l.split()[-1])
                if "Initial trust radius" in l:
                    data["init_trust_radius"] = convert_float(l.split()[-1])
                if "GDIIS used if gradient norm" in l:
                    data["min_grad_norm_for_gdiis"] = convert_float(l.split()[-1])
                if "Number of previous steps for GDIIS" in l:
                    data["prev_steps_for_gdiis"] = convert_int(l.split()[-1])
                if "Hessian update method" in l:
                    data["hessian_update_method"] = l.split()[-1]
                if "Threshold for energy change" in l:
                    data["thr_energy_change"] = convert_float(l.split()[-1])
                if "Threshold for max displacement element" in l:
                    data["thr_max_displ"] = convert_float(l.split()[-1])
                if "Threshold for max gradient element" in l:
                    data["thr_max_grad"] = convert_float(l.split()[-1])
                if "Threshold for RMS of displacement" in l:
                    data["thr_rms_displ"] = convert_float(l.split()[-1])
                if "Threshold for RMS of gradient" in l:
                    data["thr_rms_grad"] = convert_float(l.split()[-1])

        data["use_defaults"] = "Keyword $statpt not found - using default options" in self.string

        r = r"$\s+Final\s+step\s+radius\s*:\s+(" + float_number_e_re + r")"
        match = re.search(r, self.string, re.MULTILINE)
        if match:
            data["final_step_radius"] = convert_float(match.group(1))

        if all(not v for v in data.values()):
            return None

        return data

    @lazy_property
    def relax_info(self):
        """
        Initial information provided in relax.
        Valid only for relax.

        Returns:
            dict with "optimizations" list of optimized quantities and "thr_int_coord".
        """
        data = dict(optimizations=None,
                    thr_int_coord=None)

        # Parses for example:
        #  optimization will be performed with respect to
        #  - INTERNAL coordinates
        r = r"optimization will be performed with respect to\s+(.*?)-{10,}"
        match = re.search(r, self.string, re.DOTALL)
        if match:
            opt = []
            for l in match.group(1).splitlines():
                l = l.strip()
                if l.startswith("-"):
                    opt.append(l.replace("-", "", 1).strip())
            data["optimizations"] = opt

        r = r"$\s+convergence\s+criterion\s+for\s+internal\s+coordinates\s*:\s+(" + float_number_e_re + r")"
        match = re.search(r, self.string, re.MULTILINE)
        if match:
            data["thr_int_coord"] = convert_float(match.group(1))

        if all(v is None for v in data.values()):
            return None

        return data

    @lazy_property
    def relax_gradient_values(self):
        """
        Gradient values extracted from the relax/stapt output that taking into account
        the frozen coordinates.
        Valid for relax and statpt.

        Returns:
            dict with "norm_cartesian", "norm_internal", "max_internal"
        """
        #Parses:
        #   norm of actual CARTESIAN gradient : .54369
        #   norm of actual  INTERNAL gradient : .55998
        #   maximum norm of actual  INTERNAL gradient : .46843

        data = dict(norm_cartesian=None,
                    norm_internal=None,
                    max_internal=None)

        r = r"$\s+norm\s+of\s+actual\s+CARTESIAN\s+gradient\s*:\s+("+float_number_e_re+r")"
        match = re.search(r, self.string, re.MULTILINE)
        if match:
            data["norm_cartesian"] = convert_float(match.group(1))

        r = r"$\s+norm\s+of\s+actual\s+INTERNAL\s+gradient\s*:\s+("+float_number_e_re+r")"
        match = re.search(r, self.string, re.MULTILINE)
        if match:
            data["norm_internal"] = convert_float(match.group(1))

        r = r"maximum\s+norm\s+of\s+actual\s+INTERNAL\s+gradient\s*:\s+("+float_number_e_re+r")"
        match = re.search(r, self.string)
        if match:
            data["max_internal"] = convert_float(match.group(1))

        if all(v is None for v in data.values()):
            return None

        return data

    @lazy_property
    def relax_conv_info(self):
        """
        Final information about convergence according to what is defined in the
        control file.
        Valid for relax and statpt.

        Returns:
            dict with "energy_change", "rms_displ", "rms_grad", "max_displ", "max_grad".
            For each one a dict with "value", "thr" and "conv".
        """
        # Parses:
        #       ******************************************************************
        #                           CONVERGENCE INFORMATION
        #
        #                                Converged?     Value      Criterion
        #              Energy change         no      76.3430059   0.0000010
        #              RMS of displacement   no       0.0517230   0.0005000
        #              RMS of gradient       no       0.0340771   0.0005000
        #              MAX displacement      no       0.0731468   0.0010000
        #              MAX gradient          no       0.0481922   0.0010000
        #       ******************************************************************
        r = r"CONVERGENCE INFORMATION.*?Criterion\s*(.+?)\*{60,}"
        match = re.search(r, self.string, re.DOTALL)

        if not match:
            return None

        def get_single_data(line):
            s = line.split()
            converged = "yes" == s[-3]
            return {"value": convert_float(s[-2]), "thr": convert_float(s[-1]), "conv": converged}

        data = dict(energy_change=None,
                    rms_displ=None,
                    rms_grad=None,
                    max_displ=None,
                    max_grad=None)
        for l in match.group(1).splitlines():
            if "Energy change" in l:
                data["energy_change"] = get_single_data(l)
            if "RMS of displacement" in l:
                data["rms_displ"] = get_single_data(l)
            if "RMS of gradient" in l:
                data["rms_grad"] = get_single_data(l)
            if "MAX displacement" in l:
                data["max_displ"] = get_single_data(l)
            if "MAX geom. grad." in l or "MAX gradient" in l:
                data["max_grad"] = get_single_data(l)

        return data

    @lazy_property
    def aoforce_numerical_integration(self):
        """
        Information about the numerical integration in aoforce.
        Will be present only if a proper aoforce is run. Absent if running after numforce.

        Returns:
            dict with "memory" (a dict with "core_memory_dft", "memory_per_atom", "atoms_per_loop") and
            "construction_timings" (a list of lists with construction timings. For each element:
                [str with description of the timing, cpu time in seconds, wall time in seconds]).
        """

        # This section may not be present if aoforce runs after a numforce but also
        #  if the calculation is not dft.
        r = r"PREPARING NUMERICAL INTEGRATION.*atoms per loop"
        match = re.search(r, self.string, re.DOTALL)

        if match:
            memory = {}
            for l in match.group().splitlines():
                if "Remaining core memory for DFT" in l:
                    memory["core_memory_dft"] = convert_int(l.split()[-2])
                if "Memory needed per atom" in l:
                    memory["memory_per_atom"] = convert_int(l.split()[-2])
                if "atoms per loop" in l:
                    memory["atoms_per_loop"] = convert_int(l.split()[-4])
        else:
            memory = None

        timings = []
        # this matches a series a of "CONSTRUCTING xxx to get their timings. Example:
        #  CONSTRUCTING first deriv. of <mu|x,y,z|nu> -> Dip. deriv.
        #   ...
        #       ...terminated. cpu:       0.00       wall:       0.00
        # it skips the CONSTRUCTING integral bounds, that does not contain any timing.
        r = r"CONSTRUCTING((?! integral bounds).*?wall:\s+\d\.\d+)\s+(?=CONSTRUCTING|SOLVING|-------)"
        for m in re.findall(r, self.string, re.DOTALL):
            # consider a different case for the second deriv. of 2e energy
            # since if it a dft calculations there are two timings there.
            if "second deriv. of 2e energy" in m:
                gl = m.splitlines()
                treating_line = ""
                for i in range(1, len(gl)):
                    if "treating" in gl[i]:
                        treating_line = gl[i]
                    if "...terminated" in gl[i]:
                        gls = gl[i].split()
                        timings.append((gl[0] + "\n" + treating_line, convert_float(gls[-3]), convert_float(gls[-1])))

            else:
                lines = []
                cpu = None
                wall = None
                for l in m.splitlines():
                    if "...terminated" in l:
                        s = l.split()
                        cpu = convert_float(s[-3])
                        wall = convert_float(s[-1])
                    elif "->" in l or ("<" and ">" in l) or ("[" and "]" in l):
                        lines.append(l)
                timings.append(("\n".join(lines), cpu, wall))

        if not memory and not timings:
            return None

        data = dict(memory=memory, construction_timings=timings)

        return data

    @lazy_property
    def aoforce_analysis(self):
        """
        Output of rotational and vibrational analysis in aoforce.

        Returns:
            dict with "rotational" (a dict with "dipole_moment", "b", "intensities", "m") and "vibrational".
            Vibrational, contains dict with keys "frequencies", "symmetries", "ir", "dDIP_dQ", "intensities",
            "intensities_perc", "raman", "eigenvectors", "reduced_masses". One value for each frequency.
        """
        # See test files for examples. Parses the section starting with:
        #              -----------------------------------
        #              rotational and vibrational analysis
        #              -----------------------------------
        split1 = self.string.split("dipole moment in principle axis system")
        if len(split1) != 2:
            return None

        split2 = split1[1].split("NORMAL MODES and VIBRATIONAL FREQUENCIES")
        if len(split2) != 2:
            return None

        # Rotational section
        # Parses:
        #  dipole moment in principle axis system (a.u.) :
        #     -0.0000000000     0.0000000000     0.8436248828
        #  norm :  0.843624882772472
        #
        #  rotational constants b for rotations around axis of inertia
        #  and optical intensities for (1 <-- 0) transition
        #
        #    b   :    25.5620992978     9.0755697626    14.0715251781   (cm**(-1))
        #    b   :      766332.4580      272078.7367      421853.7121     (MHz)
        #   int. :     0.0000000000     0.0000000000     0.7117029428     (a.u.)
        #
        #    x   :     1.0000000000     0.0000000000     0.0000000000
        #    y   :     0.0000000000     1.0000000000     0.0000000000
        #    z   :     0.0000000000     0.0000000000     1.0000000000


        rot = {}
        rot_lines = split2[0].splitlines()

        rot["dipole_moment"] = [convert_float(f) for f in rot_lines[1].split()]
        b = None
        rot_intensities = None
        m = None
        for l in rot_lines:
            l = l.strip()
            if l.startswith("b ") and "cm" in l:
                s = l.split()
                b = [convert_float(s[i]) for i in (2,3,4)]
            if l.startswith("int."):
                s = l.split()
                rot_intensities = [convert_float(s[i]) for i in (2,3,4)]
            if l.startswith("x "):
                # the x :   y :   z : section may be missing.
                # Initialize m here only if it is actually filled
                m = [[0., 0., 0.]] * 3
                s = l.split()
                m[0] = [convert_float(s[i]) for i in (2,3,4)]
            if l.startswith("y "):
                s = l.split()
                m[1] = [convert_float(s[i]) for i in (2,3,4)]
            if l.startswith("z "):
                s = l.split()
                m[2] = [convert_float(s[i]) for i in (2,3,4)]

        #TODO verify the meaning of int and m
        rot["b"] = b
        rot["intensities"] = rot_intensities
        rot["m"] = m

        # Vibrational section
        # Parses groups of the form:
        #        mode               7        8        9
        #
        #      frequency        1603.56  3527.41  3641.97
        #
        #      symmetry            a1       a1       b1
        #
        #         IR               YES      YES      YES
        # |dDIP/dQ|   (a.u.)     0.0059   0.0002   0.0031
        # intensity (km/mol)      62.07     0.11    16.54
        # intensity (  %   )     100.00     0.17    26.64
        #
        #        RAMAN             YES      YES      YES
        #
        #   1   o           x  -0.00000  0.00000  0.06983
        #                   y   0.00000  0.00000  0.00000
        #                   z   0.07060  0.04999 -0.00000
        #   2   h           x   0.42845 -0.58424 -0.55423
        #                   y   0.00000  0.00000  0.00000
        #                   z  -0.56031 -0.39676 -0.43634
        #   3   h           x  -0.42845  0.58424 -0.55423
        #                   y  -0.00000  0.00000  0.00000
        #                   z  -0.56031 -0.39676  0.43634
        #
        # reduced mass(g/mol)     1.083    1.045    1.081

        frequencies = []
        symmetries = []
        ir = []
        dDIP_dQ = []
        vib_intensities = []
        vib_intensities_perc = []
        raman = []
        eigenvectors = []
        reduced_masses = []

        regex_eigenvectors = re.compile(r"(.*?)(^[ \t]+\d+\s+[a-z]{1,2}.*?)(reduced mass.*)", re.DOTALL|re.MULTILINE)

        def get_active_value(v):
            """
            Converts the value to a bool. If "-" is set to None.
            """
            if v == "YES":
                return True
            elif v == "NO":
                return False
            elif v == "-":
                return None
            else:
                raise RuntimeError("cannot correctly parse active value {}".format(v))

        # splits the block of lines with one set of frequencies and loops over them
        r = r"(mode\s{5}.*?)(?=mode |\*{10,})"
        for block in re.findall(r, split2[1], re.DOTALL):
            match = regex_eigenvectors.search(block)
            if not match:
                continue
            block_freqs = match.group(1) + match.group(3)
            block_eigenvectors = match.group(2).strip()
            freqs_lines = block_freqs.splitlines()
            n_freqs = len(freqs_lines[0].split()) - 1

            # loop over the lines of the block, extract all the properties
            for l in freqs_lines:
                s = l.split()
                if "frequency" in l:
                    # put imaginary frequencies to negative values.
                    frequencies.extend(convert_float(f.replace("i", "-")) for f in s[1:])
                elif "symmetry" in l:
                    # put the symmetries that are present
                    symmetries.extend(sym for sym in s[1:])
                    # add None values for those that are not present
                    symmetries.extend(None for i in range(n_freqs-len(s)+1))
                elif "IR " in l:
                    ir.extend(get_active_value(v) for v in s[1:])
                elif "|dDIP/dQ|" in l:
                    dDIP_dQ.extend(convert_float(f) for f in s[2:])
                elif "intensity (km/mol)" in l:
                    vib_intensities.extend(convert_float(f) for f in s[2:])
                elif "intensity" in l and "%" in l:
                    vib_intensities_perc.extend(convert_float(f) for f in s[4:])
                elif "RAMAN " in l:
                    raman.extend(get_active_value(v) for v in s[1:])
                elif "reduced mass" in l:
                    reduced_masses.extend(convert_float(f) for f in s[2:])

            eig_lines = block_eigenvectors.splitlines()
            if len(eig_lines) % 3 != 0:
                raise RuntimeError("number of line to be parsed for eigenvectors seems wrong: {}".format(len(eig_lines)))

            n_atoms = len(eig_lines) // 3
            eigv_block = np.zeros((n_freqs, n_atoms, 3))
            for i in range(n_atoms):
                for j in range(3):
                    s = eig_lines[3*i+j].split()
                    eigv_block[:, i, j] = [convert_float(f) for f in s[-n_freqs:]]

            eigenvectors.extend(eigv_block.tolist())

        energies = None
        re_zpve = r"zero point VIBRATIONAL energy.*?\*{10,}"
        match = re.search(re_zpve, self.string, re.DOTALL)
        if match:
            energies = {"zpve": None,
                        "scf": None,
                        "total": None}
            for l in match.group().splitlines():
                if "zero point VIBRATIONAL energy" in l:
                    energies["zpve"] = convert_float(l.split()[-3])
                if "SCF-energy" in l:
                    energies["scf"] = convert_float(l.split()[-2])
                if "SCF + E(vib0)" in l:
                    energies["total"] = convert_float(l.split()[-2])

        vib = dict(frequencies=frequencies, symmetries=symmetries, ir=ir, dDIP_dQ=dDIP_dQ,
                   intensities=vib_intensities, intensities_perc=vib_intensities_perc,
                   raman=raman, eigenvectors=eigenvectors, reduced_masses=reduced_masses,
                   energies=energies)

        tot_n_freqs = len(frequencies)
        if any(k != "energies" and len(v) != tot_n_freqs for k, v in vib.items()):
            raise RuntimeError("Error parsing the vibrational frequencies, for some quantity the list of "
                               "values does not contain the correct number of items")

        data = dict(rotational=rot, vibrational=vib)

        return data

    @lazy_property
    def mp2_data(self):
        """
        MP2 data: information of MP2 calculation.

        Returns:
            dict with "energy_only".
        """

        # TODO: implement or mpgrad and rimp2/ricc2

        return None

    @lazy_property
    def mp2_results(self):
        """
        MP2 results.

        Returns:
            dict with "energy".
        """

        energy = None
        # Try to parse from rimp2/ricc2 programs
        r = r"\*{62,62}\s+\*\s+\*\s+\*<{10,10}\s+"
        r += r"GROUND\s+STATE\s+FIRST-ORDER\s+PROPERTIES"
        r += r"\s+>{11,11}\*\s+\*\s+\*\s+\*{62,62}\s+"
        r += r"-{48,48}\s+Method\s+:\s+MP2\s+Total\s+Energy\s+:\s+("
        r += float_number_all_re
        r += r")\s+"
        m = re.findall(r, string=self.string)
        if len(m) == 1:
            energy = convert_float(m[0])
        elif len(m) > 1:
            raise RuntimeError("Error parsing the MP2 results. Multiple occurrences of MP2 results found.")

        # Try to parse from mpgrad program
        r = r"\*{53,53}\s+\*\s+\*\s+\*\s+"
        r += r"SCF-energy\s+:\s+(" + float_number_all_re + r")\s+\*\s+\*\s+"
        r += r"MP2-energy\s+:\s+(" + float_number_all_re + r")\s+\*\s+\*\s+"
        r += r"total\s+:\s+(" + float_number_all_re + r")\s+\*\s+\*\s+"
        r += r"\*\s+\*\s+\(MP2-energy\s+evaluated\s+from\s+T2\s+amplitudes\)\s+\*\s+\*\s+\*\s+\*{53,53}"
        m = re.findall(r, string=self.string)
        if len(m) == 1:
            if energy is not None:
                raise RuntimeError("Error parsing the MP2 results. "
                                   "Found MP2 results from both mpgrad- and rimp2/ricc2-like calculations.")
            energy = convert_float(m[0][2])
        elif len(m) > 1:
            raise RuntimeError("Error parsing the MP2 results. Multiple occurrences of MP2 results found.")

        return dict(energy=energy)

    def get_split_jobex_parsers(self):
        """
        Given the string of the "job.last" output file of jobex generates the
        parsers for each of the step: energy, gradient and relax.

        Returns:
            namedtuple:
                exec_en: executable of the energy step.
                parser_en: parser with the output of the energy step.
                exec_grad: executable of the gradient step.
                parser_grad: parser with the output of the gradient step.
                exec_relax: executable of the energy step.
                parser_relax: parser with the output of the relax step.
        """
        tot_split = self.string.split("next step =")
        p_en = p_grad = p_relax = None
        exec_en = exec_grad = exec_relax = None

        for s in tot_split:
            p = Parser(s)
            header = p.header
            if not header:
                continue

            executable = header["executable"]

            if executable in ("dscf", "ridft"):
                p_en = p
                exec_en = executable
            elif executable in ("grad", "rdgrad", "egrad", "mpgrad", "ricc2"):
                p_grad = p
                exec_grad = executable
            elif executable in ("relax", "statpt", "frog"):
                p_relax = p
                exec_relax = executable

        JobexParsers = namedtuple("JobexParsers", ["exec_en", "parser_en", "exec_grad", "parser_grad", "exec_relax", "parser_relax"])

        return JobexParsers(exec_en, p_en, exec_grad, p_grad, exec_relax, p_relax)

    def grep_line(self, text, nlines=0):
        """
        Searches the specified text in a single line and return the full line plus
        the "nlines" following the matching one. Returns the first match.

        Args:
            text (str): the text to be searched.
            nlines (int): the number of lines

        Returns:
            str: the text found
        """

        r = r"^.*" + text + r".*$(\n.*){0," + str(nlines) + r"}"

        m = re.search(r, self.string, re.MULTILINE)

        if not m:
            return None

        return m.group(0)

    def get_value(self, text, chunk, occurrence=0, converter=None):
        """
        Finds the i-th occurrence of a line containing the specified text,
        splits the line and takes the specified chunk obtained from the split.
        Optionally converts it to a specific format.

        Args:
            text (str): text to search
            chunk (int): the chunk to take after the split of the line.
                Can be negative as it identifies the index in the list.
            occurrence (int): the 0-based index of the occurrence of the "text"
                in the string. Can be negative (e.g. -1 is the last occurrence).
            converter (callable): called to convert the selected chunk to some
                other type. It can be simply int, but for a float it would be better
                to use the convert_float function instead.

        Returns:
            the converted selected chunk.
        """

        r = r"^.*" + text + r".*$"

        m = re.findall(r, self.string, re.MULTILINE)

        if not m:
            return None

        line = m[occurrence]

        v = line.split()[chunk]

        if converter:
            v = converter(v)

        return v
