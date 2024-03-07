# -*- coding: utf-8 -*-
# The turbomoleio package, a python interface to Turbomole
# for preparing inputs, parsing outputs and other related tools.
#
# Copyright (C) 2018-2022 BASF SE, Matgenix SRL.
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

import pytest

from turbomoleio.input.utils import get_define_template
from turbomoleio.output.files import ScfOutput
from turbomoleio.testing import run_itest

# structures = ['aceton', 'ch4', 'h2o', 'h3cbr', 'methanol',
#               'nh3', 'phenol', 'sf4', 'sih4']
structures = ["h2o", "nh3"]

# disp_corrections = ["DFT-D1", "DFT-D2", "DFT-D3", "DFT-D3 BJ"]
# DFT-D1 is not working anymore in Turbomole 7.5! Turbomole developers
# have been made aware of this fact. Even though DFT-D1 is not used
# anymore by many people, it is always good to keep it to be able
# to reproduce previous results. An issue has been opened on
# the turbomoleio to keep track of this problem. If DFT-D1 dispersion
# corrections are fixed in a subsequent Turbomole version (e.g. 7.6),
# Enable again DFT-D1.
disp_corrections = ["DFT-D2", "DFT-D3", "DFT-D3 BJ"]


@pytest.mark.integration
class TestDscf:
    @pytest.mark.parametrize("structure", structures)
    def test_run_dscf(self, structure):

        assert run_itest(
            "dscf",
            get_define_template("dscf"),
            structure,
            "dscf_{}_std".format(structure),
            ScfOutput,
        )

    @pytest.mark.parametrize("structure", structures)
    def test_run_dscf_ired(self, structure):
        define_opt = get_define_template("dscf")
        define_opt["ired"] = True

        assert run_itest(
            "dscf", define_opt, structure, "dscf_{}_ired".format(structure), ScfOutput
        )

    @pytest.mark.parametrize("structure", structures)
    def test_run_dscf_hf(self, structure):

        define_opt = get_define_template("dscf")
        define_opt["method"] = "hf"
        define_opt["desy"] = False
        define_opt["sym"] = "c1"
        define_opt["sym_eps"] = 0.001
        define_opt["scfiterlimit"] = 300
        define_opt["scfconv"] = 6
        define_opt["title"] = None

        assert run_itest(
            "dscf", define_opt, structure, "dscf_{}_hf".format(structure), ScfOutput
        )

    def test_run_dscf_usemo(self, test_data):
        """
        Tests the usemo functionalities with mos file
        """

        define_opt = get_define_template("dscf")

        define_opt["basis"] = "def2-SV(P)"
        define_opt["usemo"] = test_data / "mo" / "mos_nh3_nosym"

        assert run_itest("dscf", define_opt, "nh3", "dscf_nh3_usemo", ScfOutput)

    def test_run_dscf_ecp_atom(self):
        """
        Tests the usemo functionalities with mos file
        """

        define_opt = get_define_template("dscf")

        define_opt["basis"] = "def2-SV(P)"
        define_opt["ecp_atom"] = {"N": "ecp-2-sdf"}

        assert run_itest("dscf", define_opt, "nh3", "dscf_nh3_ecp_atom", ScfOutput)

    # TODO check if alpha/beta for usemo are needed
    # def test_run_dscf_usemo_alpha(self):
    #     """
    #     Tests the usemo functionalities with alpha and
    #     beta files in the case of unpaired electrons
    #     """
    #
    #     define_opt = get_define_template("dscf")
    #
    #     define_opt["basis"] = "def2-SV(P)"
    #     define_opt["charge"] = 1
    #     define_opt["unpaired_electrons"] = 1
    #     define_opt["usemo"] = os.path.join(get_tfp(), "mo", "alpha_beta_nh3")
    #
    #     assert run_itest("dscf", define_opt, "nh3", "dscf_nh3_usemo_alpha", {})

    def test_run_dscf_copymo(self, test_data):
        """
        Tests the copy functionalities with mos file
        """

        define_opt = get_define_template("dscf")

        define_opt["basis"] = "def2-SV(P)"
        define_opt["copymo"] = test_data / "mo" / "mos_nh3"

        assert run_itest("dscf", define_opt, "nh3", "dscf_nh3_copymo_mos", ScfOutput)

    def test_run_dscf_copymo_mos_unpaired(self, test_data):
        """
        Tests the copy functionalities with mos file
        """

        define_opt = get_define_template("dscf")

        define_opt["basis"] = "def2-SV(P)"
        define_opt["charge"] = 1
        define_opt["unpaired_electrons"] = 1
        define_opt["copymo"] = test_data / "mo" / "mos_nh3"

        assert run_itest(
            "dscf", define_opt, "nh3", "dscf_nh3_copymo_mos_unpaired", ScfOutput
        )

    def test_run_dscf_copymo_alpha_unpaired(self, test_data):
        """
        Tests the copy functionalities with mos file
        """

        define_opt = get_define_template("dscf")

        define_opt["basis"] = "def2-SV(P)"
        define_opt["charge"] = 1
        define_opt["unpaired_electrons"] = 1
        define_opt["copymo"] = test_data / "mo" / "alpha_beta_nh3"

        assert run_itest(
            "dscf", define_opt, "nh3", "dscf_nh3_copymo_alpha_unpaired", ScfOutput
        )

    @pytest.mark.parametrize("structure", structures)
    def test_run_dscf_cosmo(self, structure):

        define_opt = get_define_template("dscf")
        define_opt["use_cosmo"] = True
        define_opt.update(
            {
                "epsilon": 60.0,
                "nppa": 92,
                "nspa": None,
                "disex": 0,
                "rsolv": 1.3,
                "routf": None,
                "cavity": "open",
            }
        )
        assert run_itest(
            "dscf", define_opt, structure, "dscf_{}_cosmo".format(structure), ScfOutput
        )

    @pytest.mark.parametrize("disp", disp_corrections)
    def test_run_dscf_disp(self, disp):

        define_opt = get_define_template("dscf")
        define_opt["disp"] = disp
        disp_fn = "dscf_nh3_disp_{}".format("_".join(disp.split()))
        assert run_itest("dscf", define_opt, "nh3", disp_fn, ScfOutput)
