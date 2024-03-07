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
from turbomoleio.output.files import EscfOutput, ScfOutput
from turbomoleio.testing import run_itest

# structures = ['aceton', 'ch4', 'h2o', 'h3cbr', 'methanol',
#               'nh3', 'phenol', 'sf4', 'sih4']
structures = ["h2o", "nh3"]


@pytest.mark.integration
class TestEscf:
    @pytest.mark.parametrize("structure", structures)
    def test_run_ridft_escf(self, structure):
        dp = get_define_template("ridft_escf")
        dp["desy"] = True
        dp["ex_all_states"] = 12

        assert run_itest(
            ["ridft", "escf"],
            dp,
            structure,
            "ridft_escf_{}_std".format(structure),
            [ScfOutput, EscfOutput],
        )

    @pytest.mark.parametrize("structure", structures)
    def test_run_dscf_escf(self, structure):
        dp = get_define_template("dscf_escf")
        dp["desy"] = False
        dp["ex_all_states"] = None
        dp["ex_irrep_states"] = {"a": 8}

        assert run_itest(
            ["dscf", "escf"],
            dp,
            structure,
            "dscf_escf_{}_std".format(structure),
            [ScfOutput, EscfOutput],
        )

    @pytest.mark.parametrize("structure", structures)
    def test_run_dscf_escf_triplet(self, structure):
        dp = get_define_template("dscf_escf")
        dp["desy"] = True
        dp["ex_all_states"] = 10
        dp["ex_multi"] = "triplet"

        assert run_itest(
            ["dscf", "escf"],
            dp,
            structure,
            "dscf_escf_{}_triplet".format(structure),
            [ScfOutput, EscfOutput],
        )

    @pytest.mark.parametrize("structure", structures)
    def test_run_dscf_escf_charged(self, structure):
        dp = get_define_template("dscf_escf")
        dp["desy"] = True
        dp["ex_all_states"] = 10
        dp["charge"] = 1

        assert run_itest(
            ["dscf", "escf"],
            dp,
            structure,
            "dscf_escf_{}_charged".format(structure),
            [ScfOutput, EscfOutput],
        )

    def test_run_dscf_escf_set_dict(self):
        """
        Tests some properties that are set with a dictionary.

        These tests are structure specific.
        It also sets some specific parameters to add coverage.
        The structure is nh3_wrong_internal, which is the same as nh3, but with
        internal coordinates not matching with the cartesian, to check that the
        code will handle that as well.
        """

        define_opt = get_define_template("dscf_escf")
        define_opt["metric"] = 2
        define_opt["sym"] = None
        define_opt["desy"] = True
        define_opt["desy_eps"] = 1e-5
        define_opt["basis"] = None
        define_opt["basis_atom"] = {"n": "def2-SV(P)", "h": "def2-SV(P)"}
        define_opt["ex_all_states"] = None
        define_opt["ex_irrep_states"] = {"a1": 10, "a2": 8}

        assert run_itest(
            ["dscf", "escf"],
            define_opt,
            "nh3_wrong_internal",
            "dscf_escf_nh3_set_dict",
            ScfOutput,
        )
