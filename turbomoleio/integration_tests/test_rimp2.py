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
from turbomoleio.output.files import Ricc2Output, ScfOutput
from turbomoleio.testfiles.utils import run_itest

# structures = ['aceton', 'ch4', 'h2o', 'h3cbr', 'methanol',
#               'nh3', 'phenol', 'sf4', 'sih4']
structures = ["h2o", "nh3"]


@pytest.mark.integration
class TestRimp2:
    @pytest.mark.parametrize("structure", structures)
    def test_run_ridft_rimp2(self, structure):

        assert run_itest(
            ["ridft", "rimp2"],
            get_define_template("ridft_rimp2"),
            structure,
            "ridft_rimp2_{}_std".format(structure),
            [ScfOutput, Ricc2Output],
        )

    def test_run_ridft_adc2(self):

        define_opt = get_define_template("ridft_rimp2")
        define_opt["method"] = "adc(2)"
        define_opt["ex_mp2"] = {"a1": [1, 10], "a2": [3, 10]}
        define_opt["maxiter"] = 100
        define_opt["maxcor"] = 400

        assert run_itest(
            ["ridft", "rimp2"],
            define_opt,
            "nh3",
            "ridft_rimp2_nh3_adc2",
            [ScfOutput, None],
        )

    # FIXME according to TM: "CCSD(T) is no longer available in ricc2 : use ccsdf12"
    # move it to another module?
    # @pytest.mark.parametrize("structure", structures)
    # def test_run_ridft_ccsdt(self, structure):
    #
    #     define_opt = get_define_template("ridft_rimp2")
    #     define_opt["method"] = "ccsdt"
    #     define_opt["mp2energy"] = True
    #
    #     assert run_itest(["ridft", "rimp2"], define_opt, structure,
    #                      "ridft_rimp2_{}_ccsdt".format(structure), [ScfOutput, None])

    @pytest.mark.parametrize("structure", structures)
    def test_run_ridft_f12(self, structure):

        define_opt = get_define_template("ridft_rimp2")
        define_opt["use_f12"] = True
        # define_opt["use_f12*"] = True

        assert run_itest(
            ["ridft", "ricc2"],
            define_opt,
            structure,
            "ridft_ricc2_{}_f12".format(structure),
            [ScfOutput, Ricc2Output],
            datagroups_options={
                "ricc2": "\n   mp2 energy only",
                "rir12": "\n   comaprox T+V",
            },
        )

    @pytest.mark.parametrize("structure", structures)
    def test_run_ridft_f12x(self, structure):

        define_opt = get_define_template("ridft_rimp2")
        define_opt["use_f12"] = False
        define_opt["use_f12*"] = True

        assert run_itest(
            executables=["ridft", "rimp2"],
            define_options=define_opt,
            coord_filename=structure,
            control_reference_filename="ridft_rimp2_{}_f12*".format(structure),
            file_classes=[ScfOutput, Ricc2Output],
        )


# executables, define_options, coord_filename, control_reference_filename, file_classes
