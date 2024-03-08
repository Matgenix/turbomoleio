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


@pytest.mark.integration
class TestRidft:
    @pytest.mark.parametrize("structure", structures)
    def test_run_ridft(self, structure, test_data):

        assert run_itest(
            "ridft",
            get_define_template("ridft"),
            structure,
            "ridft_{}_std".format(structure),
            ScfOutput,
            test_data=test_data,
        )

    @pytest.mark.parametrize("structure", structures)
    def test_run_ridft_rijk(self, structure, test_data):
        define_opt = get_define_template("ridft")
        define_opt["ri"] = False
        define_opt["rijk"] = True

        assert run_itest(
            "ridft",
            define_opt,
            structure,
            "ridft_{}_rik".format(structure),
            ScfOutput,
            test_data=test_data,
        )

    @pytest.mark.parametrize("structure", structures)
    def test_run_ridft_marij(self, structure, test_data):
        define_opt = get_define_template("ridft")
        define_opt["ri"] = True
        define_opt["marij"] = True

        assert run_itest(
            "ridft",
            define_opt,
            structure,
            "ridft_{}_marij".format(structure),
            ScfOutput,
            test_data=test_data,
        )

    @pytest.mark.parametrize("structure", structures)
    def test_run_ridft_fermi(self, structure, test_data):
        define_opt = get_define_template("ridft")
        define_opt["ri"] = False
        define_opt["rijk"] = True
        datagroups_options = {"fermi": "tmstrt=500 tmend=50 tmfac=0.9 hlcrt=0.2"}

        assert run_itest(
            "ridft",
            define_opt,
            structure,
            "ridft_{}_fermi".format(structure),
            ScfOutput,
            datagroups_options=datagroups_options,
            test_data=test_data,
        )
