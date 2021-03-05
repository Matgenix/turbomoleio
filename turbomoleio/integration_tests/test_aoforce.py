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

import pytest

from turbomoleio.testfiles.utils import run_itest
from turbomoleio.input.utils import get_define_template
from turbomoleio.output.files import ScfOutput, AoforceOutput

structures = ['h2o', 'nh3']


@pytest.mark.integration
class TestAoforce:

    @pytest.mark.parametrize("structure", structures)
    def test_run_ridft_aoforce_nosym(self, structure):
        dp = get_define_template("ridft")
        dp["desy"] = False

        assert run_itest(["ridft", "aoforce"], dp,
                         structure, "ridft_aoforce_{}_nosym".format(structure), [ScfOutput, AoforceOutput])

    @pytest.mark.parametrize("structure", structures)
    def test_run_dscf_aoforce_sym(self, structure):
        dp = get_define_template("dscf")
        dp["desy"] = True

        assert run_itest(["dscf", "aoforce"], dp,
                         structure, "dscf_aoforce_{}_sym".format(structure), [ScfOutput, AoforceOutput])
