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

import os

import pytest

from turbomoleio.core.control import sdg
from turbomoleio.input.utils import get_define_template
from turbomoleio.output.files import JobexOutput
from turbomoleio.testing import run_itest

structures = ["h2o", "nh3"]


@pytest.mark.integration
class TestJobex:
    @pytest.mark.parametrize("structure", structures)
    def test_run_jobex_dscf(self, structure, test_data):
        dp = get_define_template("dscf")
        dp["desy"] = True
        dp["ired"] = True

        assert run_itest(
            "jobex",
            dp,
            structure,
            "jobex_dscf_{}_sym".format(structure),
            JobexOutput,
            arguments="-c 2",
            test_data=test_data,
        )

    @pytest.mark.parametrize("structure_filename", ["graphene"])
    def test_run_jobex_riper(self, structure_filepath, test_data):
        structure_filename = os.path.basename(structure_filepath)
        dp = get_define_template("ridft")
        periodic = sdg("periodic", structure_filepath)
        cell = sdg("cell", structure_filepath)
        assert run_itest(
            "jobex",
            dp,
            structure_filename,
            "jobex_riper_{}".format(structure_filename),
            JobexOutput,
            arguments="-c 2",
            datagroups_options={"periodic": periodic, "cell": cell, "optcell": ""},
            test_data=test_data,
        )
