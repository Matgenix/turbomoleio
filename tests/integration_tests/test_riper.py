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
from turbomoleio.output.files import ScfOutput
from turbomoleio.testing import run_itest

structures = ["graphene"]


@pytest.mark.integration
class TestRiper:
    @pytest.mark.parametrize("structure_filename", structures)
    def test_run_riper(self, structure_filepath):
        structure_filename = os.path.basename(structure_filepath)
        # Turbomole does not seem to recognize periodic and cell within the coord file
        # even when they are reference in the control file with "$periodic file=coord"
        # and "$cell file=coord"
        periodic = sdg("periodic", structure_filepath)
        cell = sdg("cell", structure_filepath)
        assert run_itest(
            ["riper"],
            get_define_template("ridft"),
            structure_filename,
            "ridft_riper_{}_std".format(structure_filename),
            [ScfOutput],
            datagroups_options={"periodic": periodic, "cell": cell},
        )
