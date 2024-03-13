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
from turbomoleio.output.files import (
    EgradOutput,
    GradOutput,
    RelaxOutput,
    ScfOutput,
    StatptOutput,
)
from turbomoleio.testing import run_itest

structures = ["h2o", "nh3"]


@pytest.mark.integration
class TestRelax:
    @pytest.mark.parametrize("structure", structures)
    def test_run_dscf_grad_statpt(self, structure):
        dp = get_define_template("dscf")
        dp["desy"] = True
        dp["ired"] = True

        assert run_itest(
            ["dscf", "grad", "statpt"],
            dp,
            structure,
            "dscf_grad_statpt_{}_sym".format(structure),
            [ScfOutput, GradOutput, StatptOutput],
        )

    @pytest.mark.parametrize("structure", structures)
    def test_run_ridft_rdgrad_relax(self, structure):
        dp = get_define_template("ridft")
        dp["desy"] = True
        dp["ired"] = True

        assert run_itest(
            ["ridft", "rdgrad", "relax"],
            dp,
            structure,
            "ridft_rdgrad_relax_{}_sym".format(structure),
            [ScfOutput, GradOutput, RelaxOutput],
        )

    @pytest.mark.parametrize("structure", structures)
    def test_run_dscf_egrad_relax(self, structure):
        dp = get_define_template("dscf_escf")
        dp["desy"] = False
        dp["ired"] = False

        assert run_itest(
            ["dscf", "egrad", "relax"],
            dp,
            structure,
            "dscf_escf_relax_{}_nosym".format(structure),
            [ScfOutput, EgradOutput, RelaxOutput],
        )

    @pytest.mark.parametrize("structure", structures)
    def test_run_ridft_egrad_statpt_ex_irrep(self, structure):
        dp = get_define_template("ridft_escf")
        dp["desy"] = True
        dp["ired"] = True
        dp["ex_all_states"] = None
        dp["ex_irrep_states"] = {"a1": 1}

        assert run_itest(
            ["ridft", "egrad", "statpt"],
            dp,
            structure,
            "ridft_egrad_statpt_{}_ex_irrep".format(structure),
            [ScfOutput, EgradOutput, StatptOutput],
        )
