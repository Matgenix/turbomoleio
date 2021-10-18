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
The parsing of all the Data objects is done in the test_files.py module and
is not repeated here. In this module only specific functions of the Data
objects are tested.
"""
import os
import pytest

from turbomoleio.output.data import TurbomoleData, ScfIterationData, AoforceVibrationalData
from turbomoleio.testfiles.utils import has_matplotlib
from turbomoleio.testfiles.utils import TM_VERSIONS


@pytest.mark.parametrize("tm_version", TM_VERSIONS)
def test_str(testdir, tm_version):

    path = os.path.join(testdir, "outputs", tm_version, "dscf", "h2o_std", "dscf.log")
    with open(path) as f:
        log = f.read()

    td = TurbomoleData.from_string(log)
    s = str(td)
    assert "Content of the TurbomoleData object" in s
    assert "dscf" in s


@pytest.mark.parametrize("tm_version", TM_VERSIONS)
def test_ScfIterationData(testdir, tm_version):
    path = os.path.join(testdir, "outputs", tm_version, "dscf", "h2o_std", "dscf.log")
    sid = ScfIterationData.from_file(path)

    if has_matplotlib():
        assert sid.plot_energies(show=False)


def test_AoforceVibrationalData(testdir):
    path = os.path.join(testdir, "outputs", "TM_v7.3", "aoforce", "aceton_full", "aoforce.log")
    avd = AoforceVibrationalData.from_file(path)

    assert avd.n_negative_freqs(tol=0.1) == 1
    assert avd.n_zero_freqs(tol=0.1) == 6
    assert avd.n_positive_freqs(tol=0.1) == 23

    df = avd.get_freqs_df()
    assert df["frequency"][0] == pytest.approx(-113.79)
