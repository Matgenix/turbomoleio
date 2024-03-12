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

"""
The parsing of all the Data objects is done in the test_files.py module and
is not repeated here. In this module only specific functions of the Data
objects are tested.
"""

import pytest

from turbomoleio.output.data import (
    AoforceRotationalData,
    AoforceVibrationalData,
    EscfData,
    EscfIterationData,
    GeometryData,
    MP2Data,
    MP2Results,
    RelaxConvergenceData,
    RelaxData,
    RunData,
    ScfEnergiesData,
    ScfIterationData,
    StatptData,
    TurbomoleData,
)
from turbomoleio.output.parser import Parser
from turbomoleio.testing import TM_VERSIONS, has_matplotlib


@pytest.mark.parametrize("tm_version", TM_VERSIONS)
def test_str(test_data, tm_version):

    path = test_data / "outputs" / tm_version / "dscf" / "h2o_std" / "dscf.log"
    with open(path) as f:
        log = f.read()

    td = TurbomoleData.from_string(log)
    s = str(td)
    assert "Content of the TurbomoleData object" in s
    assert "dscf" in s


@pytest.mark.parametrize("tm_version", TM_VERSIONS)
def test_ScfIterationData(test_data, tm_version):
    path = test_data / "outputs" / tm_version / "dscf" / "h2o_std" / "dscf.log"
    sid = ScfIterationData.from_file(path)

    if has_matplotlib():
        assert sid.plot_energies(show=False)


def test_AoforceVibrationalData(test_data):
    path = test_data / "outputs" / "TM_v7.3" / "aoforce" / "aceton_full" / "aoforce.log"
    avd = AoforceVibrationalData.from_file(path)

    assert avd.n_negative_freqs(tol=0.1) == 1
    assert avd.n_zero_freqs(tol=0.1) == 6
    assert avd.n_positive_freqs(tol=0.1) == 23

    df = avd.get_freqs_df()
    assert df["frequency"][0] == pytest.approx(-113.79)


def test_TurbomoleData():
    p = Parser("noheader")
    td = TurbomoleData.from_parser(p)
    assert td is None


def test_RunData():
    p = Parser("noheadernotimings")
    rd = RunData.from_parser(p)
    assert rd is None


@pytest.mark.parametrize("tm_version", TM_VERSIONS)
def test_RunData_only_header(test_data, tm_version):
    path = test_data / "outputs" / tm_version / "dscf" / "h2o_std" / "dscf.log"
    with open(path) as f:
        mystring = f.read(1000)
    p = Parser(mystring)
    rd = RunData.from_parser(p)

    assert rd is not None
    assert rd.host
    assert rd.start_time
    assert not rd.end_time
    assert not rd.cpu_time
    assert not rd.wall_time


@pytest.mark.parametrize("tm_version", TM_VERSIONS)
def test_RunData_only_timings(test_data, tm_version):
    path = test_data / "outputs" / tm_version / "dscf" / "h2o_std" / "dscf.log"
    with open(path) as f:
        lines = f.readlines()
    mystring = "".join(lines[-10:])
    p = Parser(mystring)
    rd = RunData.from_parser(p)

    assert rd is not None
    assert not rd.host
    assert not rd.start_time
    assert rd.end_time
    assert rd.cpu_time
    assert rd.wall_time


def test_ScfEnergiesData_empty():
    p = Parser("nothing here")
    sed = ScfEnergiesData.from_parser(p)
    assert sed is None


@pytest.mark.parametrize("tm_version", TM_VERSIONS)
def test_ScfEnergiesData(test_data, tm_version):
    # Skip tm version 7.3 as riper references were not there yet.
    if tm_version in ["TM_v7.3"]:
        return
    dscf_path = test_data / "outputs" / tm_version / "dscf" / "h2o_std" / "dscf.log"
    riper_path = test_data / "outputs" / tm_version / "riper" / "bulk_LiH" / "riper.log"
    with open(dscf_path, "r") as f:
        dscf_string = f.read()
    with open(riper_path, "r") as f:
        riper_string = f.read()
    twologs_string = dscf_string + riper_string
    with pytest.raises(
        RuntimeError,
        match=r"Found scf energy data from dscf/ridft as well as from riper.",
    ):
        p = Parser(twologs_string)
        ScfEnergiesData.from_parser(p)


def test_GeometryData_empty():
    p = Parser("nothing here")
    gd = GeometryData.from_parser(p)
    assert gd is None


def test_EscfIterationData_empty():
    p = Parser("nothing here")
    eid = EscfIterationData.from_parser(p)
    assert eid is None


def test_EscfData_empty():
    p = Parser("nothing here")
    ed = EscfData.from_parser(p)
    assert ed is None


def test_StatptData_empty():
    p = Parser("nothing here")
    sd = StatptData.from_parser(p)
    assert sd is None


def test_RelaxData_empty():
    p = Parser("nothing here")
    rd = RelaxData.from_parser(p)
    assert rd is None


def test_RelaxConvergenceData_empty():
    p = Parser("nothing here")
    rcd = RelaxConvergenceData.from_parser(p)
    assert rcd is None


def test_AoforceRotationalData_empty():
    p = Parser("nothing here")
    ard = AoforceRotationalData.from_parser(p)
    assert ard is None


def test_AoforceVibrationalData_empty():
    p = Parser("nothing here")
    avd = AoforceVibrationalData.from_parser(p)
    assert avd is None


def test_MP2Results_empty():
    p = Parser("nothing here")
    mp2r = MP2Results.from_parser(p)
    assert mp2r is None


def test_MP2Data_init():
    mp2d = MP2Data(energy_only=True)
    assert mp2d.energy_only is True
