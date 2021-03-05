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

import os
import re
from unittest import mock
from fractions import Fraction
import shutil


import pytest

from turbomoleio.output.states import State, States, EigerOutput, EigerRunner, get_mos_energies
from turbomoleio.testfiles.utils import temp_dir, assert_MSONable


def test_state():
    s = State(0.05, "a1", 1, Fraction(1), "a")
    s1 = State(0.05, "a1", 1, Fraction(1), "a")
    s2 = State(0.05, "t1", 1, Fraction(2))
    assert s != 1
    assert s == s1
    assert s.max_occupation == 1
    assert re.match(r"Eigenvalue: [.\d]+, irrep: a1, index: 1, spin: a, occupation: 1", str(s))
    assert s2.has_fractional_occ
    assert s2.max_occupation == 6
    assert_MSONable(s)
    assert_MSONable(s2)


def test_parsing_failures():
    """Tests only the failing cases, the working are tested in the other functions"""
    with pytest.raises(RuntimeError):
        get_mos_energies("wrong string")

    with pytest.raises(RuntimeError):
        EigerOutput.from_string("Nr. \n wrong format")


@pytest.mark.parametrize('dir_name', ['alpha_beta_nh3'])
def test_uhf(mo_dirpath, delete_tmp_dir):
    with temp_dir(delete_tmp_dir):
        for f in ["alpha", "beta", "control"]:
            shutil.copy2(os.path.join(mo_dirpath, f), f)

        s = States.from_file()

        assert len(s) == 28
        assert s.n_states == 28
        assert s[0].eigenvalue == pytest.approx(-14.457427310212)

        assert_MSONable(s)

        eiger_out = EigerOutput.from_file(os.path.join(mo_dirpath, "eiger"))
        eiger_wrong_val = EigerOutput.from_file(os.path.join(mo_dirpath, "eiger_wrong_val"))
        eiger_wrong_length = EigerOutput.from_file(os.path.join(mo_dirpath, "eiger_wrong_length"))
        eiger_wrong_gap = EigerOutput.from_file(os.path.join(mo_dirpath, "eiger_wrong_gap"))
        eiger_wrong_nel = EigerOutput.from_file(os.path.join(mo_dirpath, "eiger_wrong_nelec"))
        eiger_wrong_occ = EigerOutput.from_file(os.path.join(mo_dirpath, "eiger_wrong_occ"))
        eiger_wrong_match = EigerOutput.from_file(os.path.join(mo_dirpath, "eiger_wrong_match"))

        assert eiger_out.compare_states(s) is None
        assert "gap" in eiger_wrong_gap.compare_states(s)
        assert "states" in eiger_wrong_length.compare_states(s)
        assert "eigenvalue" in eiger_wrong_val.compare_states(s)
        assert "electrons" in eiger_wrong_nel.compare_states(s)
        assert "occupation" in eiger_wrong_occ.compare_states(s)
        assert "match" in eiger_wrong_match.compare_states(s)

        assert s.is_uhf
        assert not s.has_hole

        assert s.total_electrons == 9
        assert s.irreps[0] == "a1"
        assert s.occupations[0] == 1
        assert s.spins[0] == "a"
        assert s.eigenvalues[-1] == pytest.approx(2.1425788573507)
        assert s.irrep_indices[0] == 1

        assert str(s)

        s.pop(0)
        assert len(s) == 27

        s[0] = s[1]
        assert s[0] == s[1]
        s.insert(0, s[-1])
        assert s[0] == s[-1]

        # test filter
        s = States.from_file()
        assert len(s.filter_states(irrep="z")) == 0
        filtered = s.filter_states(spin="b", irrep="a1")
        assert filtered[0].spin == "b"
        assert filtered[0].irrep == "a1"

        # test get_shells
        s = States.from_file()
        with pytest.raises(ValueError):
            s.get_shells(spin=None)

        shell = s.get_shells(spin="a")
        assert shell.occupations[0] == 1
        assert shell.states[0] == ("a1", 1)


@pytest.mark.parametrize('dir_name', ['mos_nh3'])
def test_closed(mo_dirpath, delete_tmp_dir):
    with temp_dir(delete_tmp_dir):
        for f in ["mos", "control"]:
            shutil.copy2(os.path.join(mo_dirpath, f), f)

        s = States.from_file()

        assert len(s) == 14
        assert s[0].eigenvalue == pytest.approx(-13.968369531218)

        assert_MSONable(s)

        eiger_out = EigerOutput.from_file(os.path.join(mo_dirpath, "eiger"))

        assert eiger_out.compare_states(s) is None

        assert not s.is_uhf
        assert not s.has_hole

        assert s.total_electrons == 10
        assert s.irreps[0] == "a1"
        assert s.occupations[0] == 2
        assert s.spins[0] is None
        assert s.eigenvalues[-1] == pytest.approx(2.55164050421)
        assert s.irrep_indices[0] == 1

        assert str(s)

        # test filter
        s = States.from_file()
        with pytest.raises(ValueError):
            s.filter_states(spin="a", irrep="a1")

        assert len(s.filter_states(irrep="z")) == 0
        filtered = s.filter_states(irrep="a1")
        assert filtered[0].irrep == "a1"

        # test get_shells
        s = States.from_file()
        with pytest.raises(ValueError):
            s.get_shells(spin="a")

        shell = s.get_shells()
        assert shell.occupations[0] == 2
        assert shell.states[0] == ("a1", 1)

@pytest.mark.parametrize('dir_name', ['mos_nh3'])
def test_hole(mo_dirpath, delete_tmp_dir):
    with temp_dir(delete_tmp_dir):
        for f in ["mos", "control"]:
            shutil.copy2(os.path.join(mo_dirpath, f), f)

        s = States.from_file()

        # create a fake hole
        high_state = s.pop(-1)
        high_state.eigenvalue = -0.5
        s.insert(2, high_state)
        assert s.has_hole

        with pytest.raises(RuntimeError):
            s.generate_lowest_filled_states()

        lowest = s.generate_lowest_filled_states(allow_fractional=True)
        assert not lowest.has_hole
        assert lowest[2].irrep_index == 6

        lowest = s.generate_lowest_filled_states(only_occupied=True, allow_fractional=True)
        assert len(lowest) == 4

        lowest = s.generate_lowest_filled_states(allow_fractional=True, reorder_irrep_index=True)
        assert lowest[2].irrep_index == 1


def test_error_homo():
    s = States([])
    with pytest.raises(RuntimeError):
        s.homo_index


def test_no_empty_states():
    s = States([])

    assert s.lumo_index is None
    assert s.lumo is None
    assert s.gap is None


@pytest.mark.parametrize('dir_name', ['mos_nh3'])
def test_eiger_runner(mo_dirpath, delete_tmp_dir):
    with mock.patch("subprocess.run") as mock_run:
        er = EigerRunner()

        with pytest.raises(ValueError):
            er.to_file("test")
        with pytest.raises(ValueError):
            er.get_eiger_output()

        with open(os.path.join(mo_dirpath, "eiger")) as f:
            out = f.read()
        m = mock.MagicMock()
        m.stdout = out.encode("utf-8")
        mock_run.return_value = m
        er.run()
        assert er.get_eiger_output()

        with temp_dir(delete_tmp_dir):
            er.to_file("out_test")
            assert os.path.isfile("out_test")

