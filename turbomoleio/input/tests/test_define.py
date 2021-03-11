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
import pexpect
import signal
import subprocess
import pytest
from pexpect.popen_spawn import PopenSpawn
from unittest import mock
from monty.os import makedirs_p
from monty.collections import dict2namedtuple

from turbomoleio.testfiles.utils import temp_dir, touch_file
from turbomoleio.input.define import DefineError, DefineParameterError, DefineExpectError, UnsupportedDefineStateError
from turbomoleio.input.define import DefineIredError
from turbomoleio.input.define import DefineRunner
from turbomoleio.core.control import Control


class TestExceptions(object):
    """
    Tests the Define specific exceptions
    """

    def test_create(self):
        de = DefineError("test error")
        dpe = DefineParameterError("test error")
        die = DefineIredError("test error")
        dee = DefineExpectError("test error", pattern=["pattern1", "pattern2"])
        udse = UnsupportedDefineStateError("test error")


@pytest.fixture(scope="function", params=[True, False], ids=["use_popen", "use_pty"])
def dr_data(request):
    use_popen = request.param

    dr = DefineRunner({}, use_popen=use_popen)

    if use_popen:
        spawn_class = pexpect.popen_spawn.PopenSpawn
        dr.define = mock.Mock(pexpect.popen_spawn.PopenSpawn)
        dr.define.proc = mock.Mock(subprocess.Popen)
    else:
        spawn_class = pexpect.spawn
        dr.define = mock.Mock(pexpect.spawn)

    expect_mock = dr.define.expect
    sendline_mock = dr.define.sendline

    data = dict2namedtuple(dict(use_popen=use_popen, dr=dr, expect_mock=expect_mock,
                                sendline_mock=sendline_mock, spawn_class=spawn_class))

    return data


# unrealistic example of output from the state selection menu.
# designed to trigger the different ifs and exceptions in the code
str_state_selection_menu = b"""STATE SELECTION MENU

 IRREP |  #STATES  | #SELECTED |
 ----------------------------------------------------------------
  a1   |       24  |        0  | z xx+yy+zz
  a2   |        9  |        0  | Lz
  b1   |       19  |        0  | x Ly xz
  b2   |       13  |        0  | y Lx yz
  a    |    22923  |        0  | x y z Lx Ly Lz 1.85xx+1.85yy+zz
       |           |           | xz yz

 SELECT IRREP AND NUMBER OF STATES
 ENTER ? FOR HELP, * OR Q TO QUIT, & TO GO BACK
"""

match_ex_state_selection = re.search(b"----.*SELECT IRREP AND NUMBER OF STATES", str_state_selection_menu, re.DOTALL)

str_response_menu_u_off = b""" MAIN MENU FOR RESPONSE CALCULATIONS

 OPTION | STATUS | DESCRIPTION
 -------------------------------------------------------------------
 urpa   | off    | SPIN-UNRESTRICTED RPA EXCITATIONS (TDHF OR TDDFT)
 ucis   | off    | SPIN-UNRESTRICTED TDA EXCITATIONS (CI SINGLES)
 polly  | off    | STATIC POLARIZABILITY
 dynpol | off    | DYNAMIC POLARIZABILITY

 ENTER <OPTION> TO SWITCH ON/OFF OPTION, * OR q TO QUIT
"""

match_ex_response_u_off = re.search(b"----.*ENTER <OPTION> TO SWITCH ON/OFF OPTION", str_response_menu_u_off, re.DOTALL)

str_response_menu_urpa_on = bytes(str_response_menu_u_off).replace(b"urpa   | off", b"urpa   | on ")

match_ex_response_urpa_on = re.search(b"----.*ENTER <OPTION> TO SWITCH ON/OFF OPTION", str_response_menu_urpa_on, re.DOTALL)

str_response_menu_ucis_on = bytes(str_response_menu_u_off).replace(b"ucis   | off", b"ucis   | on ")

match_ex_response_ucis_on = re.search(b"----.*ENTER <OPTION> TO SWITCH ON/OFF OPTION", str_response_menu_ucis_on, re.DOTALL)

str_response_menu_st_off = b""" MAIN MENU FOR RESPONSE CALCULATIONS

 OPTION | STATUS | DESCRIPTION
 -------------------------------------------------------------------
 rpas   | off    | RPA SINGLET EXCITATIONS (TDHF OR TDDFT)
 ciss   | off    | TDA SINGLET EXCITATIONS (CI SINGLES)
 rpat   | off    | RPA TRIPLET EXCITATIONS (TDHF OR TDDFT)
 cist   | off    | TDA TRIPLET EXCITATIONS (CI SINGLES)
 polly  | off    | STATIC POLARIZABILITY
 dynpol | off    | DYNAMIC POLARIZABILITY
 single | off    | SINGLET STABILITY ANALYSIS
 triple | off    | TRIPLET STABILITY ANALYSIS
 nonrel | off    | NON-REAL STABILITY ANALYSIS

 ENTER <OPTION> TO SWITCH ON/OFF OPTION, * OR q TO QUIT
"""

match_ex_response_st_off = re.search(b"----.*ENTER <OPTION> TO SWITCH ON/OFF OPTION", str_response_menu_st_off, re.DOTALL)

str_response_menu_rpas_on = bytes(str_response_menu_st_off).replace(b"rpas   | off", b"rpas   | on ")

match_ex_response_rpas_on = re.search(b"----.*ENTER <OPTION> TO SWITCH ON/OFF OPTION", str_response_menu_rpas_on, re.DOTALL)

str_response_menu_polly_on = bytes(str_response_menu_st_off).replace(b"polly  | off", b"polly  | on ")

match_ex_response_polly_on = re.search(b"----.*ENTER <OPTION> TO SWITCH ON/OFF OPTION", str_response_menu_polly_on, re.DOTALL)

str_response_menu_wrong = b""" MAIN MENU FOR RESPONSE CALCULATIONS

 OPTION | STATUS | DESCRIPTION
 -------------------------------------------------------------------
 polly  | off    | STATIC POLARIZABILITY
 dynpol | off    | DYNAMIC POLARIZABILITY

 ENTER <OPTION> TO SWITCH ON/OFF OPTION, * OR q TO QUIT
"""

match_ex_response_wrong = re.search(b"----.*ENTER <OPTION> TO SWITCH ON/OFF OPTION", str_response_menu_wrong, re.DOTALL)


class TestDefineRunner(object):
    """
    Unit tests for DefineRunner. Mainly based on mocking.
    All methods are tested separately mocking the define and expect command
    """

    @staticmethod
    def assert_sendline_calls(sendline_mock, calls):
        """
        Checks if the calls to sendline were performed and no else.
        All the calls and the order should match exactly.

        Args:
            sendline_mock: the Mock object for sendline
            calls: a list of strings with the expected calls to sendline

        Returns:
            None
        """
        assert sendline_mock.mock_calls == [mock.call(c) for c in calls]

    @staticmethod
    def assert_sendline_has_calls(sendline_mock, calls):
        """
        Checks if at least these calls to sendline were performed. Other calls
        may have been done, but this will not be checked.
        Based on mock method 'assert_has_calls'

        Args:
            sendline_mock: the Mock object for sendline
            calls: a list of strings with the expected calls to sendline

        Returns:
            None
        """
        sendline_mock.assert_has_calls([mock.call(c) for c in calls], any_order=False)

    def test_init(self, dr_data):
        dr = DefineRunner({}, use_popen=dr_data.use_popen)
        assert dr._spawn == dr_data.spawn_class

    def test_is_alive(self, dr_data):
        dr = dr_data.dr
        dr_data.dr._is_alive

        # the called function depends on use_popen. Need to distinguish
        if dr_data.use_popen:
            dr.define.proc.poll.assert_called_once_with()
        else:
            dr.define.isalive.assert_called_once_with()

    def test_close(self, dr_data):
        dr = dr_data.dr
        dr.define = mock.Mock(dr_data.spawn_class)
        with mock.patch("turbomoleio.input.define.DefineRunner._is_alive", new_callable=mock.PropertyMock) as is_alive:
            is_alive.return_value = True
            dr._close()

            # the called function depends on use_popen. Need to distinguish
            if dr_data.use_popen:
                dr.define.kill.assert_called_once_with(sig=signal.SIGKILL)
            else:
                dr_data.dr.define.close.assert_called_once_with(force=True)

    def test_expect(self, dr_data):
        dr_data.expect_mock.return_value = 0
        case = dr_data.dr._expect(["pattern"], timeout=2, action="test action")
        assert case == 0

        assert len(dr_data.dr.history) == 1
        message = dr_data.dr.history[0]
        assert "Expect" in message
        assert "Action" in message
        assert "pattern" in message
        assert "test action" in message
        dr_data.expect_mock.assert_called_once_with(["pattern"], timeout=2)

    def test_expect_raise_timeout(self, dr_data):
        dr_data.expect_mock.side_effect = pexpect.exceptions.TIMEOUT("timeout")
        with mock.patch("turbomoleio.input.define.DefineRunner._is_alive", new_callable=mock.PropertyMock) as is_alive:
            is_alive.return_value = True
            with pytest.raises(DefineExpectError):
                dr_data.dr._expect(["pattern"], timeout=2, action="test action")

    def test_expect_raise_eof(self, dr_data):
        dr_data.expect_mock.side_effect = pexpect.exceptions.EOF("eof")
        with mock.patch("turbomoleio.input.define.DefineRunner._is_alive", new_callable=mock.PropertyMock) as is_alive:
            is_alive.return_value = False
            with pytest.raises(DefineExpectError):
                dr_data.dr._expect(["pattern"], timeout=2, action="test action")

    def test_sendline(self, dr_data):
        dr_data.sendline_mock.return_value = None
        dr_data.dr._sendline("command", action="test action")

        assert len(dr_data.dr.history) == 1
        message = dr_data.dr.history[0]
        assert "Sendline" in message
        assert "Action" in message
        assert "command" in message
        assert "test action" in message
        dr_data.sendline_mock.assert_called_once_with("command")

    @pytest.mark.parametrize('dir_name', ['mos_nh3'])
    def test_run_full_1(self, dr_data, mo_dirpath, delete_tmp_dir):
        with temp_dir(delete_tmp_dir) as tmp_dir:
            dr_data.dr.workdir = tmp_dir
            other_control = os.path.join(tmp_dir, "control2")
            touch_file(other_control)
            dr_data.dr.parameters = {"method": "dft", "ex_method": "rpa", "usemo": other_control,
                                     "copymo": mo_dirpath, "use_cosmo": True}
            dr_data.expect_mock.side_effect = iter(int, 1)
            # mock these functions to not having failures due to case == 0
            dr_data.dr._copy_mo_files = mock.MagicMock()
            dr_data.dr._define_symmetry = mock.MagicMock()
            dr_data.dr._add_cosmo = mock.MagicMock()
            dr_data.dr._spawn = mock.MagicMock()
            dr_data.dr._spawn.return_value = dr_data.dr.define
            se = [match_ex_response_st_off, match_ex_response_rpas_on, match_ex_state_selection]
            type(dr_data.dr.define).match = mock.PropertyMock(side_effect=se)

            #create an empty control file for the postprocess
            Control.empty().to_file()

            assert dr_data.dr.run_full()
            self.assert_sendline_calls(dr_data.sendline_mock, dr_data.dr.sent_commands)

    def test_run_full_2(self, dr_data, delete_tmp_dir):
        with temp_dir(delete_tmp_dir) as tmp_dir:
            dr_data.dr.workdir = tmp_dir
            dr_data.dr.parameters = {"method": "dft", "ex_method": "rpa"}
            dr_data.expect_mock.side_effect = iter(int, 1)
            dr_data.dr._copy_mo_files = mock.MagicMock()
            dr_data.dr._spawn = mock.MagicMock()
            dr_data.dr._spawn.return_value = dr_data.dr.define
            se = [match_ex_response_st_off, match_ex_response_rpas_on, match_ex_state_selection]
            type(dr_data.dr.define).match = mock.PropertyMock(side_effect=se)
            if dr_data.use_popen:
                # set poll to None so that kill is called and cover all cases
                dr_data.dr.define.proc.poll.return_value = None

            # create an empty control file for the postprocess
            Control.empty().to_file()

            assert dr_data.dr.run_full()
            if dr_data.use_popen:
                dr_data.dr.define.kill.assert_called_once_with(sig=signal.SIGKILL)
            else:
                dr_data.dr.define.close.assert_called_once_with(force=True)
            self.assert_sendline_calls(dr_data.sendline_mock, dr_data.dr.sent_commands)

    @pytest.mark.parametrize('dir_name', ['mos_nh3'])
    def test_run_update_internal_coords(self, dr_data, mo_dirpath, delete_tmp_dir):
        with temp_dir(delete_tmp_dir) as tmp_dir:
            dr_data.dr.workdir = tmp_dir
            dr_data.dr.parameters = {"copymo": mo_dirpath}
            dr_data.expect_mock.side_effect = [0, 0, 1, 0, 1, 0, 0, 0]
            dr_data.dr._spawn = mock.MagicMock()
            dr_data.dr._spawn.return_value = dr_data.dr.define
            if dr_data.use_popen:
                # set poll to None so that kill is called and cover all cases
                dr_data.dr.define.proc.poll.return_value = None
            assert dr_data.dr.run_update_internal_coords()
            self.assert_sendline_has_calls(dr_data.sendline_mock, ["sy c3v ", "*", "qq"])
            if dr_data.use_popen:
                dr_data.dr.define.kill.assert_called_once_with(sig=signal.SIGKILL)
            else:
                dr_data.dr.define.close.assert_called_once_with(force=True)
            self.assert_sendline_calls(dr_data.sendline_mock, dr_data.dr.sent_commands)

    @pytest.mark.parametrize('dir_name', ['mos_nh3'])
    def test_run_update_internal_coords_error(self, dr_data, mo_dirpath, delete_tmp_dir):
        with temp_dir(delete_tmp_dir) as tmp_dir:
            dr_data.dr.workdir = tmp_dir
            dr_data.dr.parameters = {"copymo": mo_dirpath}
            dr_data.expect_mock.side_effect = [0, 0, 1, 0, 1, 1, 0, 0]
            dr_data.dr._spawn = mock.MagicMock()
            dr_data.dr._spawn.return_value = dr_data.dr.define
            if dr_data.use_popen:
                # set poll to None so that kill is called and cover all cases
                dr_data.dr.define.proc.poll.return_value = None

            with pytest.raises(DefineIredError):
                dr_data.dr.run_update_internal_coords()

    def test_run_generate_mo_files(self, dr_data, delete_tmp_dir):
        with temp_dir(delete_tmp_dir) as tmp_dir:
            dr_data.dr.workdir = tmp_dir
            dr_data.dr.parameters = {}
            dr_data.expect_mock.side_effect = [0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0]
            dr_data.dr._spawn = mock.MagicMock()
            dr_data.dr._spawn.return_value = dr_data.dr.define
            if dr_data.use_popen:
                # set poll to None so that kill is called and cover all cases
                dr_data.dr.define.proc.poll.return_value = None
            dr_data.dr.run_generate_mo_files()
            self.assert_sendline_calls(dr_data.sendline_mock, ["", "n", "n", "y", "eht", "y", "0", "y", "qq"])
            if dr_data.use_popen:
                dr_data.dr.define.kill.assert_called_once_with(sig=signal.SIGKILL)
            else:
                dr_data.dr.define.close.assert_called_once_with(force=True)
            self.assert_sendline_calls(dr_data.sendline_mock, dr_data.dr.sent_commands)

    def test_run_generate_mo_files_sym(self, dr_data, delete_tmp_dir):
        with temp_dir(delete_tmp_dir) as tmp_dir:
            dr_data.dr.workdir = tmp_dir
            # use control here just to let the test pass
            dr_data.dr.parameters = {"sym": "c1", "usemo": "control"}
            touch_file("control")
            dr_data.expect_mock.side_effect = [0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0]
            dr_data.dr._spawn = mock.MagicMock()
            dr_data.dr._spawn.return_value = dr_data.dr.define
            if dr_data.use_popen:
                # set poll to None so that kill is called and cover all cases
                dr_data.dr.define.proc.poll.return_value = None
            dr_data.dr.run_generate_mo_files()
            self.assert_sendline_calls(dr_data.sendline_mock, ["", "y", "sy c1 ", "*", "n", "y", "use control", "", "qq"])
            if dr_data.use_popen:
                dr_data.dr.define.kill.assert_called_once_with(sig=signal.SIGKILL)
            else:
                dr_data.dr.define.close.assert_called_once_with(force=True)
            self.assert_sendline_calls(dr_data.sendline_mock, dr_data.dr.sent_commands)

    def test_get_bin_path(self, dr_data):
        assert dr_data.dr._get_bin_path() == "define"

    def test_set_metric(self, dr_data, delete_tmp_dir):
        with temp_dir(delete_tmp_dir) as tmp_dir:
            # first create the control with the specified metric
            dr_data.dr.parameters = {"metric": 2}
            dr_data.dr._set_metric()
            with open("control") as f:
                s = f.read()
                assert "metric 2" in s

            # now modify the existing one
            dr_data.dr.parameters = {"metric": 1}
            dr_data.dr._set_metric()
            with open("control") as f:
                s = f.read()
                assert "metric 1" in s
                assert "2" not in s

    def test_initialize_control_incompatibility(self, dr_data):
        dr_data.expect_mock.side_effect = [1]
        with pytest.raises(DefineError):
            dr_data.dr._initialize_control(update=True)

    def test_initialize_control_1(self, dr_data):
        dr_data.expect_mock.side_effect = [1, 0]
        dr_data.dr._initialize_control(update=False)

    def test_initialize_control_2(self, dr_data):
        dr_data.expect_mock.side_effect = [1, 0]
        dr_data.dr.parameters = {"title": "test title"}
        dr_data.dr._initialize_control(update=False)
        self.assert_sendline_calls(dr_data.sendline_mock, ["", "test title"])

    def test_geometry_menu_new_coords_1(self, dr_data):
        dr_data.dr.parameters = {"desy": True, "ired": True}
        dr_data.expect_mock.side_effect = [0, 0, 0, 0, 0]
        dr_data.dr._geometry_menu(new_coords=True)
        self.assert_sendline_calls(dr_data.sendline_mock, ["a coord", "c", "desy", "ired"])

    def test_geometry_menu_new_coords_2(self, dr_data):
        dr_data.dr.parameters = {"desy": True, "ired": True, "desy_eps": 0.0001}
        dr_data.expect_mock.side_effect = [0, 1, 0, 0]
        dr_data.dr._geometry_menu(new_coords=True)
        self.assert_sendline_calls(dr_data.sendline_mock, ["a coord", "desy 0.0001", "ired"])

    def test_geometry_menu_new_coords_3(self, dr_data):
        # here desy will be ignored
        dr_data.dr.parameters = {"desy": True, "ired": True, "sym": "c2v", "sym_eps": 0.01}
        dr_data.expect_mock.side_effect = [0, 1, 1, 0]
        dr_data.dr._geometry_menu(new_coords=True)
        self.assert_sendline_calls(dr_data.sendline_mock, ["a coord", "sy c2v 0.01", "ired"])

    def test_geometry_menu_no_new_coords(self, dr_data):
        dr_data.dr.parameters = {"desy": False, "ired": False}
        dr_data.expect_mock.side_effect = [1, 1, 0]
        dr_data.dr._geometry_menu(new_coords=False, update=True)
        self.assert_sendline_calls(dr_data.sendline_mock, ["y"])

    @pytest.mark.parametrize('dir_name', ['mos_nh3'])
    def test_geometry_menu_copymo(self, dr_data, mo_dirpath):
        dr_data.dr.parameters = {"desy": True, "ired": True, "sy": None,
                              "copymo": mo_dirpath}
        dr_data.expect_mock.side_effect = [0, 0, 0, 1, 0]
        dr_data.dr._geometry_menu(new_coords=True)
        self.assert_sendline_calls(dr_data.sendline_mock, ["a coord", "c", "sy c3v ", "ired"])

    def test_geometry_menu_incompatibility_1(self, dr_data):
        dr_data.expect_mock.side_effect = [1]
        with pytest.raises(DefineError):
            dr_data.dr._geometry_menu(new_coords=True)

    def test_geometry_menu_incompatibility_2(self, dr_data):
        dr_data.expect_mock.side_effect = [0]
        with pytest.raises(DefineError):
            dr_data.dr._geometry_menu(new_coords=False)

    def test_general_menu_wrong_method(self, dr_data):
        dr_data.dr.parameters = {"method": "wrong method"}
        with pytest.raises(DefineParameterError):
            dr_data.dr._general_menu()

    def test_general_menu_dft(self, dr_data):
        dr_data.dr.parameters = {"method": "DFt"}
        dr_data.expect_mock.side_effect = [0, 0, 0, 0, 0, 0]
        dr_data.dr._general_menu()
        self.assert_sendline_has_calls(dr_data.sendline_mock, ["dft", "on"])

    def test_general_menu_hf(self, dr_data):
        dr_data.dr.parameters = {"method": "hf"}
        dr_data.expect_mock.side_effect = [0, 0, 0, 0, 0]
        dr_data.dr._general_menu()
        self.assert_sendline_has_calls(dr_data.sendline_mock, ["dft", "off"])

    def test_general_menu_ccsdt(self, dr_data):
        dr_data.dr.parameters = {"method": "ccsdt"}
        dr_data.expect_mock.side_effect = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        dr_data.dr._general_menu()
        self.assert_sendline_has_calls(dr_data.sendline_mock, ["cc", "freeze"])

    def test_add_atomic_coordinates_from_file_cartesian(self, dr_data):
        dr_data.dr.parameters = {"coord_file": "another_coord"}
        dr_data.expect_mock.side_effect = [0, 0]
        dr_data.dr._add_atomic_coordinates_from_file()
        self.assert_sendline_calls(dr_data.sendline_mock, ["a another_coord", "c"])

    def test_add_atomic_coordinates_from_file(self, dr_data):
        dr_data.dr.parameters = {}
        dr_data.expect_mock.side_effect = [1]
        dr_data.dr._add_atomic_coordinates_from_file()
        self.assert_sendline_calls(dr_data.sendline_mock, ["a coord"])

    def test_existing_coordinates_error(self, dr_data):
        with pytest.raises(DefineError):
            dr_data.dr._existing_coordinates(case=0, update=True)

    def test_existing_coordinates(self, dr_data):
        dr_data.expect_mock.side_effect = [1, 0]
        dr_data.dr._existing_coordinates(case=1, update=True)
        self.assert_sendline_calls(dr_data.sendline_mock, ["y"])

    def test_existing_coordinates_no_update(self, dr_data):
        dr_data.expect_mock.side_effect = [1]
        dr_data.dr._existing_coordinates(case=1, update=False)
        self.assert_sendline_calls(dr_data.sendline_mock, ["n"])

    def test_existing_coordinates_cartesian(self, dr_data):
        dr_data.expect_mock.side_effect = [0, 0]
        dr_data.dr._existing_coordinates(case=2, update=True)
        self.assert_sendline_calls(dr_data.sendline_mock, ["c", "y"])

    def test_existing_coordinates_unsupported(self, dr_data):
        with pytest.raises(UnsupportedDefineStateError):
            dr_data.dr._existing_coordinates(case=3, update=True)

    def test_define_symmetry(self, dr_data):
        dr_data.dr.parameters = {"sym": "c2v", "sym_eps": 0.01}
        dr_data.expect_mock.side_effect = [1]
        dr_data.dr._define_symmetry()
        self.assert_sendline_calls(dr_data.sendline_mock, ["sy c2v 0.01"])

    def test_define_symmetry_error(self, dr_data):
        dr_data.dr.parameters = {"sym": "wrong_symbol"}
        dr_data.expect_mock.side_effect = [0]
        with pytest.raises(DefineParameterError):
            dr_data.dr._define_symmetry()

    def test_determine_symmetry_eps(self, dr_data):
        dr_data.dr.parameters = {"desy": True, "desy_eps": 0.01}
        dr_data.expect_mock.side_effect = [0]
        dr_data.dr._determine_symmetry()
        self.assert_sendline_calls(dr_data.sendline_mock, ["desy 0.01"])

    def test_determine_symmetry(self, dr_data):
        dr_data.dr.parameters = {"desy": True, "desy_eps": None}
        dr_data.expect_mock.side_effect = [0]
        dr_data.dr._determine_symmetry()
        self.assert_sendline_calls(dr_data.sendline_mock, ["desy"])

    def test_generate_internal_coordinates(self, dr_data):
        dr_data.expect_mock.side_effect = [0]
        dr_data.dr._generate_internal_coordinates()
        self.assert_sendline_calls(dr_data.sendline_mock, ["ired"])

    def test_generate_internal_coordinates_error(self, dr_data):
        dr_data.expect_mock.side_effect = [1]
        with pytest.raises(DefineIredError):
            dr_data.dr._generate_internal_coordinates()

    def test_switch_to_atomic_attribute_menu(self, dr_data):
        dr_data.expect_mock.side_effect = [1]
        dr_data.dr._switch_to_atomic_attribute_menu()
        self.assert_sendline_calls(dr_data.sendline_mock, ["*"])

    def test_switch_to_atomic_attribute_menu_internal(self, dr_data):
        dr_data.dr.parameters = {"ired": False}
        dr_data.expect_mock.side_effect = [0, 0]
        dr_data.dr._switch_to_atomic_attribute_menu()
        self.assert_sendline_calls(dr_data.sendline_mock, ["*", "no"])

    def test_switch_to_atomic_attribute_menu_internal_error(self, dr_data):
        dr_data.dr.parameters = {"ired": True}
        dr_data.expect_mock.side_effect = [0, 0]
        with pytest.raises(DefineIredError):
            dr_data.dr._switch_to_atomic_attribute_menu()

    def test_skip_atomic_attribute_menu(self, dr_data):
        dr_data.expect_mock.side_effect = [0]
        dr_data.dr._skip_atomic_attribute_menu()
        self.assert_sendline_calls(dr_data.sendline_mock, ["n"])

    def test_skip_atomic_attribute_menu_internal(self, dr_data):
        dr_data.expect_mock.side_effect = [2, 0]
        dr_data.dr._skip_atomic_attribute_menu()
        self.assert_sendline_calls(dr_data.sendline_mock, ["no", "n"])

    def test_skip_atomic_attribute_menu_error(self, dr_data):
        dr_data.expect_mock.side_effect = [1]
        with pytest.raises(UnsupportedDefineStateError):
            dr_data.dr._skip_atomic_attribute_menu()

    def test_set_basis(self, dr_data):

        dr_data.expect_mock.side_effect = [0, 0]
        dr_data.dr._set_basis("all", "def2-SVP")
        self.assert_sendline_calls(dr_data.sendline_mock, ["b", "all def2-SVP"])

    def test_set_basis_error(self, dr_data):
        dr_data.expect_mock.side_effect = [0, 1]
        with pytest.raises(DefineParameterError):
            dr_data.dr._set_basis("all", "wrong basis")

    def test_define_basis_sets(self, dr_data):

        dr_data.dr.parameters = {"basis": "def2-SVP", "basis_atom": {1: "def2-SVP", "Ca": "def2-SVP"}}
        dr_data.expect_mock.side_effect = [0, 0, 0, 0, 0, 0]
        dr_data.dr._define_basis_sets()
        self.assert_sendline_calls(dr_data.sendline_mock, ["b", "all def2-SVP", "b", "1 def2-SVP", "b", "\"Ca\" def2-SVP"])

    def test_switch_to_molecular_orbital_definition_menu(self, dr_data):

        dr_data.expect_mock.side_effect = [0]
        dr_data.dr._switch_to_molecular_orbital_definition_menu()
        self.assert_sendline_calls(dr_data.sendline_mock, ["*"])

    def test_go_to_general_menu(self, dr_data):

        dr_data.expect_mock.side_effect = [0]
        dr_data.dr._go_to_general_menu()
        self.assert_sendline_calls(dr_data.sendline_mock, [""])

    def test_set_mo_from_control_file(self, dr_data, delete_tmp_dir):
        with temp_dir(delete_tmp_dir) as tmp_dir:
            control_path = os.path.join(tmp_dir, "control")
            touch_file(control_path)
            dr_data.dr.parameters = {"usemo": tmp_dir}
            dr_data.dr.workdir = tmp_dir
            dr_data.expect_mock.side_effect = [0]
            dr_data.dr._set_mo_from_control_file()
            self.assert_sendline_calls(dr_data.sendline_mock, ["use {}".format(control_path)])

    def test_set_mo_from_control_file_2(self, dr_data, delete_tmp_dir):
        with temp_dir(delete_tmp_dir) as tmp_dir:
            control_path = os.path.join(tmp_dir, "control")
            touch_file(control_path)
            dr_data.dr.parameters = {"usemo": tmp_dir}
            dr_data.dr.workdir = tmp_dir
            dr_data.expect_mock.side_effect = [4, 2]
            dr_data.dr._set_mo_from_control_file()
            self.assert_sendline_calls(dr_data.sendline_mock, ["use {}".format(control_path), ""])

    def test_set_mo_from_control_file_longpath(self, dr_data, delete_tmp_dir):
        with temp_dir(delete_tmp_dir) as tmp_dir:
            dir_path = os.path.join(*(["extremelytoolongfilepath"]*3))
            makedirs_p(dir_path)
            control_path = os.path.join(tmp_dir, dir_path, "control")
            touch_file(control_path)
            dr_data.dr.parameters = {"usemo": dir_path}
            dr_data.dr.workdir = tmp_dir
            dr_data.expect_mock.side_effect = [1, 0]
            dr_data.dr._set_mo_from_control_file()
            self.assert_sendline_calls(dr_data.sendline_mock, ["use {}".format("tmp_mo/control"), ""])

    def test_set_mo_from_control_file_error(self, dr_data):
        with pytest.raises(DefineParameterError):
            dr_data.dr.parameters = {"usemo": "/fake/path/to/nothing"}
            dr_data.dr._set_mo_from_control_file()

    def test_extended_hueckel_theory_menu_nocharge(self, dr_data):
        dr_data.dr.parameters = {"charge": None, "unpaired_electrons": None}
        dr_data.expect_mock.side_effect = [0, 0, 0, 0, 0]
        dr_data.dr._extended_hueckel_theory_menu()
        self.assert_sendline_calls(dr_data.sendline_mock, ["eht", "y", "y", "0", "y"])

    def test_extended_hueckel_theory_menu_charge(self, dr_data):
        dr_data.dr.parameters = {"charge": 1, "unpaired_electrons": None}
        dr_data.expect_mock.side_effect = [1, 0, 0, 1, 0]
        dr_data.dr._extended_hueckel_theory_menu()
        self.assert_sendline_calls(dr_data.sendline_mock, ["eht", "y", "1", "y", ""])

    def test_extended_hueckel_theory_menu_unpaired(self, dr_data):
        dr_data.dr.parameters = {"charge": 1, "unpaired_electrons": 2}
        dr_data.expect_mock.side_effect = [1, 0, 0, 0, 0, 0]
        dr_data.dr._extended_hueckel_theory_menu()
        self.assert_sendline_calls(dr_data.sendline_mock, ["eht", "y", "1", "no", "u 2", "*"])

    def test_extended_hueckel_theory_menu_error(self, dr_data):
        dr_data.dr.parameters = {"charge": 1, "unpaired_electrons": None}
        dr_data.expect_mock.side_effect = [1, 0, 1]
        with pytest.raises(DefineParameterError):
            dr_data.dr._extended_hueckel_theory_menu()

    def test_extended_hueckel_theory_menu_default(self, dr_data):
        dr_data.dr.parameters = {"charge": None, "unpaired_electrons": None}
        dr_data.expect_mock.side_effect = [2, 0, 0, 0, 0]
        dr_data.dr._extended_hueckel_theory_menu()
        self.assert_sendline_calls(dr_data.sendline_mock, ["eht", "y", "y", "0", "y"])

    def test_excited_state_menu_rpas(self, dr_data):
        se = [match_ex_response_st_off, match_ex_response_rpas_on]
        type(dr_data.dr.define).match = mock.PropertyMock(side_effect=se)
        dr_data.dr.parameters = {}
        dr_data.expect_mock.side_effect = [0, 0, 0, 2, 0]
        dr_data.dr._excited_state_menu()
        self.assert_sendline_calls(dr_data.sendline_mock, ["ex", "rpas", "*", "*", "y"])

    def test_excited_state_menu_urpa_error(self, dr_data):
        se = [match_ex_response_st_off, match_ex_response_rpas_on]
        type(dr_data.dr.define).match = mock.PropertyMock(side_effect=se)
        dr_data.dr.parameters = {"ex_method": "rpa", "ex_multi": "dublet"}
        dr_data.expect_mock.side_effect = [0, 0, 0, 2, 0]
        with pytest.raises(DefineParameterError):
            dr_data.dr._excited_state_menu()

    def test_excited_state_menu_urpa(self, dr_data):
        se = [match_ex_response_u_off, match_ex_response_urpa_on]
        type(dr_data.dr.define).match = mock.PropertyMock(side_effect=se)
        dr_data.dr.parameters = {"ex_method": "rpa", "ex_multi": "dublet"}
        dr_data.expect_mock.side_effect = [0, 0, 0, 2, 0]
        dr_data.dr._excited_state_menu()
        self.assert_sendline_calls(dr_data.sendline_mock, ["ex", "urpa", "*", "*", "y"])


    def test_excited_state_menu_polly(self, dr_data):
        se = [match_ex_response_st_off, match_ex_response_polly_on]
        type(dr_data.dr.define).match = mock.PropertyMock(side_effect=se)
        dr_data.dr.parameters = {"ex_method": "polly", "ex_multi": "singlet"}
        dr_data.expect_mock.side_effect = [0, 0, 0, 2, 0]
        dr_data.dr._excited_state_menu()
        self.assert_sendline_calls(dr_data.sendline_mock, ["ex", "polly", "*", "*", "y"])

    def test_excited_state_menu_wrong_menu(self, dr_data):
        type(dr_data.dr.define).match = mock.PropertyMock(return_value=match_ex_response_wrong)
        dr_data.dr.parameters = {"ex_method": "rpa", "ex_multi": "dublet"}
        dr_data.expect_mock.side_effect = [0, 0, 0, 2, 0]
        with pytest.raises(UnsupportedDefineStateError):
            dr_data.dr._excited_state_menu()

    def test_excited_state_menu_wrong_ex_method(self, dr_data):
        type(dr_data.dr.define).match = mock.PropertyMock(return_value=match_ex_response_st_off)
        dr_data.dr.parameters = {"ex_method": "wrong", "ex_multi": "singlet"}
        dr_data.expect_mock.side_effect = [0, 0, 0, 2, 0]
        with pytest.raises(DefineParameterError):
            dr_data.dr._excited_state_menu()

    def test_excited_state_menu_error_not_activated(self, dr_data):
        se = [match_ex_response_st_off, match_ex_response_st_off]
        type(dr_data.dr.define).match = mock.PropertyMock(side_effect=se)
        dr_data.dr.parameters = {"ex_method": "rpa", "ex_multi": "singlet"}
        dr_data.expect_mock.side_effect = [0, 0, 0, 2, 0]
        with pytest.raises(DefineError):
            dr_data.dr._excited_state_menu()

    def test_excited_state_menu_irrep(self, dr_data):
        se = [match_ex_response_st_off, match_ex_response_rpas_on, match_ex_state_selection]
        type(dr_data.dr.define).match = mock.PropertyMock(side_effect=se)
        # dr_data.dr.define.match = match_ex_state_selection
        dr_data.dr.parameters = {"ex_all_states": 10, "ex_irrep_states": {"a1": 7}}
        dr_data.expect_mock.side_effect = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        dr_data.dr._excited_state_menu()
        # use has_calls here since the iteration of the states is done over a dict and the order may not
        # be assured.
        self.assert_sendline_has_calls(dr_data.sendline_mock,
                                       ["ex", "rpas", "*", "a1 7", "a2 9", "b1 10", "b2 10", "a 10", "*", "*", "y"])

    def test_excited_state_menu_irrep_wrong_sym(self, dr_data):
        se = [match_ex_response_st_off, match_ex_response_rpas_on, match_ex_state_selection]
        type(dr_data.dr.define).match = mock.PropertyMock(side_effect=se)
        dr_data.dr.parameters = {"ex_all_states": None, "ex_irrep_states": {"z": 7}}
        dr_data.expect_mock.side_effect = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        with pytest.raises(DefineParameterError, match="No excited state could be set given the input parameters"):
            dr_data.dr._excited_state_menu()

    def test_excited_state_menu_frequency(self, dr_data):
        se = [match_ex_response_st_off, match_ex_response_rpas_on]
        type(dr_data.dr.define).match = mock.PropertyMock(side_effect=se)
        dr_data.dr.parameters = {"ex_frequency": 100}
        dr_data.expect_mock.side_effect = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
        dr_data.dr._excited_state_menu()
        self.assert_sendline_calls(dr_data.sendline_mock, ["ex", "rpas", "*", "nm", "a", "100", "*", "*", "*", "y"])

    def test_set_ri_state(self, dr_data):
        dr_data.dr.parameters = {"rijk": True}
        dr_data.expect_mock.side_effect = [0, 0, 0]
        dr_data.dr._set_ri_state()
        self.assert_sendline_calls(dr_data.sendline_mock, ["rijk", "on", ""])

    def test_set_ri_state_marij(self, dr_data):
        dr_data.dr.parameters = {"ri": True, "marij": True}
        dr_data.expect_mock.side_effect = [0, 0, 0, 0, 0]
        dr_data.dr._set_ri_state()
        self.assert_sendline_calls(dr_data.sendline_mock, ["ri", "on", "", "marij", ""])

    def test_set_dft_options_use_dft(self, dr_data):
        dr_data.dr.parameters = {"gridsize": "m3", "functional": "pbe"}
        dr_data.expect_mock.side_effect = [0, 0, 0, 0, 0]
        dr_data.dr._set_dft_options(use_dft=True)
        self.assert_sendline_calls(dr_data.sendline_mock, ["dft", "on", "func pbe", "grid m3", ""])

    def test_set_dft_options_wrong_xc(self, dr_data):
        dr_data.dr.parameters = {"gridsize": "m3", "functional": "pbexxxx"}
        dr_data.expect_mock.side_effect = [0, 0, 1]
        with pytest.raises(DefineParameterError):
            dr_data.dr._set_dft_options(use_dft=True)

    def test_set_dft_options_no_dft(self, dr_data):
        dr_data.dr.parameters = {}
        dr_data.expect_mock.side_effect = [0, 0, 0, 0, 0]
        dr_data.dr._set_dft_options(use_dft=False)
        self.assert_sendline_calls(dr_data.sendline_mock, ["dft", "off", ""])

    def test_set_mp2_options(self, dr_data):

        dr_data.dr.parameters = {"maxcor": 300, "ex_mp2": {}, "use_f12": False, "use_f12*": False, "maxiter": 200}
        dr_data.expect_mock.side_effect = [0, 0, 0, 0, 0, 0, 0, 0]
        dr_data.dr._set_mp2_options(energy_only=True, method="mp2")
        self.assert_sendline_calls(dr_data.sendline_mock,
                                   ["cc", "freeze", "*", "cbas", "*", "memory", "300", "ricc2", "mp2",
                                    "maxiter 200", "*", "*"])

    def test_set_mp2_options_geoopt(self, dr_data):
        dr_data.dr.parameters = {"maxcor": 300, "ex_mp2": {}, "use_f12": False, "use_f12*": False, "maxiter": 200}
        dr_data.expect_mock.side_effect = [0, 0, 0, 0, 0, 0, 0, 0]
        dr_data.dr._set_mp2_options(energy_only=False, method="mp2")
        self.assert_sendline_calls(dr_data.sendline_mock,
                                   ["cc", "freeze", "*", "cbas", "*", "memory", "300", "ricc2", "geoopt mp2",
                                    "maxiter 200", "*", "*"])

    def test_set_mp2_options_mp2(self, dr_data):
        dr_data.dr.parameters = {"maxcor": None, "ex_mp2": {"a1": [10, 1]},
                                 "use_f12": False, "use_f12*": False, "maxiter": None}
        dr_data.expect_mock.side_effect = [0, 0, 0, 0, 0, 3, 0, 0, 0]
        dr_data.dr._set_mp2_options(energy_only=True, method="adc(2)")
        self.assert_sendline_calls(dr_data.sendline_mock,
                                   ["cc", "freeze", "*", "cbas", "*", "exci", "irrep=a1 multiplicity=10 nexc=1",
                                    "spectrum  states=all", "exprop  states=all unrelaxed", "*", "ricc2", "adc(2)",
                                    "*", "*"])

    def test_set_mp2_options_mp2_error(self, dr_data):
        # fail due at ex_mp2
        dr_data.dr.parameters = {"maxcor": None, "ex_mp2": {"a1": [10, 1]},
                                 "use_f12": False, "use_f12*": False, "maxiter": None}
        dr_data.expect_mock.side_effect = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        with pytest.raises(DefineParameterError):
            dr_data.dr._set_mp2_options(energy_only=True, method="adc(2)")

    def test_set_mp2_options_f12(self, dr_data, delete_tmp_dir):
        with temp_dir(delete_tmp_dir) as tmp_dir:
            c = Control.empty()
            c.to_file("control")
            dr_data.dr.parameters = {"maxcor": None, "ex_mp2": {}, "use_f12": True, "use_f12*": True, "maxiter": None,
                                     "method": "ccsdt"}
            dr_data.expect_mock.side_effect = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            dr_data.dr._set_mp2_options(energy_only=True, method="ccsd(t)")
            self.assert_sendline_calls(dr_data.sendline_mock, ["cc", "freeze", "*", "cbas", "*", "f12",
                                                               "*", "cabs", "*", "ricc2", "ccsd(t)","*", "*"])
            c = Control.from_file("control")
            assert c.sdgo("rir12", "ccsdapprox").strip() == "ccsd(f12*)"

    def test_switch_to_scf_menu(self, dr_data):

        dr_data.expect_mock.side_effect = [0]
        dr_data.dr._switch_to_scf_menu()
        self.assert_sendline_calls(dr_data.sendline_mock, ["scf"])

    def test_set_scf_options(self, dr_data):

        dr_data.dr.parameters = {"scfconv": 0.0001, "scfiterlimit": 100}
        dr_data.expect_mock.side_effect = [0, 0, 0, 0, 0, 0]
        dr_data.dr._set_scf_options()
        self.assert_sendline_calls(dr_data.sendline_mock, ["scf", "conv", "0.0001", "iter", "100", ""])

    def test_quit_general_menu(self, dr_data):

        dr_data.dr._quit_general_menu()
        self.assert_sendline_calls(dr_data.sendline_mock, ["*"])

    @pytest.mark.parametrize('dir_name', ['mos_nh3'])
    def test_copy_mo_files_mo(self, dr_data, mo_dirpath, delete_tmp_dir):
        with temp_dir(delete_tmp_dir) as tmp_dir:
            dr_data.dr.parameters = {"copymo": mo_dirpath}
            dr_data.dr.workdir = tmp_dir
            touch_file("mos")
            dr_data.dr._copy_mo_files()
            with open("mos") as f:
                s = f.read()
                assert len(s) > 10
            assert not os.path.isfile("alpha")

    @pytest.mark.parametrize('dir_name', ['alpha_beta_nh3'])
    def test_copy_mo_files_alpha_beta(self, dr_data, mo_dirpath, delete_tmp_dir):
        with temp_dir(delete_tmp_dir) as tmp_dir:
            dr_data.dr.parameters = {"copymo": mo_dirpath}
            dr_data.dr.workdir = tmp_dir
            touch_file("alpha")
            touch_file("beta")
            dr_data.dr._copy_mo_files()
            with open("alpha") as f:
                s = f.read()
                assert len(s) > 10
            assert not os.path.isfile("mos")

    @pytest.mark.parametrize('dir_name', ['mos_nh3'])
    def test_copy_mo_files_ab_mos(self, dr_data, mo_dirpath, delete_tmp_dir):
        with temp_dir(delete_tmp_dir) as tmp_dir:
            dr_data.dr.parameters = {"copymo": mo_dirpath}
            dr_data.dr.workdir = tmp_dir
            touch_file("alpha")
            touch_file("beta")
            dr_data.dr._copy_mo_files()
            with open("alpha") as f:
                s = f.read()
                assert len(s) > 10
            assert not os.path.isfile("mos")

    def test_copy_mo_files_error(self, dr_data, delete_tmp_dir):
        with temp_dir(delete_tmp_dir) as tmp_dir:
            dr_data.dr.parameters = {"copymo": "fake_path"}
            dr_data.dr.workdir = tmp_dir
            with pytest.raises(DefineParameterError):
                dr_data.dr._copy_mo_files()

    def test_add_cosmo(self, dr_data, delete_tmp_dir):
        with temp_dir(delete_tmp_dir) as tmp_dir:
            with open("control", "w") as f:
                f.write("$end")

            dr_data.dr.parameters = {"epsilon": 0.1, "nppa": 1, "nspa": 1, "disex": 0.1, "rsolv": 0.1,
                                     "routf": 0.1, "cavity": "closed"}
            dr_data.dr._add_cosmo()

            with open("control") as f:
                s = f.read()

            assert "epsilon" in s
            assert "nppa" in s
            assert "disex" in s
            assert "rsolv" in s
            assert "routf" in s
            assert "cavity" in s

    def test_dump_command_file(self, delete_tmp_dir):
        dr = DefineRunner({})
        dr.sent_commands = ['', 'dscf', 'a coord', 'desy', 'ired', '*', 'b', 'all def-SV(P)', '*',
                            'eht', 'y', '0', 'y', 'dft', 'on', 'func b-p', '', 'scf', 'iter', '200', '', '*']

        with temp_dir(delete_tmp_dir) as tmp_dir:
            dr.dump_command_file(filepath="test_file")
            with open("test_file") as f:
                lines = f.readlines()
            commands = [l.strip() for l in lines]

            assert dr.sent_commands == commands

    def test_postprocess(self, dr_data, delete_tmp_dir):
        with temp_dir(delete_tmp_dir) as tmp_dir:
            dr_data.dr.parameters = {"disp": "DFT-D3", "use_cosmo": True}

            with open("control", "w") as f:
                f.write("$end")

            dr_data.dr._postprocess()

            with open("control") as f:
                s = f.read()

            assert "disp3" in s
            assert "cosmo" in s
