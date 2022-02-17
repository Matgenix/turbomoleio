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
Integration tests for the DefineRunner.

These tests will only require a small interaction with the turbomole executable and
so running through the usual run_itest method is not required or even not possible.
"""

import os
import shutil

import pytest

from turbomoleio.core.datagroups import DataGroups
from turbomoleio.input.define import (
    DefineExpectError,
    DefineIredError,
    DefineParameterError,
    DefineRunner,
)
from turbomoleio.input.utils import get_define_template
from turbomoleio.output.states import States
from turbomoleio.testfiles.utils import ItestConfig, temp_dir


@pytest.mark.integration
class TestDefineRunner:
    def run_define_runner(self, define_opt, mol_path, run_type="full"):
        if mol_path:
            shutil.copyfile(mol_path, "coord")
        dr = DefineRunner(define_opt, timeout=ItestConfig.define_timeout)
        # run should return True
        if run_type == "full":
            assert dr.run_full()
        elif run_type == "internal":
            assert dr.run_update_internal_coords()
        elif run_type == "mo":
            assert dr.run_generate_mo_files()
        else:
            raise ValueError("wrong run_type {}".format(run_type))

    @pytest.mark.parametrize("molecule_filename", ["nh3"])
    def test_wrong_method(self, molecule_filepath, delete_tmp_dir):
        define_opt = get_define_template("dscf")
        define_opt["method"] = "wrong_method"

        with temp_dir(delete_tmp_dir):
            with pytest.raises(DefineParameterError):
                self.run_define_runner(define_opt, molecule_filepath)

    @pytest.mark.parametrize("molecule_filename", ["nh3"])
    def test_wrong_sym(self, molecule_filepath, delete_tmp_dir):

        define_opt = get_define_template("dscf")
        define_opt["sym"] = "wrong_sym"

        with temp_dir(delete_tmp_dir):
            with pytest.raises(DefineParameterError):
                self.run_define_runner(define_opt, molecule_filepath)

    @pytest.mark.parametrize("molecule_filename", ["nh3"])
    def test_wrong_basis(self, molecule_filepath, delete_tmp_dir):

        define_opt = get_define_template("dscf")
        define_opt["basis"] = "wrong_basis"

        with temp_dir(delete_tmp_dir):
            with pytest.raises(DefineParameterError):
                self.run_define_runner(define_opt, molecule_filepath)

    @pytest.mark.parametrize("molecule_filename", ["nh3"])
    def test_wrong_xc_func(self, molecule_filepath, delete_tmp_dir):
        define_opt = get_define_template("dscf")
        define_opt["functional"] = "pbexxxx"

        with temp_dir(delete_tmp_dir):
            with pytest.raises(DefineParameterError):
                self.run_define_runner(define_opt, molecule_filepath)

    @pytest.mark.parametrize("molecule_filename", ["ch4_fixed_ired_not_conv"])
    def test_not_converged_ired(self, molecule_filepath, delete_tmp_dir):
        define_opt = get_define_template("dscf")
        define_opt["ired"] = True

        with temp_dir(delete_tmp_dir):
            with pytest.raises(DefineIredError):
                self.run_define_runner(define_opt, molecule_filepath)

    @pytest.mark.parametrize("molecule_filename", ["nh3"])
    def test_eof(self, molecule_filepath, delete_tmp_dir):
        """
        Tests an example where expect fails and an error is raised due to EOF.
        """

        define_opt = get_define_template("ridft_rimp2")
        define_opt["method"] = "ccsdt"
        define_opt["mp2energy"] = False

        with temp_dir(delete_tmp_dir):
            with pytest.raises(DefineExpectError):
                self.run_define_runner(define_opt, molecule_filepath)

    @pytest.mark.parametrize("molecule_filename", ["nh3"])
    def test_timeout(self, molecule_filepath, delete_tmp_dir):
        """
        Tests an example where expect fails and an error is raised due to TIMEOUT.
        """

        define_opt = get_define_template("dscf")
        define_opt["scfconv"] = 1e-5
        with temp_dir(delete_tmp_dir):
            with pytest.raises(DefineExpectError):
                self.run_define_runner(define_opt, molecule_filepath)

    @pytest.mark.parametrize("molecule_filename", ["nh3"])
    def test_generate_mo(self, molecule_filepath, delete_tmp_dir):
        define_opt = get_define_template("dscf")
        with temp_dir(delete_tmp_dir):
            # run the full to generate all the files
            self.run_define_runner(define_opt, molecule_filepath)
            os.remove("mos")
            self.run_define_runner({}, None, "mo")
            assert os.path.isfile("mos")

    @pytest.mark.parametrize("molecule_filename", ["nh3"])
    def test_generate_mo_symmetry(self, molecule_filepath, delete_tmp_dir):
        define_opt = get_define_template("dscf")
        define_opt["desy"] = True
        with temp_dir(delete_tmp_dir):
            # run the full to generate all the files
            self.run_define_runner(define_opt, molecule_filepath)
            s = States.from_file("control")
            # check that there are other irreps coming from the mos
            assert len(s.irreps) != 1 or s.irreps[0] != "a"
            # do not remove mos to trigger another path in DefineRunner
            self.run_define_runner({"sym": "c1"}, None, "mo")
            assert os.path.isfile("mos")
            s_new = States.from_file("control")
            assert s_new.irreps[0] == "a"

    @pytest.mark.parametrize("molecule_filename", ["nh3"])
    def testrun_update_internal_coords(self, molecule_filepath, delete_tmp_dir):
        define_opt = get_define_template("dscf")
        define_opt["ired"] = True
        with temp_dir(delete_tmp_dir):
            # run the full to generate all the files
            self.run_define_runner(define_opt, molecule_filepath)
            States.from_file("control")
            dg = DataGroups.from_file("coord")
            dg.kdg("redundant")
            dg.to_file("coord")
            # do not remove mos to trigger another path in DefineRunner
            self.run_define_runner({"ired": True}, None, "internal")
            dg_new = DataGroups.from_file("coord")
            assert dg_new.sdg("redundant")
