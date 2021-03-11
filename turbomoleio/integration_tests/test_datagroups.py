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
Tests for objects related to DataGroups.
Complementary to what is present in the other integration tests, addresses mainly the
writing methods that are not tested elsewhere.
"""
import pytest
import shutil
import subprocess
import os

import numpy as np
from turbomoleio.testfiles.utils import temp_dir, ItestConfig
from turbomoleio.core.control import Control, Shells, mdgo, sdg, cdg, sdgo, cpc
from turbomoleio.input.define import DefineRunner
from turbomoleio.output.parser import Parser
from turbomoleio.input.utils import get_define_template


def run_define_runner(define_opt, mol_path):
    if mol_path:
        shutil.copyfile(mol_path, 'coord')
    dr = DefineRunner(define_opt, timeout=ItestConfig.define_timeout)
    assert dr.run_full()


def run_tm(executable):
    output_file = "{}.log".format(executable)
    with open(output_file, 'w') as f_out:
        process = subprocess.Popen([executable], stdin=subprocess.PIPE, stdout=f_out,
                               stderr=subprocess.PIPE, encoding='utf-8')
        program_std_out, program_std_err = process.communicate()

        ret_code = process.wait()

        if ret_code or "ended normally" not in program_std_err:
            raise RuntimeError("Executable {} has failed with return code {}".format(executable, ret_code))


@pytest.mark.integration
class TestControl:

    @pytest.mark.parametrize('molecule_filename', ["nh3"])
    def test_mod_control(self, molecule_filepath, delete_tmp_dir):
        with temp_dir(delete_tmp_dir) as tmp_dir:
            dp = get_define_template("ridft")
            dp["desy"] = True
            dp["functional"] = "b-p"
            run_define_runner(define_opt=dp, mol_path=molecule_filepath)

            assert sdgo("dft", "functional").strip() == "b-p"
            assert sdg("thime").strip() == "5"
            cdg("thime", "6")
            mdgo("dft", {"functional": "functional pbe"})

            run_tm("ridft")

            p = Parser.from_file("ridft.log")
            assert p.integral["thime"] == 6
            assert p.density_functional_data["functional_name"] == "pbe"

    @pytest.mark.parametrize('molecule_filename', ["nh3"])
    def test_cpc(self, molecule_filepath, delete_tmp_dir):
        with temp_dir(delete_tmp_dir) as tmp_dir:
            dp = get_define_template("ridft")
            run_define_runner(define_opt=dp, mol_path=molecule_filepath)
            run_tm("ridft")

            process = subprocess.Popen(["cpc", "cpc_orig"], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE, encoding='utf-8')
            process.communicate()

            ret_code = process.wait()
            if ret_code:
                raise RuntimeError("original cpc failed with return code {}".format(ret_code))

            cpc("cpc_tmio")
            s_orig = set(os.listdir("cpc_orig"))
            s_orig.remove("wherefrom")
            s_tmio = set(os.listdir("cpc_tmio"))
            assert s_orig == s_tmio

    @pytest.mark.parametrize('molecule_filename', ["nh3"])
    def test_remove_last_energy(self, molecule_filepath, delete_tmp_dir):
        """
        Tests the remove_last_energy method of Control.
        First runs twice an energy-grad-relax loop, then removes the last
        energy and runs escf. The output of escf contains the gs energy and
        it reads it from the energy file. The test checks that the value from
        the escf.log output corresponds to the first value of the energy and not
        to the one that should have been removed.
        """
        with temp_dir(delete_tmp_dir) as tmp_dir:
            dp = get_define_template("ridft_escf")
            run_define_runner(define_opt=dp, mol_path=molecule_filepath)
            for i in range(2):
                run_tm("ridft")
                run_tm("rdgrad")
                run_tm("statpt")

            c = Control.from_file()
            en = c.energy.total
            assert len(c.energy.total) == 2
            c.remove_last_energy()
            run_tm("escf")
            p = Parser.from_file("escf.log")
            gs_en_from_escf = p.escf_gs_total_en
            assert gs_en_from_escf == pytest.approx(en[0])
            assert not gs_en_from_escf == pytest.approx(en[1])

    @pytest.mark.parametrize('molecule_filename', ["nh3"])
    def test_remove_last_gradient(self, molecule_filepath, delete_tmp_dir):
        """
        Tests the remove_last_gradient method of Control.
        First runs twice an energy-grad-relax loop, then removes the last
        energy and runs escf. The output of escf contains the gs energy and
        it reads it from the energy file. The test checks that the value from
        the escf.log output corresponds to the first value of the energy and not
        to the one that should have been removed.
        """
        with temp_dir(delete_tmp_dir) as tmp_dir:
            dp = get_define_template("ridft")
            run_define_runner(define_opt=dp, mol_path=molecule_filepath)

            run_tm("ridft")
            run_tm("rdgrad")
            run_tm("statpt")
            run_tm("ridft")
            run_tm("rdgrad")

            c = Control.from_file()
            g = c.gradient
            assert len(g.gradients) == 2
            c.remove_last_gradient()
            run_tm("statpt")
            p = Parser.from_file("statpt.log")
            g_from_statp = p.relax_gradient_values
            assert np.allclose(g.norms[0], g_from_statp["norm_cartesian"])
            assert not np.allclose(g.norms[1], g_from_statp["norm_cartesian"])


@pytest.mark.integration
class TestShells:

    @pytest.mark.parametrize('molecule_filename', ["nh3"])
    def test_mod_shells(self, molecule_filepath, delete_tmp_dir):
        with temp_dir(delete_tmp_dir) as tmp_dir:
            dp = get_define_template("ridft")
            dp["desy"] = False
            dp["functional"] = "b-p"
            run_define_runner(define_opt=dp, mol_path=molecule_filepath)

            s = Shells.from_file("closed")
            s.states.pop(-1)
            s.occupations.pop(-1)
            cdg("closed shells", s.to_datagroup())

            run_tm("ridft")

            p = Parser.from_file("ridft.log")

            assert p.electrostatic_moments["charge"] == pytest.approx(2)
