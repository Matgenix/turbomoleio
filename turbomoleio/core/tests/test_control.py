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
import shutil
import tempfile
import numpy as np
from fractions import Fraction

import pytest

from turbomoleio.core.control import Control, sdg, kdg, adg, cdg, mdgo, sdgo, cpc
from turbomoleio.core.control import Energy, Gradient, Shells
from turbomoleio.testfiles.utils import temp_dir, assert_MSONable, has_matplotlib


def check_equivalent_dg(dg1, dg2):
    """
    Compares two different strings that should represent the same variables
    It will compare the lines as a set and ignore repetitions of white spaces.
    The order of the lines is not important
    Args:
        dg1 (str): the first data group
        dg2 (str): the second data group

    Returns:
        True if the data group are equivalent
    """

    dg1 = dg1.strip()
    dg2 = dg2.strip()
    s1 = set()
    for l in dg1.splitlines():
        l = l.strip()
        if not l:
            continue
        l = re.sub(r"\s+", " ", l)
        m = re.search(r"^(.+?)\((.+?)\)", l)
        f = m.group(2).replace(" ", "")
        s1.add((m.group(1), Fraction(f)))

    s2 = set()
    for l in dg2.splitlines():
        l = l.strip()
        if not l:
            continue
        l = re.sub(r"\s+", " ", l)
        m = re.search(r"^(.+?)\((.+?)\)", l)
        f = m.group(2).replace(" ", "")
        s2.add((m.group(1), Fraction(f)))

    return s1 == s2


class TestEnergy:

    def test_energy_scf(self, testdir, delete_tmp_dir):
        with temp_dir(delete_tmp_dir):
            shutil.copy2(os.path.join(testdir, 'control', 'energy_test-energy'), 'energy')
            energy = Energy.from_file()
            assert len(energy.scf) == 4
            assert energy.n_steps == 4
            assert energy.scf[0] == pytest.approx(-10.48350884610)
            assert energy.scfkin[2] == pytest.approx(37.86647200156)
            assert energy.scfpot[1] == pytest.approx(-59.35016151450)
            assert energy.delta_e == pytest.approx(9.99999999577)
            assert energy.mp2 is None
            assert np.all(energy.scf == energy.total)

            assert_MSONable(energy)

            if has_matplotlib():
                assert energy.plot(show=False)

            energy_string = '      SCF               SCFKIN            SCFPOT\n' \
                            '1   -10.48350884610    39.86550149370   -50.34901033980\n' \
                            '2   -20.48350872544    38.86665278906   -59.35016151450\n' \
                            '\n' \
                            '3   -30.48350875500    37.86647200156   -68.34998075656\n' \
                            '4   -40.48350875077    36.86649728669   -77.35000603747\n'
            energy = Energy.from_string(energy_string)
            assert len(energy.scf) == 4
            assert energy.n_steps == 4
            assert energy.scf[0] == pytest.approx(-10.48350884610)
            assert energy.scfkin[2] == pytest.approx(37.86647200156)
            assert energy.scfpot[1] == pytest.approx(-59.35016151450)
            assert energy.delta_e == pytest.approx(9.99999999577)
            assert energy.mp2 is None
            assert np.all(energy.scf == energy.total)

            assert_MSONable(energy)

            if has_matplotlib():
                assert energy.plot(show=False)

    def test_energy_mp2(self, testdir, delete_tmp_dir):
        with temp_dir(delete_tmp_dir):
            shutil.copy2(os.path.join(testdir, 'control', 'energy_test-mp2'), 'energy')
            energy = Energy.from_file()
            assert len(energy.scf) == 1
            assert energy.scf[0] == pytest.approx(-75.93827460490)
            assert energy.scfkin[0] == pytest.approx(75.9346807)
            assert energy.scfpot[0] == pytest.approx(-151.872955)
            assert energy.mp2[0] == pytest.approx(-0.1869672633349)
            assert energy.mp2[0] + energy.scf[0] == energy.total[0]
            assert energy.delta_e is None

            assert_MSONable(energy)

            if has_matplotlib():
                assert energy.plot(show=False)


class TestGradient:

    def test_gradient(self, testdir, delete_tmp_dir):
        with temp_dir(delete_tmp_dir):
            shutil.copy2(os.path.join(testdir, 'control', 'gradient_test-gradient'), 'gradient')
            g = Gradient.from_file()
            assert len(g.molecules) == 6
            assert g.gradients[0,1,0] == pytest.approx(-0.02705094211381)
            assert g.energies[-1] == pytest.approx(-76.3449460148)
            assert g.norms[1] == pytest.approx(0.013946816251742454)
            assert g.max_gradients[-1] == pytest.approx(3.8028925027729e-05)
            assert g.n_steps == 6
            assert g.last_grad_norm == pytest.approx(5.145114095400708e-05)
            assert g.last_grad_max == pytest.approx(3.8028925027729e-05)
            assert g.periodicity == 0
            assert g.lattice_vectors is None
            assert g.lattice_gradients is None

            assert_MSONable(g)

            if has_matplotlib():
                assert g.plot(show=False)

    def test_empty_gradient(self, testdir):
        g = Gradient(gradients=[], energies=[], molecules=[])
        assert g.last_grad_max is None
        assert g.last_grad_norm is None

    def test_fail_parsing(self, testdir):
        with open(os.path.join(testdir, 'control', 'gradient_test-gradient')) as f:
            grad = f.read()

        lines = grad.splitlines()
        lines[6] = "0.00000000000000D+00 "
        with pytest.raises(RuntimeError):
            Gradient.from_string("\n".join(lines))

    def test_periodic(self, testdir, delete_tmp_dir):
        # 1 dimensional
        with temp_dir(delete_tmp_dir):
            shutil.copy2(os.path.join(testdir, 'control', 'periodic_1D_C_chain.control'), "control")
            shutil.copy2(os.path.join(testdir, 'control', 'periodic_1D_C_chain.energy'), "energy")
            shutil.copy2(os.path.join(testdir, 'control', 'periodic_1D_C_chain.gradient'), "gradient")
            g = Gradient.from_file(filename="control")
            assert g.periodicity == 1
            assert g.lattice_vectors is not None
            assert g.lattice_gradients is not None
            assert len(g.lattice_vectors) == 12
            assert len(g.lattice_gradients) == 12
            assert len(g.lattice_vectors[0]) == 1
            assert len(g.lattice_vectors[0][0]) == 1
            assert g.lattice_vectors[0][0][0] == pytest.approx(3.7794522492515)
            assert len(g.lattice_vectors[11]) == 1
            assert len(g.lattice_vectors[11][0]) == 1
            assert g.lattice_vectors[11][0][0] == pytest.approx(4.8021148447)
            assert g.lattice_vectors[5][0][0] == pytest.approx(4.7512307475)
            assert len(g.lattice_gradients[0]) == 1
            assert len(g.lattice_gradients[0][0]) == 1
            assert g.lattice_gradients[0][0][0] == pytest.approx(-0.92130064839626)
            assert len(g.lattice_gradients[11]) == 1
            assert len(g.lattice_gradients[11][0]) == 1
            assert g.lattice_gradients[11][0][0] == pytest.approx(6.3810585401318e-05)
            assert g.lattice_gradients[4][0][0] == pytest.approx(-0.080897501039401)

        # 2 dimensional
        with temp_dir(delete_tmp_dir):
            shutil.copy2(os.path.join(testdir, 'control', 'periodic_2D_graphene.control'), "control")
            shutil.copy2(os.path.join(testdir, 'control', 'periodic_2D_graphene.energy'), "energy")
            shutil.copy2(os.path.join(testdir, 'control', 'periodic_2D_graphene.gradient'), "gradient")
            g = Gradient.from_file(filename="control")
            assert g.periodicity == 2
            assert g.lattice_vectors is not None
            assert g.lattice_gradients is not None
            assert len(g.lattice_vectors) == 10
            assert len(g.lattice_gradients) == 10
            assert len(g.lattice_vectors[0]) == 2
            assert len(g.lattice_vectors[0][0]) == 2
            assert g.lattice_vectors[0][0] == pytest.approx([4.66463026, 0.0])
            assert g.lattice_vectors[0][1] == pytest.approx([-2.33231519083, 4.0396882697516])
            assert len(g.lattice_vectors[9]) == 2
            assert len(g.lattice_vectors[9][0]) == 2
            assert g.lattice_vectors[9][0] == pytest.approx([4.5809809753, 0.0])
            assert g.lattice_vectors[9][1] == pytest.approx([-2.0545069412832, 4.0953069748706])
            assert g.lattice_gradients[9][0] == pytest.approx([-1.62632065e-05, -5.49909640e-06])
            assert g.lattice_gradients[9][1] == pytest.approx([1.04586218e-05, 6.97617028e-05])

        # 3 dimensional
        with temp_dir(delete_tmp_dir):
            shutil.copy2(os.path.join(testdir, 'control', 'periodic_3D_HC.control'), "control")
            shutil.copy2(os.path.join(testdir, 'control', 'periodic_3D_HC.energy'), "energy")
            shutil.copy2(os.path.join(testdir, 'control', 'periodic_3D_HC.gradient'), "gradient")
            g = Gradient.from_file(filename="control")
            assert g.periodicity == 3
            assert g.lattice_vectors is not None
            assert g.lattice_gradients is not None
            assert len(g.lattice_vectors) == 18
            assert len(g.lattice_gradients) == 18
            assert len(g.lattice_vectors[0]) == 3
            assert len(g.lattice_vectors[0][0]) == 3
            assert g.lattice_vectors[0][0] == pytest.approx([5.09738357, 0.0, 0.0])
            assert g.lattice_vectors[0][1] == pytest.approx([-2.7319912004445, 3.9245528516757, 0.0])
            assert g.lattice_vectors[0][2] == pytest.approx([0.75633542716256, -1.6073467332536, 11.481850306654])
            assert len(g.lattice_vectors[9]) == 3
            assert len(g.lattice_vectors[9][0]) == 3
            assert g.lattice_vectors[9][0] == pytest.approx([4.7817678505, 0.0, 0.0])
            assert g.lattice_vectors[9][1] == pytest.approx([-2.3276141683287, 4.1856654539183, 0.0])
            assert g.lattice_vectors[9][2] == pytest.approx([0.89126099072162, -1.3162285716079, 11.50855367142])
            assert g.lattice_vectors[17][0] == pytest.approx([4.7537059886, 0.0, 0.0])
            assert g.lattice_vectors[17][1] == pytest.approx([-2.2801435854351, 4.1747226082797, 0.0])
            assert g.lattice_vectors[17][2] == pytest.approx([0.92551779449765, -1.3738823392983, 11.500943692659])
            assert g.lattice_gradients[17][0] == pytest.approx([-2.31848114e-05, -1.33019726e-04,  8.75872366e-05])
            assert g.lattice_gradients[17][1] == pytest.approx([1.26202299e-05, -2.91666006e-05, -1.75235859e-05])
            assert g.lattice_gradients[17][2] == pytest.approx([7.47949755e-05,  9.05389095e-05, -4.96335880e-05])


class TestShells:

    def test_closed(self, testdir, delete_tmp_dir):
        with temp_dir(delete_tmp_dir):
            shutil.copy2(os.path.join(testdir, 'control', 'control_test-subfiles'), "control")
            with pytest.raises(ValueError):
                Shells.from_file("alpha")
            s = Shells.from_file("closed")

            assert len(s.states) == 3
            assert s.states[0] == ("a1", 1)
            assert s.occupations[0] == 2
            assert "a1" in s.irreps
            assert "t2" in s.irreps
            assert len(s.irreps) == 2
            assert s.total_electrons == 10

            original = Control.from_file().sdg("closed shells")
            generated = s.to_datagroup()
            assert check_equivalent_dg(original, generated)

            assert_MSONable(s)

    def test_uhf(self, testdir):
        alpha = """
 a1      1-3                                    ( 1 )
 b1      1                                      ( 1 )
 b2      1                                      ( 1 )        
        """
        s = Shells.from_string(alpha)
        assert len(s.states) == 5
        assert s.states[0] == ("a1", 1)
        assert s.states[3] == ("b1", 1)
        assert s.occupations[0] == 1
        assert s.total_electrons == 5
        assert check_equivalent_dg(alpha, s.to_datagroup())
        assert_MSONable(s)

        beta = """
 a1      1-3                                    ( 1 )
 b1      1,3                                    ( 1 )
        """

        s = Shells.from_string(beta)
        assert len(s.states) == 5
        assert s.states[4] == ("b1", 3)
        assert s.total_electrons == 5
        assert check_equivalent_dg(beta, s.to_datagroup())

        alpha = """
a1      1,3-5,7-8                              ( 1 )
a1      9                                      ( 1 / 2 )
b2      1-2,5,7,10-12                          ( 0.5 )        
       """
        s = Shells.from_string(alpha)
        assert len(s.states) == 14
        assert s.states[0] == ("a1", 1)
        assert s.states[3] == ("a1", 5)
        assert s.states[6] == ("a1", 9)
        assert s.states[9] == ("b2", 5)
        assert s.occupations[0] == 1
        assert s.occupations[6] == Fraction("1/2")
        assert s.occupations[9] == Fraction("1/2")
        assert s.total_electrons == 10
        assert check_equivalent_dg(alpha, s.to_datagroup())

    def test_fail_parse(self):
        bad_shells="""
 a1      1-3                                    ( 1 )
 b1      1                                      ( 1 )
 a2        
        """
        with pytest.raises(RuntimeError):
            Shells.from_string(bad_shells)

    def test_fractional(self):
        alpha_frac = """
a1      1-3                                    ( 1 )
b1      1                                      ( 1 )
b2      1                                      ( 1 / 2 )        
        """
        s = Shells.from_string(alpha_frac)

        assert check_equivalent_dg(alpha_frac, s.to_datagroup())

        alpha_float = """
        a1      1-3                                    ( 1 )
        b1      1                                      ( 1 )
        b2      1                                      ( 0.5 )        
                """
        s = Shells.from_string(alpha_float)

        # here test again with alpha_frac since the fraction will be writtes as 1/2
        assert check_equivalent_dg(alpha_frac, s.to_datagroup())


class TestControl(object):
    """Testing of the Control object."""

    @pytest.mark.parametrize('control_filename', ['control_test-Control'])
    def test_from_file(self, control):
        """Testing instantiation."""

        # Check that we are actually getting an instance of Control
        assert isinstance(control, Control)

        dg_list_ref = ['$title\n',
                       '$operating system unix\n',
                       '$symmetry c1\n',
                       '$user-defined bonds    file=coord\n',
                       '$coord    file=coord\n',
                       '$optimize\n'
                       ' internal   off\n'
                       ' redundant  off\n'
                       ' cartesian  on\n'
                       ' global     off\n'
                       ' basis      off\n',
                       '$end\n']
        assert control.dg_list == dg_list_ref
        assert str(control) == ''.join(dg_list_ref)
        assert_MSONable(control)

        # Test to_file and from_file
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = os.path.join(tmpdir, 'control_test')
            control.to_file(filename=fname)
            control_from_file = Control.from_file(filename=fname)
            assert control.dg_list == control_from_file.dg_list

    @pytest.mark.parametrize('control_filename', ['control_test-Control'])
    def test_add_cosmo(self, control):
        """Testing addition of Cosmo"""
        control.add_cosmo(epsilon=0.1, nppa=1, nspa=1, disex=0.1,
                          rsolv=0.1, routf=0.1, cavity="closed")

        dg = control.show_data_group("$cosmo")
        assert "epsilon" in dg
        assert "nppa" in dg
        assert "disex" in dg
        assert "rsolv" in dg
        assert "routf" in dg
        assert "cavity" in dg

    def test_energy(self, testdir, delete_tmp_dir):
        with temp_dir(delete_tmp_dir):
            for fname in ['control', 'energy']:
                source = os.path.join(testdir, 'control', '{}_test-subfiles'.format(fname))
                shutil.copy2(source, fname)
            cc = Control.from_file()

            ee = cc.energy
            assert isinstance(ee, Energy)
            assert ee.scf[0] == pytest.approx(-40.48328750375)
            assert ee.scfkin[0] == pytest.approx(39.86406749278)
            assert ee.scfpot[0] == pytest.approx(-80.34735499653)

            os.unlink("energy")
            assert cc.energy is None

    def test_gradient(self, testdir, delete_tmp_dir):
        with temp_dir(delete_tmp_dir):
            shutil.copy2(os.path.join(testdir, 'control', 'control_test-subfiles'), "control")
            shutil.copy2(os.path.join(testdir, 'control', 'gradient_test-gradient'), 'gradient')
            c = Control.from_file()
            g = c.gradient
            assert len(g.molecules) == 6
            assert g.gradients[0,1,0] == pytest.approx(-0.02705094211381)

            os.unlink("gradient")
            assert c.gradient is None

    def test_remove_energy(self, testdir, delete_tmp_dir):
        with temp_dir(delete_tmp_dir):
            shutil.copy2(os.path.join(testdir, 'control', 'control_test-subfiles'), "control")
            shutil.copy2(os.path.join(testdir, 'control', 'energy_test-energy'), "energy")

            c = Control.from_file()
            c.remove_last_energy()
            e = Energy.from_file("energy")
            assert len(e.scf) == 3
            assert e.scf[-1] == pytest.approx(-30.48350875500)

            c.remove_last_energy(backup_suffix="_backup")
            assert os.path.isfile("energy_backup")
            e = Energy.from_file("energy")
            assert len(e.scf) == 2
            assert e.scf[-1] == pytest.approx(-20.48350872544)

            c.remove_last_energy()
            c.remove_last_energy()

            with pytest.raises(RuntimeError):
                c.remove_last_energy()

    def test_remove_energy_in_control(self, testdir, delete_tmp_dir):
        with temp_dir(delete_tmp_dir):
            shutil.copy2(os.path.join(testdir, 'control', 'control_test-energy'), "control")

            c = Control.from_file()

            with pytest.raises(ValueError):
                c.remove_last_energy(filename=None)

            c.remove_last_energy(backup_suffix="_backup")
            assert os.path.isfile("control_backup")
            e = c.energy
            assert len(e.scf) == 0

            with pytest.raises(RuntimeError):
                c.remove_last_energy()

    def test_remove_gradient(self, testdir, delete_tmp_dir):
        with temp_dir(delete_tmp_dir):
            shutil.copy2(os.path.join(testdir, 'control', 'control_test-subfiles'), "control")
            shutil.copy2(os.path.join(testdir, 'control', 'gradient_test-gradient'), "gradient")

            c = Control.from_file()
            c.remove_last_gradient()
            g = Gradient.from_file("gradient")
            assert len(g.gradients) == 5
            assert g.gradients[-1,0,2] == pytest.approx(-0.00021835002325094)

            c.remove_last_gradient(backup_suffix="_backup")
            assert os.path.isfile("gradient_backup")
            g = Gradient.from_file("gradient")
            assert len(g.gradients) == 4
            assert g.gradients[-1,0,2] == pytest.approx(0.00047140511636532)

            for i in range(4):
                c.remove_last_gradient()

            with pytest.raises(RuntimeError):
                c.remove_last_gradient()

    def test_remove_gradient_in_control(self, testdir, delete_tmp_dir):
        with temp_dir(delete_tmp_dir):
            shutil.copy2(os.path.join(testdir, 'control', 'control_test-gradient'), "control")

            c = Control.from_file()

            with pytest.raises(ValueError):
                c.remove_last_gradient(filename=None)

            c.remove_last_gradient(backup_suffix="_backup")
            assert os.path.isfile("control_backup")
            g = c.gradient
            assert len(g.gradients) == 0

            with pytest.raises(RuntimeError):
                c.remove_last_gradient()

    def test_from_metric(self):
        c = Control.from_metric(2)
        dg = c.show_data_group("$redund_inp")
        assert dg is not None
        assert "metric 2" in dg

    @pytest.mark.parametrize('control_filename', ['control_test-Control'])
    def test_disp(self, control):

        control.set_disp("DFT-D1")
        string_c = str(control)
        assert "$olddisp" in string_c
        assert "$disp" not in string_c

        control.set_disp("DFT-D2")
        string_c = str(control)
        assert "$olddisp" not in string_c
        assert "$disp" in string_c

        control.set_disp("DFT-D3")
        string_c = str(control)
        assert "$olddisp" not in string_c
        assert control.sdg("disp3").strip() == ""

        control.set_disp("DFT-D3 BJ")
        string_c = str(control)
        assert "$olddisp" not in string_c
        assert control.sdg("disp3").strip() == "bj"

        control.set_disp(None)
        string_c = str(control)
        assert "$olddisp" not in string_c
        assert "$disp" not in string_c

        with pytest.raises(ValueError):
            control.set_disp("wrong_value")

    @pytest.mark.parametrize('control_filename', ['control_test-Control'])
    def test_is_uhf(self, control):
        assert not control.is_uhf

    @pytest.mark.parametrize('control_filename', ['control_test-subfiles'])
    def test_get_shells(self, control):
        with pytest.raises(ValueError):
            control.get_shells("alpha")
        s = control.get_shells("closed")
        assert s.states[0] == ("a1", 1)
        assert s.occupations[0] == 2

    def test_get_charge(self, testdir, delete_tmp_dir):
        c_empty = Control.empty()
        assert c_empty.get_charge() is None
        with temp_dir(delete_tmp_dir):
            shutil.copy2(os.path.join(testdir, 'control', 'control_test-subfiles'), "control")
            c = Control.from_file()
            assert c.get_charge() == 0

    def test_get_subfiles_list(self, testdir):

        c = Control.from_file(os.path.join(testdir, 'control', 'control_test-subfiles'))
        l = c.get_subfiles_list()
        assert set(l) == {"basis", "twoint2", "forceapprox", "mos", "coord", "energy", "gradient", "twoint1"}

    def test_cpc(self, testdir, delete_tmp_dir):
        with temp_dir(delete_tmp_dir):
            c_empty = Control.empty()
            c_empty.cpc("test_dir")
            assert os.path.isfile(os.path.join("test_dir", "control"))

            shutil.copy2(os.path.join(testdir, 'control', 'control_test-subfiles'), "control")
            shutil.copy2(os.path.join(testdir, 'control', 'gradient_test-gradient'), "gradient")

            c = Control.from_file("control")

            c.cpc("test_dir")
            c_test = Control.from_file(os.path.join("test_dir", "control"))
            # check that did not overwrite
            assert len(c_test.dg_list) == 1
            assert os.path.isfile(os.path.join("test_dir", "gradient"))

            c.cpc("test_dir", force_overwrite=True)
            c_test = Control.from_file(os.path.join("test_dir", "control"))
            assert len(c_test.dg_list) > 1


@pytest.mark.parametrize('control_filename', ['control_test-Control'])
class TestDatagroupFunctions(object):

    def test_sdg(self, control_filepath, delete_tmp_dir):

        with temp_dir(delete_tmp_dir):
            shutil.copy2(control_filepath, "control")
            assert sdg("operating system") == " unix\n"

    def test_kdg(self, control_filepath, delete_tmp_dir):

        with temp_dir(delete_tmp_dir):
            shutil.copy2(control_filepath, "control")
            with open("control") as f:
                s = f.read()
                matches = re.findall(r"\$operating system unix",
                                     s, re.MULTILINE)
                assert matches == ["$operating system unix"]

            kdg("operating system", backup_file="control_backup")
            assert os.path.isfile("control_backup")

            with open("control") as f:
                s = f.read()
                matches = re.findall(r"\$operating system unix",
                                     s, re.MULTILINE)
                assert matches == []

            with open("control_backup") as f:
                s = f.read()
                matches = re.findall(r"\$operating system unix",
                                     s, re.MULTILINE)
                assert matches == ["$operating system unix"]

    def test_adg(self, control_filepath, delete_tmp_dir):

        with temp_dir(delete_tmp_dir):
            shutil.copy2(control_filepath, "control")
            with open("control") as f:
                s = f.read()
                matches = re.findall(r"\$scfconv 6", s, re.MULTILINE)
                assert matches == []

            adg("scfconv", "6", backup_file="control_backup")
            assert os.path.isfile("control_backup")

            with open("control") as f:
                s = f.read()
                matches = re.findall(r"\$scfconv 6", s, re.MULTILINE)
                assert matches == ["$scfconv 6"]

            with open("control_backup") as f:
                s = f.read()
                matches = re.findall(r"\$scfconv 6", s, re.MULTILINE)
                assert matches == []

    def test_cdg(self, control_filepath, delete_tmp_dir):

        with temp_dir(delete_tmp_dir):
            shutil.copy2(control_filepath, "control")
            with open("control") as f:
                s = f.read()
                matches = re.findall(r"\$operating system unix",
                                     s, re.MULTILINE)
                assert matches == ["$operating system unix"]

            cdg("operating system", "windows", backup_file="control_backup")
            assert os.path.isfile("control_backup")

            with open("control") as f:
                s = f.read()
                matches = re.findall(r"\$operating system windows",
                                     s, re.MULTILINE)
                assert matches == ["$operating system windows"]
                matches = re.findall(r"\$operating system unix",
                                     s, re.MULTILINE)
                assert matches == []

            with open("control_backup") as f:
                s = f.read()
                matches = re.findall(r"\$operating system unix",
                                     s, re.MULTILINE)
                assert matches == ["$operating system unix"]
                matches = re.findall(r"\$operating system windows",
                                     s, re.MULTILINE)
                assert matches == []

    def test_mdgo(self, control_filepath, delete_tmp_dir):

        with temp_dir(delete_tmp_dir):
            shutil.copy2(control_filepath, "control")
            options = {"internal": None,
                       "redundant": None,
                       "cartesian": "cartesian off",
                       "global": None}
            mdgo("optimize", options, backup_file="control_backup")
            assert os.path.isfile("control_backup")

            with open("control") as f:
                s = f.read()
                matches = re.findall(r"cartesian off",
                                     s, re.MULTILINE)
                assert matches == ["cartesian off"]
                matches = re.findall(r"\$redundant",
                                     s, re.MULTILINE)
                assert matches == []

    def test_sdgo(self, control_filepath, delete_tmp_dir):
        with temp_dir(delete_tmp_dir):
            shutil.copy2(control_filepath, "control")
            assert sdgo("$optimize", "internal") == "   off"

    def test_cpc(self, control_filepath, testdir, delete_tmp_dir):
        with temp_dir(delete_tmp_dir):

            os.makedirs("orig_dir")

            shutil.copy2(os.path.join(testdir, 'control', 'control_test-subfiles'),
                         os.path.join("orig_dir", "control"))
            shutil.copy2(os.path.join(testdir, 'control', 'gradient_test-gradient'),
                         os.path.join("orig_dir", "gradient"))

            cpc("test_dir", control_dir="orig_dir")
            assert os.path.isfile(os.path.join("test_dir", "control"))
            assert os.path.isfile(os.path.join("test_dir", "gradient"))

            shutil.copy2(os.path.join(testdir, 'control', 'control_test-subfiles'), "control")

            cpc("test_dir2")
            assert os.path.isfile(os.path.join("test_dir2", "control"))
