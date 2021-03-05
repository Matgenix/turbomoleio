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
Tests for the MoleculeSystem object.
Mainly testing that int_def and user-defined bonds are taken into account correctly.
Relies on DefineRunner.
"""
import pytest
import shutil
import re

from turbomoleio.testfiles.utils import temp_dir, assert_MSONable
from turbomoleio.core.molecule import Distance, BondAngle, DihedralAngle, InverseDistance, OutOfPlaneAngle
from turbomoleio.core.molecule import CollinearBendingAngle, PerpendicularBendingAngle, InternalDefinition
from turbomoleio.core.molecule import MoleculeSystem
from turbomoleio.core.datagroups import DataGroups
from turbomoleio.input.define import DefineRunner


dr_parameters = dict(ired=True)

@pytest.mark.integration
@pytest.mark.parametrize('control_filename', ['control_itest-molecule'])
@pytest.mark.parametrize('molecule_filename', ['ch4.json'])
class TestMoleculeSystem:

    def test_distance(self, molecule, control_filepath, delete_tmp_dir):
        ms = MoleculeSystem(molecule, frozen_indices={2,3})
        ms.add_distance(0, 1, weights=1.0, value=2.5)
        assert ms.has_inconsistencies()

        with temp_dir(delete_tmp_dir) as tmp_dir:
            ms.to_file("coord", fmt="coord")
            shutil.copy2(control_filepath, "control")
            dr = DefineRunner(parameters=dr_parameters)
            assert dr.run_update_internal_coords()

            ms_new = MoleculeSystem.from_file("coord")
            assert not ms_new.has_inconsistencies()

            # check that intdef is present and that the value has been updated
            assert len(ms_new.int_def) == 1
            assert ms_new.int_def[0].value == pytest.approx(ms.int_def[0].value, abs=1e-4)

            # check that there are the same user-defined bonds
            assert ms.user_defined_bonds == ms_new.user_defined_bonds

            # check that in redundant the coordinate has been taken into account
            dg = DataGroups.from_file("coord")
            redundant = dg.sdg("redundant", strict=True)
            assert re.search(r"\d+\s+f\s+1\.0+\s+stre\s+1\s+2\s+val=\s*2.5", redundant)

            # check the presence of the same frozen atoms
            assert ms.frozen_indices == ms_new.frozen_indices

            assert str(ms_new.int_def[0])

            assert_MSONable(ms_new)

    def test_distance_no_val(self, molecule, control_filepath, delete_tmp_dir):
        ms = MoleculeSystem(molecule)
        ms.add_distance(0, 1, weights=1.0, value=None)
        assert not ms.has_inconsistencies()

        with temp_dir(delete_tmp_dir) as tmp_dir:
            ms.to_file("coord", fmt="coord")
            shutil.copy2(control_filepath, "control")
            dr = DefineRunner(parameters=dr_parameters)
            assert dr.run_update_internal_coords()

            ms_new = MoleculeSystem.from_file("coord")
            assert not ms_new.has_inconsistencies()

            # check that intdef is present and that the value has been updated
            assert len(ms_new.int_def) == 1
            assert ms_new.int_def[0].value == pytest.approx(2.05846, abs=1e-4)

            # check that in redundant the coordinate has been taken into account
            dg = DataGroups.from_file("coord")
            redundant = dg.sdg("redundant", strict=True)
            assert re.search(r"\d+\s+f\s+1\.0+\s+stre\s+1\s+2\s+val=\s*2\.05", redundant)

    def test_distance_linear_combination(self, molecule, control_filepath, delete_tmp_dir):
        ms = MoleculeSystem(molecule)
        ms.add_distance([0, 1], [1, 2], weights=[1.0, 1.0], value=None)
        assert not ms.has_inconsistencies()

        with temp_dir(delete_tmp_dir) as tmp_dir:
            ms.to_file("coord", fmt="coord")
            shutil.copy2(control_filepath, "control")
            dr = DefineRunner(parameters=dr_parameters)
            assert dr.run_update_internal_coords()

            ms_new = MoleculeSystem.from_file("coord")
            assert not ms_new.has_inconsistencies()

            # check that intdef is present and that the value has been updated
            assert len(ms_new.int_def) == 1
            assert ms_new.int_def[0].value == pytest.approx(2.70995, abs=1e-4)

            # check that in redundant the coordinate has been taken into account
            dg = DataGroups.from_file("coord")
            redundant = dg.sdg("redundant", strict=True)
            assert re.search(r"\d+\s+f\s+1\.0+\s+stre\s+1\s+2\s+val=\s*2\.70\d+\s+1\.0+\s+stre\s+2\s+3\s*$",
                             redundant, re.DOTALL|re.MULTILINE)

    def test_bond_angle(self, molecule, control_filepath, delete_tmp_dir):
        ms = MoleculeSystem(molecule)
        ms.add_bond_angle(0, 1, 2, value=30, weights=[1.0])
        assert ms.has_inconsistencies()

        with temp_dir(delete_tmp_dir) as tmp_dir:
            ms.to_file("coord", fmt="coord")
            shutil.copy2(control_filepath, "control")
            dr = DefineRunner(parameters=dr_parameters)
            assert dr.run_update_internal_coords()

            ms_new = MoleculeSystem.from_file("coord")
            assert not ms_new.has_inconsistencies()

            # check that intdef is present and that the value has been updated
            assert len(ms_new.int_def) == 1
            assert ms_new.int_def[0].value == pytest.approx(ms.int_def[0].value, abs=1e-4)

            # check that there are the same user-defined bonds
            assert ms.user_defined_bonds == ms_new.user_defined_bonds

            # check that in redundant the coordinate has been taken into account
            dg = DataGroups.from_file("coord")
            redundant = dg.sdg("redundant", strict=True)
            assert re.search(r"\d+\s+f\s+1\.0+\s+bend\s+1\s+2\s+3\s+val=\s*30", redundant)

            assert str(ms_new.int_def[0])

    def test_dihedral_angle(self, molecule, control_filepath, delete_tmp_dir):
        ms = MoleculeSystem(molecule)
        ms.add_dihedral(0, 1, 2, 3, value=30, weights=1.0)
        assert ms.has_inconsistencies()

        with temp_dir(delete_tmp_dir) as tmp_dir:
            ms.to_file("coord", fmt="coord")
            shutil.copy2(control_filepath, "control")
            dr = DefineRunner(parameters=dr_parameters)
            assert dr.run_update_internal_coords()

            ms_new = MoleculeSystem.from_file("coord")
            assert not ms_new.has_inconsistencies()

            # check that intdef is present and that the value has been updated
            assert len(ms_new.int_def) == 1
            assert ms_new.int_def[0].value == pytest.approx(ms.int_def[0].value, abs=1e-4)

            # check that there are the same user-defined bonds
            assert ms.user_defined_bonds == ms_new.user_defined_bonds

            # check that in redundant the coordinate has been taken into account
            dg = DataGroups.from_file("coord")
            redundant = dg.sdg("redundant", strict=True)
            assert re.search(r"\d+\s+f\s+1\.0+\s+tors\s+1\s+2\s+3\s+4\s+val=\s*30", redundant)

            assert str(ms_new.int_def[0])

    def test_inverse_distance(self, molecule, control_filepath, delete_tmp_dir):
        ms = MoleculeSystem(molecule)
        ms.int_def.append(InverseDistance(indices=[0, 1], value=0.5, weights=[1.0], status="f"))
        assert ms.has_inconsistencies()

        with temp_dir(delete_tmp_dir) as tmp_dir:
            ms.to_file("coord", fmt="coord")
            shutil.copy2(control_filepath, "control")
            dr = DefineRunner(parameters=dr_parameters)
            assert dr.run_update_internal_coords()

            ms_new = MoleculeSystem.from_file("coord")
            assert not ms_new.has_inconsistencies()

            # check that intdef is present and that the value has been updated
            assert len(ms_new.int_def) == 1
            assert ms_new.int_def[0].value == pytest.approx(ms.int_def[0].value, abs=1e-4)

            # check that there are the same user-defined bonds
            assert ms.user_defined_bonds == ms_new.user_defined_bonds

            # check that in redundant the coordinate has been taken into account
            dg = DataGroups.from_file("coord")
            redundant = dg.sdg("redundant", strict=True)
            assert re.search(r"\d+\s+f\s+1\.0+\s+invr\s+1\s+2\s+val=\s*0.5", redundant)

            assert str(ms_new.int_def[0])

    def test_out_of_plane(self, molecule, control_filepath, delete_tmp_dir):
        ms = MoleculeSystem(molecule)
        ms.int_def.append(OutOfPlaneAngle(indices=[[0, 1, 2, 3]], value=-25, weights=None, status="f"))
        assert ms.has_inconsistencies()

        with temp_dir(delete_tmp_dir) as tmp_dir:
            ms.to_file("coord", fmt="coord")
            shutil.copy2(control_filepath, "control")
            dr = DefineRunner(parameters=dr_parameters)
            assert dr.run_update_internal_coords()

            ms_new = MoleculeSystem.from_file("coord")
            assert not ms_new.has_inconsistencies()

            # check that intdef is present and that the value has been updated
            assert len(ms_new.int_def) == 1
            assert ms_new.int_def[0].value == pytest.approx(ms.int_def[0].value, abs=1e-4)

            # check that there are the same user-defined bonds
            assert ms.user_defined_bonds == ms_new.user_defined_bonds

            # check that in redundant the coordinate has been taken into account
            dg = DataGroups.from_file("coord")
            redundant = dg.sdg("redundant", strict=True)
            assert re.search(r"\d+\s+f\s+1\.0+\s+outp\s+1\s+2\s+3\s+4\s+val=\s*-25", redundant)

            assert str(ms_new.int_def[0])

    def test_linc(self, molecule, control_filepath, delete_tmp_dir):
        ms = MoleculeSystem(molecule)
        ms.int_def.append(CollinearBendingAngle(indices=[[0, 1, 2, 3]], value=80, status="f"))
        assert ms.has_inconsistencies()

        with temp_dir(delete_tmp_dir) as tmp_dir:
            ms.to_file("coord", fmt="coord")
            shutil.copy2(control_filepath, "control")
            dr = DefineRunner(parameters=dr_parameters)
            assert dr.run_update_internal_coords()

            ms_new = MoleculeSystem.from_file("coord")
            assert not ms_new.has_inconsistencies()

            # check that intdef is present and that the value has been updated
            assert len(ms_new.int_def) == 1
            assert ms_new.int_def[0].value == pytest.approx(ms.int_def[0].value, abs=1e-4)

            # check that there are the same user-defined bonds
            assert ms.user_defined_bonds == ms_new.user_defined_bonds

            # check that in redundant the coordinate has been taken into account
            dg = DataGroups.from_file("coord")
            redundant = dg.sdg("redundant", strict=True)
            assert re.search(r"\d+\s+f\s+1\.0+\s+linc\s+1\s+2\s+3\s+4\s+val=\s*80", redundant)

            assert str(ms_new.int_def[0])

    def test_linp(self, molecule, control_filepath, delete_tmp_dir):
        ms = MoleculeSystem(molecule)
        ms.int_def.append(PerpendicularBendingAngle(indices=[[0, 1, 2, 3]], value=20, status="f"))
        assert ms.has_inconsistencies()

        with temp_dir(delete_tmp_dir) as tmp_dir:
            ms.to_file("coord", fmt="coord")
            shutil.copy2(control_filepath, "control")
            dr = DefineRunner(parameters=dr_parameters)
            assert dr.run_update_internal_coords()

            ms_new = MoleculeSystem.from_file("coord")
            assert not ms_new.has_inconsistencies()

            # check that intdef is present and that the value has been updated
            assert len(ms_new.int_def) == 1
            assert ms_new.int_def[0].value == pytest.approx(ms.int_def[0].value, abs=1e-4)

            # check that there are the same user-defined bonds
            assert ms.user_defined_bonds == ms_new.user_defined_bonds

            # check that in redundant the coordinate has been taken into account
            dg = DataGroups.from_file("coord")
            redundant = dg.sdg("redundant", strict=True)
            assert re.search(r"\d+\s+f\s+1\.0+\s+linp\s+1\s+2\s+3\s+4\s+val=\s*20", redundant)

            assert str(ms_new.int_def[0])
