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
import pytest
import re

from pymatgen.core.periodic_table import DummySpecie
from pymatgen.core.structure import Molecule
from pymatgen.core.structure import Structure
from turbomoleio.testfiles.utils import temp_dir, assert_MSONable
from turbomoleio.core.molecule import Distance, BondAngle, DihedralAngle, InverseDistance, OutOfPlaneAngle
from turbomoleio.core.molecule import CollinearBendingAngle, PerpendicularBendingAngle, InternalDefinition
from turbomoleio.core.molecule import MoleculeSystem
from turbomoleio.core.datagroups import DataGroups

try:
    from pymatgen.core.periodic_table import DummySpecies
except ImportError:
    DummySpecies = None


replace_comments_re = re.compile(r"^\s*#.*?$", flags=re.MULTILINE)


def check_dg(to_test, ref):
    """
    Helper function to check if the string of $intdef or $coord datagroup is
    equivalent to a reference.
    Treats all the numbers as floats and compares them with approx.
    Comments are removed.
    In case of $intedef assumes that the value is given with the string "val=",
    with no spaces.

    Exceptions raised through assertions if comparison fails.

    Args:
        to_test (str): the string to be tested.
        ref (str): the reference string.
    """
    # remove the val= to avoid problems if spaces are left before the actual value
    to_test = to_test.replace("val=", "")
    ref = ref.replace("val=", "")

    # comments (empty lines are fine for the rest)
    to_test = replace_comments_re.sub("", to_test)
    ref = replace_comments_re.sub("", ref)

    to_test = to_test.strip()
    ref = ref.strip()

    t_split = to_test.split()
    r_split = ref.split()
    assert len(t_split) == len(r_split)
    for s, r in zip(t_split, r_split):
        # if a float/int compare the value, otherwise they should be the same string
        try:
            s_f = float(s)
            r_f = float(r)
            assert s_f == pytest.approx(r_f)
        except ValueError:
            assert s == r


def check_user_defined_bonds_dg(to_test, ref):
    """
    Helper function to check if the string of $user-defined bonds datagroup is
    equivalent to a reference.
    Comments are removed.
    Assumes that the elements are properly sorted and compare strings removing
    spaces.

    Exceptions raised through assertions if comparison fails.

    Args:
        to_test (str): the string to be tested.
        ref (str): the reference string.
    """

    # comments (empty lines are fine for the rest)
    to_test = replace_comments_re.sub("", to_test)
    ref = replace_comments_re.sub("", ref)

    to_test = re.sub(r"\s", "", to_test)
    ref = re.sub(r"\s", "", ref)

    assert to_test == ref


@pytest.mark.parametrize('molecule_filename', ['ch4.json'])
class TestInternalDefinition:

    def test_distance(self, molecule):
        d = Distance(status="k", indices=[0,1])
        assert d.is_valid(molecule)
        d.value = 1.0
        assert not d.is_valid(molecule)
        d.value = 2.058455873547
        assert d.is_valid(molecule)

        d = Distance(status="k", indices=[[0,1], [1,2]], weights=[1, 0.5], value=2.49279)
        assert d.is_valid(molecule)

        d = Distance(status="f", indices=[[0,1], [1,2]], weights=[1, -0.5], value=0.25182)
        assert d.is_valid(molecule)
        # wrong status
        with pytest.raises(ValueError, match="Unsupported internal definition status.*z"):
            Distance(status="z", indices=[0,1])

        assert str(d)

        s_out = """1 f 1 stre 1 2 val=0.25182
    -0.5 stre 2 3
"""
        check_dg(d.to_string(), s_out)

    def test_bond_angle(self, molecule):
        # wrong number of indices
        with pytest.raises(ValueError, match="Wrong number of indices.*2"):
            BondAngle(status="k", indices=[0,1])

        ba = BondAngle(status="k", indices=[0,1,2])
        assert ba.is_valid(molecule)
        ba.value = 90
        assert not ba.is_valid(molecule)
        ba.value = 35.26438972
        assert ba.is_valid(molecule)

        ba = BondAngle(status="k", indices=[[0,1,2], [1,2,3]], weights=[1, 0.2], value=39.38699)
        assert ba.is_valid(molecule)

        assert str(ba)

        s_out = """1 k 1 bend 1 2 3 val=39.38699
    0.2 bend 2 3 4
"""
        check_dg(ba.to_string(), s_out)

    def test_dihedral_angle(self, molecule):
        # wrong number of indices
        with pytest.raises(ValueError, match="Wrong number of indices.*3"):
            DihedralAngle(status="k", indices=[0,1,2])

        da = DihedralAngle(status="k", indices=[0,1,2,3])
        assert da.is_valid(molecule)
        da.value = 90
        assert not da.is_valid(molecule)
        da.value = 35.26438972
        assert da.is_valid(molecule)

        da = DihedralAngle(status="k", indices=[[0,1,2,3], [1,2,3,4], [0,2,3,4]],
                           weights=[-0.1, 0.7, 0.7], value=-51.72110)
        assert da.is_valid(molecule)

        assert str(da)

        s_out = """1 k -0.1 tors 1 2 3 4 val=-51.7211
    0.7 tors 2 3 4 5
    0.7 tors 1 3 4 5
"""
        check_dg(da.to_string(), s_out)

    def test_inverse_distance(self, molecule):
        i = InverseDistance(status="k", indices=[[0,1]])
        assert i.is_valid(molecule)
        i.value = 1.0
        assert not i.is_valid(molecule)
        i.value = 0.48580103797
        assert i.is_valid(molecule)

        assert str(i)

        s_out = """1 k 1.0 invr 1 2 val=0.48580103797"""
        check_dg(i.to_string(), s_out)

    def test_out_of_plane_angle(self, molecule):
        outp = OutOfPlaneAngle(status="k", indices=[0,1,2,3])
        assert outp.is_valid(molecule)
        outp.value = 1.0
        assert not outp.is_valid(molecule)
        outp.value = -19.47122076
        assert outp.is_valid(molecule)

        assert str(outp)

        s_out = """1 k 1.0 outp 1 2 3 4 val=-19.47122076"""
        check_dg(outp.to_string(), s_out)

    def test_collinear_bending_angle(self, molecule):
        linc = CollinearBendingAngle(status="k", indices=[0,1,2,3])
        assert linc.is_valid(molecule)
        linc.value = 1.0
        assert not linc.is_valid(molecule)
        linc.value = 84.73561
        assert linc.is_valid(molecule)

        assert str(linc)

        s_out = """1 k 1.0 linc 1 2 3 4 val=84.73561"""
        check_dg(linc.to_string(), s_out)

    def test_perpendicular_bending_angle(self, molecule):
        linp = PerpendicularBendingAngle(status="k", indices=[0,1,2,3])
        # is_valid always returns True
        assert linp.is_valid(molecule)
        linp.value = 1.0
        assert not linp.is_valid(molecule)
        linp.value = 24.73561
        assert linp.is_valid(molecule)

        assert str(linp)

        s_out = """1 k 1.0 linp 1 2 3 4 val=24.73561"""
        check_dg(linp.to_string(), s_out)


    def test_from_string(self, molecule):
        # Distance
        s_distance = """1 k  1.0000000000000 stre    1    2           val=   0.25182
       -0.5000000000000 stre    2    3
"""
        i = InternalDefinition.from_string(s_distance)
        assert isinstance(i, Distance)
        assert i.is_valid(molecule)
        assert len(i.indices) == 2
        assert i.value == pytest.approx(0.25182)
        assert i.status == "k"

        # Bond angle
        s_bond_angle = """   1 f  1.0000000000000 bend    1    2    3      val=  19.38699
       -0.2000000000000 bend    2    3    4
"""
        i = InternalDefinition.from_string(s_bond_angle)
        assert isinstance(i, BondAngle)
        assert i.is_valid(molecule)
        assert len(i.indices) == 2
        assert i.value == pytest.approx(19.38699)
        assert i.status == "f"

        # Dihedral angle
        s_dihedral_angle = """   1 i -0.1000000000000 tors    1    2    3    4 val= -51.72110
        0.7000000000000 tors    2    3    4    5
        0.7000000000000 tors    1    3    4    5
"""
        i = InternalDefinition.from_string(s_dihedral_angle)
        assert isinstance(i, DihedralAngle)
        assert i.is_valid(molecule)
        assert len(i.indices) == 3
        assert i.value == pytest.approx(-51.72110)
        assert i.status == "i"

        # Inverse distance
        s_inverse_distance = """# definitions of internal coordinates
  1 k  1.0000000000000 invr    1    2 
"""
        i = InternalDefinition.from_string(s_inverse_distance)
        assert isinstance(i, InverseDistance)
        assert i.is_valid(molecule)
        assert len(i.indices) == 1
        assert i.value == None
        assert i.status == "k"

        # Out of plane angle
        s_out_of_plane = """1 d  1.0000000000000 outp    1    2    3    4 val= -19.47122
"""
        i = InternalDefinition.from_string(s_out_of_plane)
        assert isinstance(i, OutOfPlaneAngle)
        assert i.is_valid(molecule)
        assert len(i.indices) == 1
        assert i.value == pytest.approx(-19.47122)
        assert i.status == "d"

        # collinear bending angle
        s_linc = """   1 k  1.0000000000000 linc    1    2    3    4 val=  84.73561
"""
        i = InternalDefinition.from_string(s_linc)
        assert isinstance(i, CollinearBendingAngle)
        assert i.is_valid(molecule)
        assert len(i.indices) == 1
        assert i.value == pytest.approx(84.73561)
        assert i.status == "k"

        # perpendicular bending angle
        s_linp = """   1 k  1.0000000000000 linp    1    2    3    4 val=  24.73561
"""
        i = InternalDefinition.from_string(s_linp)
        assert isinstance(i, PerpendicularBendingAngle)
        assert i.is_valid(molecule)
        assert len(i.indices) == 1
        assert i.value == pytest.approx(24.73561)
        assert i.status == "k"

        # check some possible errors
        # wrong number of atoms
        s = """   1 f  1.0000000000000 bend    1    2      val=  19.38699
"""
        with pytest.raises(ValueError, match="invalid literal for int.*val=.*"):
            InternalDefinition.from_string(s)

        # wrong type
        s = """   1 z  1.0000000000000 bend    1    2  3    val=  19.38699
        """
        with pytest.raises(ValueError, match="Unsupported internal definition status.*z"):
            InternalDefinition.from_string(s)

        # wrong type of coordinate
        s = """   1 f  1.0000000000000 zzzz    1    2  3    val=  19.38699
        """
        with pytest.raises(ValueError, match="Could not find a subclass with coord_str matching zzzz"):
            InternalDefinition.from_string(s)

        # misformatted string
        s = """    1.0000000000000 bend    1    2  3    val=  19.38699
        """
        with pytest.raises(ValueError, match="Could not find a subclass with coord_str matching 2"):
            InternalDefinition.from_string(s)


class TestMoleculeSystem:

    def test_from_string(self):
        # basic test
        string = """
$coord
 .00000000000000       .00000000000000       .00000000000000      n       
-1.15103063747470     -1.99364354517457       .00000000000000      o       
2.30206127494940       .00000000000000       .00000000000000      o       
-1.15103063747470      1.99364354517457       .00000000000000      o       
$end        
"""

        ms = MoleculeSystem.from_string(string=string, fmt="coord")

        mol = ms.molecule
        assert mol[1].coords[0] == pytest.approx(-0.6090991821345737)
        assert len(mol) == 4

        assert len(ms.frozen_indices) == 0

        assert_MSONable(ms)

        # no coord
        with pytest.raises(ValueError, match=r'^The string does not contain \$coord!$'):
            MoleculeSystem.from_string(string="$end", fmt="coord")

        # with frozen and internal definitions
        string = """
$coord
 .00000000000000       .00000000000000       .00000000000000      n  f
-1.15103063747470     -1.99364354517457       .00000000000000      o       
2.30206127494940       .00000000000000       .00000000000000      o       
-1.15103063747470      1.99364354517457       .00000000000000      o f       
$intdef
# definitions of internal coordinates

1 k  1.0000000000000 stre    1    2           val=   2.43987
2 f  1.0000000000000 bend    1    2  3        
     -0.5000000000000 bend    2    3  4        
$end        
"""
        ms = MoleculeSystem.from_string(string=string, fmt="coord")
        mol = ms.molecule
        assert mol[1].coords[0] == pytest.approx(-0.6090991821345737)
        assert len(mol) == 4

        assert ms.frozen_indices == {0, 3}

        assert len(ms.int_def) == 2
        assert ms.int_def[0].value == pytest.approx(2.43987)
        assert ms.int_def[1].value is None
        assert ms.int_def[0].status == "k"
        assert ms.int_def[1].status == "f"
        assert ms.int_def[0].indices[0] == [0, 1]
        assert len(ms.int_def[1].indices) == 2
        assert ms.int_def[1].weights[1] == pytest.approx(-0.5)

        dg = DataGroups(ms.to_coord_string())
        dg_ref = DataGroups(string)
        assert len(dg.dg_list) == 3
        check_dg(dg.sdg("coord", strict=True), dg_ref.sdg("coord", strict=True))
        check_dg(dg.sdg("intdef", strict=True), dg_ref.sdg("intdef", strict=True))

        # with user-defined bonds
        string = """
$coord
 .00000000000000       .00000000000000       .00000000000000      n 
-1.15103063747470     -1.99364354517457       .00000000000000      o       
2.30206127494940       .00000000000000       .00000000000000      o       
-1.15103063747470      1.99364354517457       .00000000000000      o       
$user-defined bonds
1-2, 2 - 3,3|4      
$end        
"""
        ms = MoleculeSystem.from_string(string=string, fmt="coord")
        mol = ms.molecule
        assert mol[1].coords[0] == pytest.approx(-0.6090991821345737)
        assert len(mol) == 4

        assert ms.user_defined_bonds == {(0, "-", 1), (1, "-", 2), (2, "|", 3)}
        dg = DataGroups(ms.to_coord_string())
        dg_ref = DataGroups(string)
        assert len(dg.dg_list) == 3
        check_dg(dg.sdg("coord", strict=True), dg_ref.sdg("coord", strict=True))
        check_user_defined_bonds_dg(dg.sdg("user-defined bonds", strict=True),
                                    dg_ref.sdg("user-defined bonds", strict=True))

        # malformed user-defined bonds
        string = """
$coord
 .00000000000000       .00000000000000       .00000000000000      n 
-1.15103063747470     -1.99364354517457       .00000000000000      o       
2.30206127494940       .00000000000000       .00000000000000      o       
-1.15103063747470      1.99364354517457       .00000000000000      o       
$user-defined bonds
1-2, 2 3,3|4      
$end        
"""
        with pytest.raises(ValueError, match="Cannot parse user-defined bonds.*"):
            MoleculeSystem.from_string(string=string, fmt="coord")

        # from xyz format
        ms = MoleculeSystem.from_string(mol.to(fmt="xyz"), fmt="xyz")
        assert ms.molecule[1].coords[0] == pytest.approx(-0.6090991821345737)

    @pytest.mark.parametrize('molecule_filename', ['co2.json'])
    def test_to_coord_string(self, molecule):
        ms = MoleculeSystem(molecule)
        test_value = """
0.00000000000000 0.00000000000000 0.00000000000000 c
0.00000000000000 0.00000000000000 2.86118897312869 o
0.00000000000000 0.00000000000000 -2.86118897312869 o
"""
        dg = DataGroups(ms.to_coord_string())
        assert len(dg.dg_list) == 2
        check_dg(dg.sdg("coord", strict=True), test_value)

        ms.frozen_indices = {0, 1}
        test_value = """
0.00000000000000 0.00000000000000 0.00000000000000 c f
0.00000000000000 0.00000000000000 2.86118897312869 o f
0.00000000000000 0.00000000000000 -2.86118897312869 o
"""
        dg = DataGroups(ms.to_coord_string())
        assert len(dg.dg_list) == 2
        check_dg(dg.sdg("coord", strict=True), test_value)

    @pytest.mark.parametrize('molecule_filename', ['h2o_dummy_atom'])
    def test_dummy_atoms(self, molecule_filepath):
        ms = MoleculeSystem.from_file(molecule_filepath, fmt="coord")
        mol = ms.molecule
        # Pymatgen's Specie and DummySpecie have been changed to Species and
        # DummySpecies in v2020.10.9. We keep testing both for backward compatibility.
        assert isinstance(mol[-1].specie, (DummySpecies, DummySpecie))
        assert mol[-1].specie.symbol == "Q"

        test_value2 = """
0.00000000000000 0.00000000000000 -0.12178983933899 o
1.41713420892173 0.00000000000000 0.96657854674257 h
-1.41713420892173 0.00000000000000 0.96657854674257 h
0.00000000000000 0.00000000000000 0.00000000000000 q
"""
        dg = DataGroups(ms.to_coord_string())
        assert len(dg.dg_list) == 2
        check_dg(dg.sdg("coord", strict=True), test_value2)

        assert_MSONable(ms)

    @pytest.mark.parametrize('molecule_filename', ['disordered_mol.json'])
    def test_disordered(self, molecule_filepath):
        m = Molecule.from_file(molecule_filepath)
        assert not m.is_ordered

        with pytest.raises(ValueError, match=r'^Turbomoleio and turbomole do not handle disordered structures.$'):
            MoleculeSystem(m).to_coord_string()

    @pytest.mark.parametrize('molecule_filename', ['co2.json'])
    def test_to_file(self, molecule, delete_tmp_dir):
        ms = MoleculeSystem(molecule)

        with temp_dir(delete_tmp_dir) as tmp_dir:
            fname = os.path.join(tmp_dir, 'coord_test')
            ms.to_file(filepath=fname, fmt="coord")
            assert os.path.isfile("coord_test")
            dg = DataGroups.from_file("coord_test")
            assert dg.show_data_group("coord") is not None

            ms.to_file("mol_test.xyz")
            assert os.path.isfile("mol_test.xyz")
            assert ms.from_file("mol_test.xyz") is not None

    @pytest.mark.parametrize('molecule_filename', ['ch4.json'])
    def test_add_distance(self, molecule):
        ms = MoleculeSystem(molecule)
        ms.add_distance(1, 2)
        assert not ms.has_inconsistencies()
        assert len(ms.user_defined_bonds) == 1
        ms.add_distance(0, 1, status="k", add_user_def_bonds=False, value=2.058455873547)
        assert not ms.has_inconsistencies()
        assert len(ms.user_defined_bonds) == 1

        ms.add_distance(2, 3, value=10, weights=1.0)
        assert ms.has_inconsistencies()

        ms.add_distance([0, 1], [2, 3])
        assert len(ms.int_def[-1].indices) == 2

        with pytest.raises(ValueError):
            ms.add_distance([1, 2, 3], [1, 2])

        assert_MSONable(ms)

    @pytest.mark.parametrize('molecule_filename', ['ch4.json'])
    def test_add_bond_angle(self, molecule):
        ms = MoleculeSystem(molecule)
        ms.add_bond_angle(1, 2, 3)
        assert not ms.has_inconsistencies()
        assert len(ms.user_defined_bonds) == 2
        ms.add_bond_angle(0, 1, 2, status="k", add_user_def_bonds=False, value=35.26438972)
        assert not ms.has_inconsistencies()
        assert len(ms.user_defined_bonds) == 2

        ms.add_bond_angle(2, 3, 4, value=10, weights=1.0)
        assert ms.has_inconsistencies()

        ms.add_bond_angle([0,1], [2,3], [3,4])
        assert len(ms.int_def[-1].indices) == 2

        with pytest.raises(ValueError):
            ms.add_bond_angle([1,2,3], [1,2], [2,3])

    @pytest.mark.parametrize('molecule_filename', ['ch4.json'])
    def test_add_dihedral(self, molecule):
        ms = MoleculeSystem(molecule)
        ms.add_dihedral(1, 2, 3, 4)
        assert not ms.has_inconsistencies()
        assert len(ms.user_defined_bonds) == 3
        ms.add_dihedral(0, 1, 2, 3, status="k", add_user_def_bonds=False, value=35.26438972)
        assert not ms.has_inconsistencies()
        assert len(ms.user_defined_bonds) == 3

        ms.add_dihedral(2, 3, 4, 1, value=10)
        assert ms.has_inconsistencies()

        ms.add_dihedral([0, 1], [2, 3], [3, 4], [1, 0])
        assert len(ms.int_def[-1].indices) == 2

        with pytest.raises(ValueError):
            ms.add_dihedral([1, 2, 3], [1, 2], [2, 3], [3, 4])

    @pytest.mark.parametrize('molecule_filename', ['ch4.json'])
    def test_check_index(self, molecule):
        ms = MoleculeSystem(molecule)
        with pytest.raises(ValueError, match="One of the indices representing the atoms is negative or larger then the number of sites"):
            ms._check_index([5])

        with pytest.raises(ValueError, match="One of the indices representing the atoms is negative or larger then the number of sites"):
            ms._check_index([-1])

        assert ms._check_index([4]) is None
