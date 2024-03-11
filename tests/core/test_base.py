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


import os

import pytest
from pymatgen.core.structure import Structure

from turbomoleio.core.base import get_mol_and_indices_frozen
from turbomoleio.core.control import DataGroups
from turbomoleio.core.molecule import MoleculeSystem
from turbomoleio.testing import temp_dir


def test_get_mol_and_indices_frozen_cell_lattice(test_data):
    with pytest.raises(
        ValueError, match=r"Only cell_string or lattice_string should be set."
    ):
        get_mol_and_indices_frozen(
            "fakestring",
            cell_string="fakecellstring",
            lattice_string="fakelatticestring",
        )


@pytest.mark.parametrize("molecule_filename", ["co2.json"])
def test_get_mol_and_indices_frozen_cell(molecule):
    ms = MoleculeSystem(molecule)
    dg = DataGroups(ms.to_coord_string())
    coord_string = dg.show_data_group("coord")
    with pytest.raises(
        ValueError,
        match=r"The \$cell data group should contain 1 number for 1D periodic systems",
    ):
        get_mol_and_indices_frozen(
            coord_string, cell_string="1.0 1.0 1.0", periodic_string="1"
        )
    with pytest.raises(
        ValueError,
        match=r"The \$cell data group should contain 3 numbers for 2D periodic systems",
    ):
        get_mol_and_indices_frozen(coord_string, cell_string="1.0", periodic_string="2")
    with pytest.raises(
        ValueError,
        match=r"The \$cell data group should contain 6 numbers for 3D periodic systems",
    ):
        get_mol_and_indices_frozen(coord_string, cell_string="1.0", periodic_string="3")
    with pytest.raises(RuntimeError, match=r"Periodicity should be one of 1, 2 or 3."):
        get_mol_and_indices_frozen(coord_string, cell_string="1.0", periodic_string="4")

    coord_system = get_mol_and_indices_frozen(
        coord_string, cell_string="10.0", periodic_string="1"
    )
    struct = coord_system.molecule_or_structure
    assert isinstance(struct, Structure)
    assert struct.lattice.abc == pytest.approx([5.29177210903, 5.0, 8.028152])
    assert struct.lattice.is_orthogonal

    coord_system = get_mol_and_indices_frozen(
        coord_string, cell_string="10.0", periodic_string="1", periodic_extension=8.0
    )
    struct = coord_system.molecule_or_structure
    assert struct.lattice.abc == pytest.approx([5.29177210903, 8.0, 11.028152])

    coord_system = get_mol_and_indices_frozen(
        coord_string, cell_string="10.0 12.0 14.0 65 70 80", periodic_string="3"
    )
    struct = coord_system.molecule_or_structure
    assert isinstance(struct, Structure)
    assert struct.lattice.abc == pytest.approx(
        [5.29177210903, 6.350126531, 7.408480953]
    )
    assert struct.lattice.angles == pytest.approx([65, 70, 80])


@pytest.mark.parametrize("molecule_filename", ["co2.json"])
def test_get_mol_and_indices_frozen_lattice(molecule):
    ms = MoleculeSystem(molecule)
    dg = DataGroups(ms.to_coord_string())
    coord_string = dg.show_data_group("coord")
    with pytest.raises(
        ValueError,
        match=r"The lattice_numbers should contain 1 number for 1D periodic systems",
    ):
        get_mol_and_indices_frozen(
            coord_string, lattice_string="1.0 1.0 1.0", periodic_string="1"
        )
    with pytest.raises(
        ValueError,
        match=r"The lattice_numbers should contain 1 number for 1D periodic systems",
    ):
        get_mol_and_indices_frozen(
            coord_string, lattice_string="1.0\n1.0", periodic_string="1"
        )
    with pytest.raises(
        ValueError,
        match=r"The lattice_numbers should contain 2x2 numbers for 2D periodic systems",
    ):
        get_mol_and_indices_frozen(
            coord_string, lattice_string="1.0\n1.0", periodic_string="2"
        )
    with pytest.raises(
        ValueError,
        match=r"The lattice_numbers should contain 3x3 numbers for 3D periodic systems",
    ):
        get_mol_and_indices_frozen(
            coord_string, lattice_string="1.0\n1.0", periodic_string="3"
        )
    with pytest.raises(RuntimeError, match=r"Periodicity should be one of 1, 2 or 3."):
        get_mol_and_indices_frozen(
            coord_string, lattice_string="1.0\n1.0", periodic_string="4"
        )


@pytest.mark.parametrize("molecule_filename", ["co2.json"])
def test_to_file_no_format(molecule, delete_tmp_dir):
    ms = MoleculeSystem(molecule)
    with temp_dir(delete_tmp_dir) as tmp_dir:
        fname = os.path.join(tmp_dir, "co2.xyz")
        ms.to_file(filepath=fname, fmt=None)
        assert os.path.isfile("co2.xyz")
