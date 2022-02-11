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

import pytest

from pymatgen.core.structure import Structure
from turbomoleio.core.periodic import PeriodicSystem
from turbomoleio.core.datagroups import DataGroups


class TestPeriodicSystem:

    @pytest.mark.parametrize('structure_filename', ['SiO2.json'])
    def test_periodic_3d(self, structure):
        ps = PeriodicSystem(structure=structure)
        assert ps.periodicity == 3
        coord_string = ps.to_coord_string()
        assert '$periodic 3' in coord_string
        assert '$cell' in coord_string
        lines = coord_string.split('\n')
        iline_cell = lines.index('$cell')
        abc_alphabetagamma = [float(nn) for nn in lines[iline_cell+1].split()]
        assert abc_alphabetagamma == pytest.approx(
            [9.489264667, 9.489264667, 10.41346809,
             90.0, 90.0, 120.0]
        )

    @pytest.mark.parametrize('structure_filename', ['SiO2_rotated.json'])
    def test_periodic_wrong_rotation(self, structure):
        ps = PeriodicSystem(structure=structure)
        with pytest.raises(ValueError,
                           match=r'Lattice should be oriented such that first vector '
                                 r'is aligned with x cartesian direction and second vector '
                                 r'is in the xy cartesian plane.'):
            ps.to_coord_string()

    @pytest.mark.parametrize('structure_filename', ['graphene.json'])
    def test_periodic_2d(self, structure):
        ps = PeriodicSystem(structure=structure, periodicity=2)
        assert ps.periodicity == 2
        coord_string = ps.to_coord_string()
        assert '$periodic 2' in coord_string
        assert '$cell' in coord_string
        lines = coord_string.split('\n')
        iline_cell = lines.index('$cell')
        ab_gamma = [float(nn) for nn in lines[iline_cell + 1].split()]
        assert ab_gamma == pytest.approx(
            [4.664630258, 4.664630258, 120.0]
        )

    @pytest.mark.parametrize('structure_filename', ['graphene'])
    def test_read_periodic_2d(self, structure_filepath):
        ps = PeriodicSystem.from_file(structure_filepath, fmt='coord', periodic_extension=5.0)
        assert isinstance(ps.structure, Structure)
        assert ps.periodicity == 2
        current_coord_dg = DataGroups(string=ps.to_coord_string())
        ref_coord_dg = DataGroups.from_file(structure_filepath)
        compare_coord = current_coord_dg.compare(ref_coord_dg, tol=1e-6)
        assert compare_coord is None, compare_coord
        assert ps.structure.lattice.c == pytest.approx(5.0)
        ps = PeriodicSystem.from_file(structure_filepath, fmt='coord', periodic_extension=7.0)
        assert ps.structure.lattice.c == pytest.approx(7.0)
