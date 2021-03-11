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
from turbomoleio.input.utils import get_define_template, validate_parameters


def test_get_define_template():
    """Testing get define template util function."""

    dscf_dict = get_define_template("dscf")
    assert dscf_dict["basis"] == "def-SV(P)"

    with pytest.raises(ValueError,
                       match=r'^Could not find template file '
                             r'\S*non_existing_template.yaml$'):
        get_define_template("non_existing_template.yaml")


def test_validate_parameters():

    assert validate_parameters({})

    assert validate_parameters(get_define_template("dscf"))
    assert validate_parameters(get_define_template("ridft"))
    assert validate_parameters(get_define_template("dscf_escf"))
    assert validate_parameters(get_define_template("ridft_escf"))
    assert validate_parameters(get_define_template("ridft_rimp2"))

    assert not validate_parameters({"fake_key": 1})

    # check dependency
    d = {"use_f12*": True}
    assert not validate_parameters(d)
    d["use_f12"] = True
    assert not validate_parameters(d)
    d["method"] = "mp2"
    assert validate_parameters(d)

    #wrong type
    assert not validate_parameters({"method": 1})

    assert validate_parameters({"disp": "DFT-D1"})
    # value not in possible options
    assert not validate_parameters({"disp": "wrong_value"})
