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

import json
import os
import shutil

import numpy as np
import pytest
from monty.json import MSONable
from monty.os import makedirs_p
from monty.serialization import MontyEncoder
from pymatgen.core.structure import Molecule

from turbomoleio.testing import (
    ARRAYS_DIFFER,
    DICT_DIFFERENT_KEYS,
    NUMBERS_DIFFER,
    OBJECTS_DIFFER,
    REF_DICT_TEST_OTHER,
    REF_NUMBER_TEST_OTHER,
    REF_SEQUENCE_TEST_OTHER,
    REF_STRING_TEST_OTHER,
    SEQUENCE_DIFFERENT_SIZES,
    STRINGS_DIFFER,
    TEST_NUMBER_REF_OTHER,
    TEST_STRING_REF_OTHER,
    assert_almost_equal,
    assert_MSONable,
    compare_differences,
    get_control_integration,
    get_sp,
    get_test_data_dir,
    get_tfp,
    gisnan,
    temp_dir,
    touch_file,
)


class MSONableExample(MSONable):
    def __init__(self, a, b):
        self.a = a
        self.b = b


class ExplicitAsFromDictExample:
    def __init__(self, d, e):
        self.d = d
        self.e = e

    def as_dict(self):
        return {
            "@module": "tests.test_testing",
            "@class": "ExplicitAsFromDictExample",
            "@version": None,
            "d": self.d,
            "e": self.e,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(d=d["d"], e=d["e"])

    def to_json(self) -> str:
        return json.dumps(self, cls=MontyEncoder)


class TestFunctions(object):
    def test_temp_dir(self):
        # test that the file is not deleted if delete=False
        # manually delete the file after testing
        with temp_dir(delete=False) as tmp_dir:
            assert tmp_dir == os.getcwd()
            shutil.rmtree(tmp_dir)

        # test without changing to the directory
        with temp_dir(delete=True, changedir=False) as tmp_dir:
            assert os.path.exists(tmp_dir)
            assert tmp_dir != os.getcwd()

    def test_get_test_data_dir(self, delete_tmp_dir):
        with temp_dir(delete_tmp_dir) as tmp_dir:
            with pytest.raises(RuntimeError, match=r"test_data directory not found."):
                get_test_data_dir(tmp_dir)

        with temp_dir(delete_tmp_dir) as tmp_dir:
            module_dir = os.path.join(
                tmp_dir, "tests", "test_data", "tests", "test_data"
            )
            makedirs_p(module_dir)
            with pytest.raises(
                RuntimeError, match=r"Found multiple turbomoleio directories."
            ):
                get_test_data_dir(module_dir)

    def test_touch_file(self, delete_tmp_dir):
        with temp_dir(delete_tmp_dir):
            touch_file("test_file")
            assert os.path.isfile("test_file")

    def test_assert_MSONable(self):
        m = MSONableExample(1, 2)
        assert_MSONable(m)
        m = ExplicitAsFromDictExample(5, 12)
        assert not isinstance(m, MSONable)
        assert_MSONable(m, test_if_subclass=False)

    def test_gisnan(self):
        assert gisnan(np.nan)
        assert not gisnan(1)

    def test_get_tfp(self, test_data):
        """
        tesing the test file provider tool
        """
        fpath = get_tfp(os.path.join("structures", "methanol"), test_data=test_data)
        assert os.path.isfile(fpath)
        spath = get_sp("methanol", test_data=test_data)
        assert os.path.isfile(spath)
        assert fpath == spath
        dpath = get_tfp(test_data=test_data)
        assert os.path.isdir(dpath)
        cpath = get_control_integration("dscf", test_data=test_data)
        assert cpath == get_tfp(
            os.path.join("integration", "control", "dscf"), test_data=test_data
        )

    def test_assert_almost_equal(self):
        d1 = {"a": 1}
        d2 = {"b": 1}
        with pytest.raises(AssertionError):
            assert_almost_equal(d1, d2)

        d1 = {"a": 1}
        d2 = {"a": 1}
        assert_almost_equal(d1, d2)

        d1 = {"a": {"b": 1.000001}}
        d2 = {"a": {"b": 1}}
        with pytest.raises(AssertionError):
            assert_almost_equal(d1, d2, rtol=1e-8)
        assert_almost_equal(d1, d2, rtol=1e-4)

        d1 = {"a": {"b": np.array([1.000001, 2, 3])}}
        d2 = {"a": {"b": [1, 2, 3]}}
        with pytest.raises(AssertionError):
            assert_almost_equal(d1, d2, rtol=1e-8)
        assert_almost_equal(d1, d2, rtol=1e-4)

        d1 = {"a": "x"}
        d2 = {"a": "y"}
        with pytest.raises(AssertionError):
            assert_almost_equal(d1, d2)
        d2 = {"a": "x"}
        assert_almost_equal(d1, d2)

        d1 = {"a": "x"}
        d2 = {"a": 1}
        with pytest.raises(AssertionError):
            assert_almost_equal(d1, d2)

        d1 = {"a": 1.01}
        d2 = {"a": "x"}
        with pytest.raises(AssertionError):
            assert_almost_equal(d1, d2)

        d1 = {"a": 1.0001 + 1.00001j}
        d2 = {"a": 1 + 1j}
        with pytest.raises(AssertionError):
            assert_almost_equal(d1, d2, rtol=1e-8)
        assert_almost_equal(d1, d2, rtol=1e-3)

        d1 = 1
        d2 = {"b": 1}
        with pytest.raises(AssertionError):
            assert_almost_equal(d1, d2)

        d1 = 1
        d2 = complex(1, 0)
        assert_almost_equal(d1, d2)
        assert_almost_equal(d2, d1)

        d1 = complex(3, 5)
        d2 = complex(3, 5)
        assert_almost_equal(d1, d2)
        assert_almost_equal(d2, d1)

        d1 = [np.nan]
        d2 = np.nan
        with pytest.raises(AssertionError):
            assert_almost_equal(d1, d2)
        with pytest.raises(AssertionError):
            assert_almost_equal(d2, d1)

        d1 = 0.0
        d2 = -0.0
        assert_almost_equal(d1, d2)
        assert_almost_equal(d2, d1)

        d1 = [np.datetime64("NaT")]
        d2 = [np.datetime64("NaT")]
        assert_almost_equal(d1, d2)

        d1 = [np.timedelta64("NaT")]
        d2 = [np.timedelta64("NaT")]
        assert_almost_equal(d1, d2)

        d1 = [np.timedelta64("NaT")]
        d2 = [np.datetime64("NaT")]
        with pytest.raises(AssertionError):
            assert_almost_equal(d1, d2)
        with pytest.raises(AssertionError):
            assert_almost_equal(d2, d1)

    def test_compare_differences(self):
        diffs = compare_differences({}, {})
        assert len(diffs) == 0

        diffs = compare_differences([], [])
        assert len(diffs) == 0

        diffs = compare_differences({}, [])
        assert len(diffs) == 1
        assert diffs[0][0] == [
            ("root", "<class 'list'>"),
        ]
        assert diffs[0][1].startswith(f">>>{REF_SEQUENCE_TEST_OTHER}<<<")

        diffs = compare_differences([], {})
        assert len(diffs) == 1
        assert diffs[0][0] == [
            ("root", "<class 'dict'>"),
        ]
        assert diffs[0][1].startswith(f">>>{REF_DICT_TEST_OTHER}<<<")

        diffs = compare_differences(
            {1: [{}, {3: 5, 4: {"hello": [0, 1, 2]}}], 2: 3},
            {1: [{2: 4}, {3: 5, 4: {"hello": [0, 1, 5]}}]},
        )
        assert len(diffs) == 3
        k1 = [
            ("root", "<class 'dict'>"),
        ]
        k2 = [("root", "<class 'dict'>"), (1, "<class 'list'>"), (0, "<class 'dict'>")]
        k3 = [
            ("root", "<class 'dict'>"),
            (1, "<class 'list'>"),
            (1, "<class 'dict'>"),
            (4, "<class 'dict'>"),
            ("hello", "<class 'list'>"),
        ]
        diffs_levels = [diff[0] for diff in diffs]
        assert k1 in diffs_levels
        assert k2 in diffs_levels
        assert k3 in diffs_levels
        assert diffs[diffs_levels.index(k1)][1].startswith(
            f">>>{DICT_DIFFERENT_KEYS}<<<"
        )
        assert diffs[diffs_levels.index(k2)][1].startswith(
            f">>>{DICT_DIFFERENT_KEYS}<<<"
        )
        assert diffs[diffs_levels.index(k3)][1].startswith(f">>>{ARRAYS_DIFFER}<<<")

        diffs = compare_differences({"hi": "tata"}, {"hi": "titi"})
        assert len(diffs) == 1
        assert diffs[0][0] == [("root", "<class 'dict'>"), ("hi", "<class 'str'>")]
        assert diffs[0][1].startswith(f">>>{STRINGS_DIFFER}<<<")

        diffs = compare_differences({"hi": []}, {"hi": [1]})
        assert len(diffs) == 1
        assert diffs[0][0] == [("root", "<class 'dict'>"), ("hi", "<class 'list'>")]
        assert diffs[0][1].startswith(f">>>{ARRAYS_DIFFER}<<<")

        diffs = compare_differences(
            [1, {"a": 0.5, "b": 1.0}], [1, {"a": 0.5, "b": 1.01}]
        )
        assert len(diffs) == 1
        assert diffs[0][0] == [
            ("root", "<class 'list'>"),
            (1, "<class 'dict'>"),
            ("b", "<class 'float'>"),
        ]
        assert diffs[0][1].startswith(f">>>{NUMBERS_DIFFER}<<<")

        diffs = compare_differences(
            [1, {"a": 0.5, "b": 10.0}], [1, {"a": 0.5, "b": 10.1}], rtol=0.02
        )
        assert len(diffs) == 0

        diffs = compare_differences(
            [1, {"a": 0.5, "b": np.array([0.5, 0.9, 1.2])}],
            [1, {"a": 0.5, "b": np.array([0.5, 0.9, 1.21])}],
        )
        assert len(diffs) == 1
        assert diffs[0][0] == [
            ("root", "<class 'list'>"),
            (1, "<class 'dict'>"),
            ("b", "<class 'numpy.ndarray'>"),
        ]
        assert diffs[0][1].startswith(f">>>{ARRAYS_DIFFER}<<<")

        diffs = compare_differences(
            [1, {"a": 0.5, "b": np.array([0.5, 0.9, 1.2])}],
            [1, {"a": 0.5, "b": np.array([0.5, 0.9, 1.21])}],
            rtol=0.01,
        )
        assert len(diffs) == 0

        diffs = compare_differences(
            [
                1,
                {
                    "a": 0.5,
                    "b": Molecule(
                        species=["H", "H"], coords=[[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]]
                    ),
                },
            ],
            [
                1,
                {
                    "a": 0.5,
                    "b": Molecule(
                        species=["H", "H"], coords=[[-0.6, 0.0, 0.0], [0.6, 0.0, 0.0]]
                    ),
                },
            ],
        )
        assert len(diffs) == 1
        assert diffs[0][0] == [
            ("root", "<class 'list'>"),
            (1, "<class 'dict'>"),
            ("b", "<class 'pymatgen.core.structure.Molecule'>"),
        ]
        assert diffs[0][1].startswith(f">>>{OBJECTS_DIFFER}<<<")

        diffs = compare_differences(
            [
                1,
                {
                    "a": 0.5,
                    "b": Molecule(
                        species=["H", "H"], coords=[[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]]
                    ),
                },
            ],
            [
                1,
                {
                    "a": 0.5,
                    "b": Molecule(
                        species=["H", "H"], coords=[[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]]
                    ),
                },
            ],
        )
        assert len(diffs) == 0

        diffs = compare_differences(
            [
                1,
                Molecule(
                    species=["H", "H"], coords=[[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]]
                ),
            ],
            [1],
        )
        assert len(diffs) == 1
        assert diffs[0][0] == [("root", "<class 'list'>")]
        assert diffs[0][1].startswith(f">>>{SEQUENCE_DIFFERENT_SIZES}<<<")

        diffs = compare_differences({"hi": 3}, {"hi": "titi"})
        assert len(diffs) == 1
        assert diffs[0][0] == [("root", "<class 'dict'>"), ("hi", "<class 'str'>")]
        assert diffs[0][1].startswith(f">>>{REF_STRING_TEST_OTHER}<<<")

        diffs = compare_differences({"hi": "titi"}, {"hi": 3})
        assert len(diffs) == 1
        assert diffs[0][0] == [("root", "<class 'dict'>"), ("hi", "<class 'int'>")]
        assert diffs[0][1].startswith(f">>>{REF_NUMBER_TEST_OTHER}<<<")

        diffs = compare_differences(
            {"hi": "titi"},
            {
                "hi": Molecule(
                    species=["H", "H"], coords=[[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]]
                )
            },
        )
        assert len(diffs) == 1
        assert diffs[0][0] == [
            ("root", "<class 'dict'>"),
            ("hi", "<class 'pymatgen.core.structure.Molecule'>"),
        ]
        assert diffs[0][1].startswith(f">>>{TEST_STRING_REF_OTHER}<<<")

        diffs = compare_differences(
            {"hi": 3.5},
            {
                "hi": Molecule(
                    species=["H", "H"], coords=[[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]]
                )
            },
        )
        assert len(diffs) == 1
        assert diffs[0][0] == [
            ("root", "<class 'dict'>"),
            ("hi", "<class 'pymatgen.core.structure.Molecule'>"),
        ]
        assert diffs[0][1].startswith(f">>>{TEST_NUMBER_REF_OTHER}<<<")
