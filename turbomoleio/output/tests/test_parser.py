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
Test module for the Parser object.
The testing of all the properties is done through pytest fixtures with parameters.
The output obtained from the parsing are compared with those stored in json files.
A helper function is provided to generate the reference json files. These files
should then be verified by the user to ensure the correctness of the parser.
"""
import pytest
import os
import json

from turbomoleio.output.parser import Parser, convert_float, convert_int, convert_time_string
from turbomoleio.testfiles.utils import assert_almost_equal, temp_dir
from turbomoleio.testfiles.utils import PARSER_METHODS
from turbomoleio.testfiles.utils import TM_VERSIONS
from turbomoleio.testfiles.utils import TESTS_CONFIGS_TM_VERSIONS

excluded_execs = ['jobex']
files_list = [
    (tm_version, tm_exec, test_name)
    for tm_version in TM_VERSIONS
    for tm_exec, exec_tests in TESTS_CONFIGS_TM_VERSIONS[tm_version]['testlist'].items()
    if tm_exec not in excluded_execs
    for test_name in exec_tests
]


@pytest.fixture(scope="function", params=files_list, ids=[os.path.join(*f) for f in files_list])
def parser_and_dict(request, testdir):
    tm_version = request.param[0]
    tm_exec = request.param[1]
    test_name = request.param[2]
    path = os.path.join(testdir, "outputs", tm_version, tm_exec, test_name)
    parser = Parser.from_file(os.path.join(path, f"{tm_exec}.log"))
    with open(os.path.join(path, "ref_parser.json")) as f:
        d = json.load(f)

    return parser, d


@pytest.fixture(scope="function", params=PARSER_METHODS)
def method(request):
    return request.param


class TestParser:

    def test_properties(self, parser_and_dict, method):
        parser, desired = parser_and_dict

        if method not in desired:
            pytest.skip("Method {} is not present in the dictionary".format(method))

        parsed_data = getattr(parser, method)

        # ignore date values since in the dictionary they are datetime, while
        # just strings in the json file.
        assert_almost_equal(parsed_data, desired[method], rtol=1e-4,
                            ignored_values=["start_time", "end_time", "@version"])

    @pytest.mark.parametrize("tm_version", TM_VERSIONS)
    def test_get_split_jobex_parsers(self, testdir, tm_version):
        path = os.path.join(testdir, "outputs", tm_version, "jobex", "h2o_dscf", "job.last")
        p = Parser.from_file(path)
        jp = p.get_split_jobex_parsers()
        assert jp.exec_en == "dscf"
        assert jp.exec_grad == "grad"
        assert jp.exec_relax == "statpt"

        path = os.path.join(testdir, "outputs", tm_version, "jobex", "no3_ridft", "job.last")
        p = Parser.from_file(path)
        jp = p.get_split_jobex_parsers()
        assert jp.exec_en == "ridft"
        assert jp.exec_grad == "rdgrad"
        assert jp.exec_relax == "statpt"

    def test_grep_line(self):
        string = """some text
some more text
following line"""
        p = Parser(string)
        assert p.grep_line("some more", nlines=1) == "some more text\nfollowing line"
        assert p.grep_line("wrong text") is None

    def test_get_value(self):
        p = Parser("line to match 0.1 xxx")

        assert p.get_value("line to", -2, 0, float) == pytest.approx(0.1)

    def test_fail_all_done_check(self, delete_tmp_dir):
        with temp_dir(delete_tmp_dir) as tmp:
            with open("test.log", "wt") as f:
                f.write("""some text
some more text
following line""")

            with pytest.raises(ValueError, match="The string does not contain data for a completed calculation"):
                Parser.from_file("test.log")


class TestFunctions:

    def test_convert_float(self):
        assert convert_float("0.1") == pytest.approx(0.1)
        assert convert_float("0.1E-01") == pytest.approx(0.01)
        assert convert_float("0.1D-02") == pytest.approx(0.001)
        assert convert_float("****") is None

        with pytest.raises(ValueError, match="could not convert string to float*"):
            convert_float("xxxx")

    def test_convert_int(self):
        assert convert_int("1") == 1
        assert convert_int("10000") == 10000
        assert convert_int("****") is None

        with pytest.raises(ValueError, match="invalid literal for int*"):
            convert_int("xxxx")

    def test_convert_time_string(self):
        assert convert_time_string("1 hours 37 minutes and 28 seconds") == 5848
        assert convert_time_string("2 days 1 hours 0 minutes and 0 seconds") == 176400


def generate_files(files=None, methods=None, overwrite=False):
    """
    Helper function to generate the reference files for the test of the parser.
    Allows to target only specific methods and specific files.

    Args:
        files (list): list of tuples with (folder name, file name without extension),
            like the one in "files_list". Only the json for the files in this list
            will be generated. If None "files_list" will be used.
        methods (list): list of string with the methods of Parser for which the
            reference data will be generated. If None "PARSER_METHODS" will be used.
        overwrite (bool): if False, in case a method has already its reference value
            in the json file it will not be changed. If True it will be overwritten.
    """
    if methods is None:
        methods = PARSER_METHODS
    if files is None:
        files = files_list

    for directory, name in files:
        path = os.path.join(os.path.split(__file__)[0], "../../testfiles", "outputs", directory, name)
        parser = Parser.from_file(path+".log")
        json_path = path+".json"
        if os.path.isfile(json_path):
            with open(json_path) as f:
                ref_data = json.load(f)
        else:
            ref_data = {}

        for m in methods:
            if m not in ref_data or overwrite:
                parsed_data = getattr(parser, m)
                ref_data[m] = parsed_data

        with open(json_path, "wt") as f:
            from monty.json import jsanitize
            json.dump(jsanitize(ref_data), f, indent=2)
