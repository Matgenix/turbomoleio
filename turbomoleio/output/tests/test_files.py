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
Module to verify the parsing of the files corresponding to the different types
of calculations. Implicitly tests the construction of the Data object. The tests
will not be repeated for the Data objects separately
"""

import pytest
import os
import json
from monty.serialization import loadfn

from turbomoleio.output.files import exec_to_out_obj, JobexOutput, EscfOnlyOutput, EscfOutput
from turbomoleio.testfiles.utils import assert_almost_equal
from turbomoleio.testfiles.utils import TM_VERSIONS
from turbomoleio.testfiles.utils import TESTS_CONFIGS_TM_VERSIONS


files_list = [
    (tm_version, tm_exec, test_name)
    for tm_version in TM_VERSIONS
    for tm_exec, exec_tests in TESTS_CONFIGS_TM_VERSIONS[tm_version]['testlist'].items()
    for test_name in exec_tests
]


@pytest.fixture(scope="function", params=files_list, ids=[os.path.join(*f) for f in files_list])
def cls_dict_path(request, testdir):
    tm_version = request.param[0]
    tm_exec = request.param[1]
    test_name = request.param[2]
    path = os.path.join(testdir, "outputs", tm_version, tm_exec, test_name)
    output_clses = []
    ref_dicts = []
    fpaths = []
    if tm_exec == "escf":
        output_clses.append(EscfOnlyOutput)
        ref_dicts.append(loadfn(os.path.join(path, "ref_escf_output.json"), cls=None))
        fpaths.append(os.path.join(path, f"{tm_exec}.log"))
        output_clses.append(EscfOutput)
        ref_dicts.append(loadfn(os.path.join(path, "ref_output.json"), cls=None))
        fpaths.append(os.path.join(path, f"{tm_exec}.log"))
    elif tm_exec == "jobex":
        output_clses.append(JobexOutput)
        ref_dicts.append(loadfn(os.path.join(path, "ref_output.json"), cls=None))
        fpaths.append(os.path.join(path, "job.last"))
    else:
        output_clses.append(exec_to_out_obj[tm_exec])
        ref_dicts.append(loadfn(os.path.join(path, "ref_output.json"), cls=None))
        fpaths.append(os.path.join(path, f"{tm_exec}.log"))

    return output_clses, ref_dicts, fpaths


class TestFiles:

    def test_properties(self, cls_dict_path):
        output_clses, desired_list, paths_list = cls_dict_path

        for output_cls, desired, path in zip(output_clses, desired_list, paths_list):
            parsed_data = output_cls.from_file(path).as_dict()

            # ignore date values since in the dictionary they are datetime, while
            # just strings in the json file.
            assert_almost_equal(parsed_data, desired, rtol=1e-4,
                                ignored_values=["start_time", "end_time", "@version"])


def generate_files(files=None, overwrite=False):
    """
    Helper function to generate the reference files for the test of the files objects.
    Allows to target only specific files.

    Args:
        files (list): list of tuples with (folder name, file name without extension),
            like the one in "files_list". Only the json for the files in this list
            will be generated. If None "files_list" will be used.
        overwrite (bool): if False, in case a method has already its reference value
            in the json file it will not be changed. If True it will be overwritten.
    """
    if files is None:
        files = files_list

    for ref_name, name in files:
        if ref_name == "escf_only":
            directory = "escf"
        else:
            directory = ref_name
        path = os.path.join(os.path.split(__file__)[0], "../../testfiles", "outputs", directory, name)
        if ref_name == "escf_only":
            json_path = path + "_outfile_escf_only.json"
        else:
            json_path = path + "_outfile.json"

        if os.path.isfile(json_path) and not overwrite:
            continue

        exec_dict = dict(exec_to_out_obj)
        exec_dict["jobex"] = JobexOutput
        exec_dict["escf_only"] = EscfOnlyOutput

        output_cls = exec_dict[ref_name]
        file_path = path + ".log" if ref_name != "jobex" else path + "_job.last"

        parsed_data = output_cls.from_file(file_path).as_dict()

        with open(json_path, "wt") as f:
            from monty.json import jsanitize
            json.dump(jsanitize(parsed_data), f, indent=2)
