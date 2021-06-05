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

from turbomoleio.output.files import exec_to_out_obj, JobexOutput, EscfOnlyOutput
from turbomoleio.testfiles.utils import assert_almost_equal
from turbomoleio.testfiles.utils import TM_VERSIONS


files_list = [("dscf", "h2o_std"), ("dscf", "h2o_uhf"), ("dscf", "nh3_cosmo_fermi"),
              ("dscf", "nh3_dftd1"), ("dscf", "aceton_dftd3_tzvp"),("ridft", "h2o_dftd2_marij"),
              ("ridft", "h2o_dftd3-bj_not_conv"), ("ridft", "nh3_rijk_xcfun_m06"),
              ("ridft", "b28_many_irreps"), ("grad", "h2o_std"), ("rdgrad", "h2o_dftd3-bj"),
              ("relax", "h2o_internal"), ("relax", "h2o_cartesian"), ("relax", "no_version_header"),
              ("statpt", "h3cbr_internal"), ("statpt", "aceton_cartesian"),
              ("escf", "Al6_columns"), ("escf", "h2o_ridft_cosmo"), ("escf", "h2o_ridft_rpat"),
              ("egrad", "h2o_sym"), ("egrad", "h3cbr_nosym"), ("aoforce", "aceton_full"),
              ("aoforce", "h2_numforce"), ("jobex", "h2o_dscf"), ("jobex", "no3_ridft"),
              ("escf_only", "h2o_ridft_cosmo")]


@pytest.fixture(scope="function", params=files_list, ids=[os.path.join(*f) for f in files_list])
def cls_dict_path(request, testdir, tm_version):
    ref_name = request.param[0]
    if ref_name == "escf_only":
        directory = "escf"
    else:
        directory = ref_name
    name = request.param[1]
    path = os.path.join(testdir, "outputs", tm_version, directory, name)

    if ref_name == "escf_only":
        json_path = os.path.join(path, "ref_escf_output.json")
    else:
        json_path = os.path.join(path, "ref_output.json")

    with open(json_path) as f:
        d = json.load(f)

    exec_dict = dict(exec_to_out_obj)
    exec_dict["jobex"] = JobexOutput
    exec_dict["escf_only"] = EscfOnlyOutput

    # file_path = path + ".log" if ref_name != "jobex" else path + "_job.last"
    file_path = os.path.join(
        path, f"{directory}.log"
    ) if ref_name != "jobex" else os.path.join(
        path, "job.last"
    )

    return exec_dict[ref_name], d, file_path


class TestParser:

    @pytest.mark.parametrize("tm_version", TM_VERSIONS)
    def test_properties(self, cls_dict_path):
        output_cls, desired, path = cls_dict_path

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
