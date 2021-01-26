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


files_list = [("dscf", "h2o_std"), ("dscf", "h2o_uhf"), ("dscf", "nh3_cosmo_fermi"),
              ("dscf", "nh3_dftd1"), ("dscf", "aceton_dftd3_tzvp"),("ridft", "h2o_dftd2_marij"),
              ("ridft", "h2o_dftd3-bj_not_conv"), ("ridft", "nh3_rijk_xcfun_m06"),
              ("ridft", "b28_many_irreps"), ("grad", "h2o_std"), ("rdgrad", "h2o_dftd3-bj"),
              ("relax", "h2o_internal"), ("relax", "h2o_cartesian"), ("relax", "no_version_header"),
              ("statpt", "h3cbr_internal"), ("statpt", "aceton_cartesian"),
              ("escf", "Al6_columns"), ("escf", "h2o_ridft_cosmo"), ("escf", "h2o_ridft_rpat"),
              ("egrad", "h2o_sym"), ("egrad", "h3cbr_nosym"), ("aoforce", "aceton_full"),
              ("aoforce", "h2_numforce")]


parser_methods = ["all_done", "header", "centers", "coordinates", "basis", "symmetry",
                  "cosmo_header", "density_functional_data", "rij_info", "dftd", "pre_scf_run",
                  "scf_iterations", "scf_energies", "cosmo_results", "electrostatic_moments",
                  "timings", "s2", "is_uhf", "fermi", "integral", "pre_escf_run", "escf_iterations",
                  "escf_gs_total_en", "escf_excitations", "rdgrad_memory", "gradient", "egrad_excited_state",
                  "statpt_info", "relax_info", "relax_gradient_values", "relax_conv_info",
                  "aoforce_numerical_integration", "aoforce_analysis"]


@pytest.fixture(scope="function", params=files_list, ids=[os.path.join(*f) for f in files_list])
def parser_and_dict(request, testdir):
    directory = request.param[0]
    name = request.param[1]
    path = os.path.join(testdir, "outputs", directory, name)
    parser = Parser.from_file(path+".log")
    with open(path+".json") as f:
        d = json.load(f)

    return parser, d


@pytest.fixture(scope="function", params=parser_methods)
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

    def test_get_split_jobex_parsers(self, testdir):
        path = os.path.join(testdir, "outputs", "jobex", "h2o_dscf_job.last")
        p = Parser.from_file(path)
        jp = p.get_split_jobex_parsers()
        assert jp.exec_en == "dscf"
        assert jp.exec_grad == "grad"
        assert jp.exec_relax == "statpt"

        path = os.path.join(testdir, "outputs", "jobex", "no3_ridft_job.last")
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
            reference data will be generated. If None "parser_methods" will be used.
        overwrite (bool): if False, in case a method has already its reference value
            in the json file it will not be changed. If True it will be overwritten.
    """
    if methods is None:
        methods = parser_methods
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
