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

"""Configuration, fixtures and utilities for pytest testing infrastructure."""

import os

import pytest
from pymatgen.core.structure import Molecule, Structure

from turbomoleio.core.control import Control

DRYRUN_FPATH = os.path.join(os.path.split(__file__)[0], "dryrun_itest.json")


@pytest.fixture(scope="session")
def test_data():
    """Get path to the test_data directory."""
    from pathlib import Path

    module_dir = Path(__file__).resolve().parent
    test_data_dir = module_dir / "test_data"
    return test_data_dir.resolve()


@pytest.fixture
def control_filepath(control_filename, test_data):
    """Get path to control file from file name."""
    return test_data / "control" / control_filename


@pytest.fixture
def control(control_filename, test_data):
    """Get Control object from file name."""
    return Control.from_file(test_data / "control" / control_filename)


@pytest.fixture
def molecule_filepath(molecule_filename, test_data):
    """Get path to molecule file."""
    return test_data / "structures" / molecule_filename


@pytest.fixture
def molecule(molecule_filename, test_data):
    """Get Molecule object from file name."""
    return Molecule.from_file(test_data / "structures" / molecule_filename)


@pytest.fixture
def structure(structure_filename, test_data):
    """Get Structure object from file name."""
    return Structure.from_file(test_data / "structures" / structure_filename)


@pytest.fixture
def structure_filepath(structure_filename, test_data):
    """Get path to structure file."""
    return test_data / "structures" / structure_filename


@pytest.fixture
def mo_dirpath(dir_name, test_data):
    """Get molecular orbital directory path from directory name."""
    return test_data / "mo" / dir_name


@pytest.fixture(scope="session")
def delete_tmp_dir(request):
    """Get option whether to delete the temporary directory or not."""
    return not request.config.getoption("--keep-tmpdir", False)


def pytest_configure(config):
    """Configure pytest."""
    from turbomoleio.testing import ItestConfig

    ItestConfig.define_timeout = config.getoption("--define-timeout")
    ItestConfig.generate_ref = config.getoption("--generate-itest-ref")
    ItestConfig.dryrun = config.getoption("--dryrun-itest")
    dryrun_fpath = config.getoption("--dryrun-fpath")
    ItestConfig.dryrun_fpath = os.path.join(os.path.split(__file__)[0], dryrun_fpath)
    ItestConfig.dryrun_use_ref_control = config.getoption(
        "--dryrun-use-reference-control"
    )
    ItestConfig.tol = config.getoption("--itest-tol")
    ItestConfig.delete_tmp_dir = not config.getoption("--keep-tmpdir")
    # When running in dry mode, only the integration tests should be run and
    # there should not be any failures.
    if ItestConfig.dryrun:
        config.option.markexpr = "integration"
        config.option.maxfail = 1


def pytest_addoption(parser):
    """Add options to pytest."""
    parser.addoption(
        "--keep-tmpdir",
        action="store_true",
        default=False,
        help="Temporary directories used during the tests will not be deleted.",
    )

    parser.addoption(
        "--define-timeout",
        default=10,
        type=int,
        help="The number of seconds used for the timeout in DefineRunner "
        "for integration tests.",
    )

    parser.addoption(
        "--generate-itest-ref",
        action="store_true",
        default=False,
        help="The output control file generated during the integration "
        "tests will be copied to the reference folder. "
        "N.B. this will overwrite previously existing files.",
    )

    parser.addoption(
        "--dryrun-itest",
        action="store_true",
        default=False,
        help="The files generated during the integration tests will be compared "
        "against reference files. A list of all differences will be "
        "provided in a file.",
    )

    parser.addoption(
        "--dryrun-use-reference-control",
        action="store_true",
        default=False,
        help="Use the reference control file for integration tests "
        "instead of the one generated using define. "
        "The control file is still generated and compared.",
    )

    parser.addoption(
        "--dryrun-fpath",
        default=DRYRUN_FPATH,
        type=str,
        help="Filepath for the results of the differences of the dryrun.",
    )

    parser.addoption(
        "--itest-tol",
        default=1e-4,
        type=float,
        help="The absolute tolerance used in the integration test to match "
        "floating point numbers when comparing to reference files.",
    )
