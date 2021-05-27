import os

import pytest

from pymatgen.core.structure import Molecule

from turbomoleio.core.control import Control

TESTDIR = os.path.join(os.path.split(__file__)[0],
                       'turbomoleio',
                       'testfiles')


@pytest.fixture(autouse=True, scope='session')
def testdir():
    return TESTDIR


@pytest.fixture
def control_filepath(control_filename):
    return os.path.join(TESTDIR, 'control', control_filename)


@pytest.fixture
def control(control_filename):
    return Control.from_file(os.path.join(TESTDIR,
                                          'control',
                                          control_filename))


@pytest.fixture
def molecule_filepath(molecule_filename):
    return os.path.join(TESTDIR, 'structures', molecule_filename)


@pytest.fixture
def molecule(molecule_filename):
    return Molecule.from_file(os.path.join(TESTDIR,
                                           'structures',
                                           molecule_filename))


@pytest.fixture
def mo_dirpath(dir_name):
    return os.path.join(TESTDIR, 'mo', dir_name)


@pytest.fixture(scope="session")
def delete_tmp_dir(request):
    return not request.config.getoption("--keep-tmpdir", False)


def pytest_configure(config):
    from turbomoleio.testfiles.utils import ItestConfig
    ItestConfig.define_timeout = config.getoption("--define-timeout")
    ItestConfig.generate_ref = config.getoption("--generate-itest-ref")
    ItestConfig.dryrun = config.getoption("--dryrun-itest")
    ItestConfig.tol = config.getoption("--itest-tol")
    ItestConfig.delete_tmp_dir = not config.getoption("--keep-tmpdir")


def pytest_addoption(parser):

    parser.addoption("--keep-tmpdir", action="store_true", default=False,
                     help="Temporary directories used during the tests will not be deleted.")

    parser.addoption("--define-timeout", default=10, type=int,
                     help="The number of seconds used for the timeout in DefineRunner for integration tests.")

    parser.addoption("--generate-itest-ref", action="store_true", default=False,
                     help="The output control file generated during the integration tests will be copied to the "
                          "reference folder. N.B. this will overwrite previously existing files.")

    parser.addoption("--dryrun-itest", action="store_true", default=False,
                     help="The files generated during the integration tests will be compared against "
                          "reference files. A list of all differences will be provided.")

    parser.addoption("--itest-tol", default=1e-4, type=float,
                     help="The absolute tolerance used in the integration test to match floating point"
                          "numbers when comparing to reference files.")
