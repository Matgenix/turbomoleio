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

from __future__ import division, print_function, unicode_literals
import os
import tempfile
import subprocess
import shlex
import shutil
import json
import inspect
import numbers
import numpy as np
from contextlib import contextmanager
from monty.os import cd, makedirs_p
from monty.json import MSONable, MontyDecoder
from monty.serialization import dumpfn, loadfn
from turbomoleio.input.define import DefineRunner
from turbomoleio.input.define import DefineError
from turbomoleio.input.utils import get_define_template
from turbomoleio.output.states import EigerRunner, States
from turbomoleio.output.files import JobexOutput, EscfOnlyOutput
from turbomoleio.core.control import Control
from turbomoleio.core.control import adg, cdg
from turbomoleio.output.files import exec_to_out_obj
from turbomoleio.output.parser import Parser


TESTDIR = os.path.split(__file__)[0]
TM_VERSIONS = ['TM_v7.3', 'TM_v7.3.1']
TESTS_CONFIGS_TM_VERSIONS = {
    tmv: loadfn(os.path.join(TESTDIR, 'outputs', tmv, 'tests_config.yaml'))
    for tmv in TM_VERSIONS
}
PARSER_METHODS = ["all_done", "header", "centers", "coordinates", "basis", "symmetry",
                  "cosmo_header", "density_functional_data", "rij_info", "dftd", "pre_scf_run",
                  "scf_iterations", "scf_energies", "cosmo_results", "electrostatic_moments",
                  "timings", "s2", "is_uhf", "fermi", "integral", "pre_escf_run", "escf_iterations",
                  "escf_gs_total_en", "escf_excitations", "rdgrad_memory", "gradient", "egrad_excited_state",
                  "statpt_info", "relax_info", "relax_gradient_values", "relax_conv_info",
                  "aoforce_numerical_integration", "aoforce_analysis", "mp2_results", "riper_scf_energies",
                  "periodicity_data"]


class ItestError(BaseException):
    """
    error to be raised if an itest fails
    """


@contextmanager
def temp_dir(delete, changedir=True):
    """
    Context manager that creates a temporary directory with tempfile.mkdtemp and cd to it with monty.os.cd.

    Args:
        delete (bool): if True the directory will be deleted at the end of the job, if False it will be preserved.
        changedir (bool): if True inside the context manager will make a cd to the temporary directoy.

    Yields:
        the path to the temporary directory created.
    """

    testdir = tempfile.mkdtemp()
    if not delete:
        print("Running folder: {}".format(testdir))

    try:
        if changedir:
            with cd(testdir):
                yield testdir
        else:
            yield testdir
    finally:
        if delete:
            shutil.rmtree(testdir, ignore_errors=True)


def get_control_filepath(filename):
    """
    The path to a reference control file in the testfiles/control folder.

    Args:
        filename (str): the name of the reference control file.

    Returns:
        str: the absolute path to the coord file.
    """
    return os.path.join(TESTDIR, 'control', filename)


def touch_file(path):
    """
    Touches a file at the specified path.

    Args:
        path (str): path to the file.

    Returns:
        None
    """
    with open(path, 'a'):
        os.utime(path, None)


def assert_MSONable(obj, test_if_subclass=True):
    """
    Tests if obj is MSONable and tries to verify whether the contract is
    fulfilled. Tries to convert an object to a dictionary and back and
    checking if the dictionaries are equivalent

    By default, the method tests whether obj is an instance of MSONable.
    This check can be deactivated by setting test_if_subclass to False.

    Args:
        obj: the object to be tests
        test_if_subclass:

    """
    if test_if_subclass:
        assert isinstance(obj, MSONable)

    # applies as_dict/from dict and compares the resulting dictionaries
    # use numpy assert_equal to also match numpy arrays
    np.testing.assert_equal(obj.as_dict(), obj.__class__.from_dict(obj.as_dict()).as_dict())
    # tests that the string is compatible with json
    assert json.loads(obj.to_json(), cls=MontyDecoder)


def gisnan(x):
    """
    Auxiliary function imported from numpy.testing._private version 16.2.
    Used in assert_almost_equal. Imported here since it is not exposed
    in the public interface.

    like isnan, but always raise an error if type not supported instead of
    returning a TypeError object.
    Notes
    -----
    isnan and other ufunc sometimes return a NotImplementedType object instead
    of raising any exception. This function is a wrapper to make sure an
    exception is always raised.
    This should be removed once this problem is solved at the Ufunc level.
    """
    st = np.core.isnan(x)
    if isinstance(st, type(NotImplemented)):
        raise TypeError("isnan not supported for this type")
    return st


def assert_almost_equal(actual, desired, rtol=1e-7, atol=0, ignored_values=None, err_msg='', verbose=True):
    """
    Function imported from numpy.testing: assert_equal. Version 16.2.
    Two key modifications compared to the original implementation:
    1) allow comparison of numbers with a tolerance on the difference (other functions
        in numpy that allow a tolerance as an argument do not support comparison
        between dictionaries).
    2) allow to skip the explicit comparison of some attributes in dictionaries.

    Raises an AssertionError if two objects are not equal within the required tolerances.
    Given two objects (scalars, lists, tuples, dictionaries or numpy arrays),
    check that all elements of these objects are almost equal. An exception is raised
    at the first conflicting values.
    Comparison for numerical values is performed with assert_allclose

    Args:
        actual: the object to check.
        desired : the expected object.
        rtol (float): relative tolerance.
        atol (float): absolute tolerance.
        ignored_values (list): if a comparison between two dictionaries, keywords contained
            in this list will not be compared (the key should still exist in both the
            dictionaries though).
        err_msg (str): the error message to be printed in case of failure.
        verbose (bool): if True, the conflicting values are appended to the error message.

    Raises:
        AssertionError: if actual and desired are not equal.
    """
    __tracebackhide__ = True  # Hide traceback for py.test
    if isinstance(desired, dict):
        if not isinstance(actual, dict):
            raise AssertionError(repr(type(actual)))
        assert_almost_equal(len(actual), len(desired), rtol, atol, ignored_values, err_msg, verbose)

        if ignored_values is None:
            ignored_values = []

        for k, i in desired.items():
            if k not in actual:
                raise AssertionError(repr(k))

            # don't check nested values if the key belong to the list to skip
            if k in ignored_values:
                continue

            assert_almost_equal(actual[k], desired[k], rtol, atol, ignored_values, 'key=%r\n%s' % (k, err_msg), verbose)
        return
    if isinstance(desired, (list, tuple)) and isinstance(actual, (list, tuple)):
        assert_almost_equal(len(actual), len(desired), rtol, atol, ignored_values, err_msg, verbose)
        for k in range(len(desired)):
            assert_almost_equal(actual[k], desired[k], rtol, atol, ignored_values, 'item=%r\n%s' % (k, err_msg), verbose)
        return
    from numpy.core import ndarray, isscalar, signbit
    from numpy.lib import iscomplexobj, real, imag
    if isinstance(actual, ndarray) or isinstance(desired, ndarray):
        return np.testing.assert_allclose(actual, desired, rtol=rtol, atol=atol, err_msg=err_msg, verbose=verbose)
    msg = np.testing.build_err_msg([actual, desired], err_msg, verbose=verbose)

    # Handle complex numbers: separate into real/imag to handle
    # nan/inf/negative zero correctly
    # XXX: catch ValueError for subclasses of ndarray where iscomplex fail
    try:
        usecomplex = iscomplexobj(actual) or iscomplexobj(desired)
    except (ValueError, TypeError):
        usecomplex = False

    if usecomplex:
        if iscomplexobj(actual):
            actualr = real(actual)
            actuali = imag(actual)
        else:
            actualr = actual
            actuali = 0
        if iscomplexobj(desired):
            desiredr = real(desired)
            desiredi = imag(desired)
        else:
            desiredr = desired
            desiredi = 0
        try:
            np.testing.assert_allclose(actualr, desiredr, rtol=rtol, atol=atol)
            np.testing.assert_allclose(actuali, desiredi, rtol=rtol, atol=atol)
        except AssertionError:
            raise AssertionError(msg)

    # isscalar test to check cases such as [np.nan] != np.nan
    if isscalar(desired) != isscalar(actual):
        raise AssertionError(msg)

    a_is_number = isinstance(actual, numbers.Number)
    d_is_number = isinstance(desired, numbers.Number)

    if a_is_number != d_is_number:
        raise AssertionError(msg)

    if a_is_number:
        if np.isclose(actual, desired, rtol=rtol, atol=atol):
            return
        else:
            raise AssertionError(msg)

    # Inf/nan/negative zero handling
    try:
        isdesnan = gisnan(desired)
        isactnan = gisnan(actual)
        if isdesnan and isactnan:
            return  # both nan, so equal

        # allow 0.0 and -0.0 to match
        if desired == 0 and actual == 0:
            return
            # if not signbit(desired) == signbit(actual):
            #     raise AssertionError(msg)

    except (TypeError, ValueError, NotImplementedError):
        pass

    try:
        isdesnat = np.core.isnat(desired)
        isactnat = np.core.isnat(actual)
        dtypes_match = np.array(desired).dtype.type == np.array(actual).dtype.type
        if isdesnat and isactnat:
            # If both are NaT (and have the same dtype -- datetime or
            # timedelta) they are considered equal.
            if dtypes_match:
                return
            else:
                raise AssertionError(msg)

    except (TypeError, ValueError, NotImplementedError):
        pass

    try:
        # First check if they are equal with ==. If not Use assert allclose instead of __eq__
        if not (desired == actual):
            raise AssertionError(msg)
            # try:
            #     np.testing.assert_allclose(actual, desired, rtol=rtol, atol=atol)
            # except (AssertionError, TypeError):
            #     raise AssertionError(msg)

    except (DeprecationWarning, FutureWarning) as e:
        # this handles the case when the two types are not even comparable
        if 'elementwise == comparison' in e.args[0]:
            raise AssertionError(msg)
        else:
            raise


REF_DICT_TEST_OTHER = 'REF_DICT_TEST_OTHER'
REF_SEQUENCE_TEST_OTHER = 'REF_SEQUENCE_TEST_OTHER'
REF_STRING_TEST_OTHER = 'REF_STRING_TEST_OTHER'
REF_NUMBER_TEST_OTHER = 'REF_NUMBER_TEST_OTHER'
TEST_NUMBER_REF_OTHER = 'TEST_NUMBER_REF_OTHER'
TEST_STRING_REF_OTHER = 'TEST_STRING_REF_OTHER'
DICT_DIFFERENT_KEYS = 'DICT_DIFFERENT_KEYS'
SEQUENCE_DIFFERENT_SIZES = 'SEQUENCE_DIFFERENT_SIZES'
NUMBERS_DIFFER = 'NUMBERS_DIFFER'
ARRAYS_DIFFER = 'ARRAYS_DIFFER'
STRINGS_DIFFER = 'STRINGS_DIFFER'
OBJECTS_DIFFER = 'OBJECTS_DIFFER'
DIFFERENCE_TYPES = [REF_DICT_TEST_OTHER, REF_SEQUENCE_TEST_OTHER,
                    REF_STRING_TEST_OTHER, REF_NUMBER_TEST_OTHER,
                    TEST_NUMBER_REF_OTHER, TEST_STRING_REF_OTHER, DICT_DIFFERENT_KEYS,
                    SEQUENCE_DIFFERENT_SIZES,
                    NUMBERS_DIFFER, ARRAYS_DIFFER, STRINGS_DIFFER, OBJECTS_DIFFER]


def _update_differences(differences, level, difference_type, message=None):
    if difference_type not in DIFFERENCE_TYPES:
        raise ValueError(f'Difference type "{difference_type}" is not valid')
    msg = f'>>>{difference_type}<<<'
    if message is not None:
        msg += f'\n{message}'
    all_levels = [tuple(item[0]) for item in differences]
    level_t = tuple(level)
    if level_t in all_levels:
        raise RuntimeError('Should not be reached here')
    differences.append((level, msg))


# TODO: We should think about how to test eigenvectors, reduced_masses etc ...
ignored_dryrun_itest_parsed_keys = ["start_time", "end_time", "wall_time", "cpu_time", "construction_timings",
                                    "host", "version", "build", "eigenvectors", "reduced_masses", "electric_dipole",
                                    "magnetic_dipole", "electric_quadrupole",
                                    "@version"]


def compare_differences(actual, desired, rtol=1e-7, atol=0, current_level=None):
    """
    Function to list all differences between test object and reference objects to within some tolerance.

    Args:
        actual: the object to check.
        desired : the expected object.
        rtol (float): relative tolerance.
        atol (float): absolute tolerance.
        current_level: current_level in a nested comparison.

    Raises:
        AssertionError: if actual and desired are not equal.
    """
    __tracebackhide__ = True  # Hide traceback for py.test
    differences = []
    if current_level is None:
        current_level = [('root', str(type(desired)))]
    else:
        current_level = list(current_level)

    if isinstance(desired, dict):
        if not isinstance(actual, dict):
            _update_differences(differences, current_level, REF_DICT_TEST_OTHER,
                                message=f'Reference object is a {dict}, '
                                        f'tested object is a {type(actual)}')
            return differences
        if set(actual) != set(desired):
            _update_differences(differences, current_level, DICT_DIFFERENT_KEYS,
                                message=f'Keys in reference dict are different from keys in tested dict:\n'
                                        f' - Reference keys: {sorted(desired.keys())}\n'
                                        f' - Test keys: {sorted(actual.keys())}\n'
                                        f' - Keys in reference but not in test: '
                                        f'{sorted(set(desired).difference(set(actual)))}\n'
                                        f' - Keys in test but not in reference: '
                                        f'{sorted(set(actual).difference(set(desired)))}')

        common_keys = sorted(set(desired).intersection(set(actual)))
        for k in common_keys:
            if k in ignored_dryrun_itest_parsed_keys:
                continue
            lvl = list(current_level)
            lvl.append((k, str(type(desired[k]))))
            differences.extend(compare_differences(actual=actual[k], desired=desired[k],
                                                   rtol=rtol, atol=atol,
                                                   current_level=lvl))

        return differences

    if isinstance(desired, str):
        if not isinstance(actual, str):
            _update_differences(differences, current_level, REF_STRING_TEST_OTHER,
                                message=f'Reference object is a {type(desired)}, '
                                        f'tested object is a {type(actual)}')
            return differences

        if desired != actual:
            _update_differences(differences, current_level, STRINGS_DIFFER,
                                message=f' - Reference string: {desired}\n'
                                        f' - Test string: {actual}\n')
        return differences

    if isinstance(desired, numbers.Number):
        if not isinstance(actual, numbers.Number):
            _update_differences(differences, current_level, REF_NUMBER_TEST_OTHER,
                                message=f'Reference object is a number ({type(desired)}), '
                                        f'tested object is not a number ({type(actual)})')
            return differences

        if not np.isclose(actual, desired, rtol=rtol, atol=atol):
            _update_differences(differences, current_level, NUMBERS_DIFFER,
                                message=f' - Reference number: {desired}\n'
                                        f' - Test number: {actual}\n')
        return differences

    actual_np, desired_np = np.asanyarray(actual), np.asanyarray(desired)
    if issubclass(actual_np.dtype.type, numbers.Number) and issubclass(desired_np.dtype.type, numbers.Number):
        if actual_np.shape != desired_np.shape:
            _update_differences(differences, current_level, ARRAYS_DIFFER,
                                message=f'Reference and test arrays do not have the same shape\n'
                                        f' - Shape of reference array: {desired_np.shape}\n'
                                        f' - Shape of test array: {actual_np.shape}')
            return differences
        if not np.allclose(actual_np, desired_np, rtol=rtol, atol=atol):
            _update_differences(differences, current_level, ARRAYS_DIFFER,
                                message=f'Reference and test arrays are not equal to tolerance '
                                        f'rtol={rtol}, atol={atol}\n'
                                        f' - Reference array: {desired}\n'
                                        f' - Test array: {actual}')
        return differences

    if isinstance(desired, (list, tuple)):
        if not isinstance(actual, (list, tuple)):
            _update_differences(differences, current_level, REF_SEQUENCE_TEST_OTHER,
                                message=f'Reference object is a {type(desired)}, '
                                        f'tested object is a {type(actual)}')
            return differences
        if len(actual) != len(desired):
            _update_differences(differences, current_level, SEQUENCE_DIFFERENT_SIZES,
                                message=f'Number of items in reference list or tuple is {len(desired)},'
                                        f'number of items in tested list or tuple is {len(actual)}')
            #TODO: decide here if we test whether the first N are the same ? or if one is a subset of the other ?
            # In any case, we should give more information of what is in there
            return differences

        for i, ref in enumerate(desired):
            lvl = list(current_level)
            lvl.append((i, str(type(ref))))
            differences.extend(compare_differences(actual=actual[i], desired=ref,
                                                   rtol=rtol, atol=atol,
                                                   current_level=lvl))
            # TODO: decide here if we test whether there is a shuffling of the items ?
        return differences

    if isinstance(actual, numbers.Number) and not isinstance(desired, numbers.Number):
        _update_differences(differences, current_level, TEST_NUMBER_REF_OTHER,
                            message=f'Reference object is a {type(desired)}, '
                                    f'tested object is a number ({type(actual)})')
        return differences

    if isinstance(actual, str) and not isinstance(desired, str):
        _update_differences(differences, current_level, TEST_STRING_REF_OTHER,
                            message=f'Reference object is a {type(desired)}, '
                                    f'tested object is a {type(actual)}')
        return differences

    # If the reference and tested objects are not a dict, list, tuple, str, number or array of numbers,
    # compare the objects directly
    if actual != desired:
        _update_differences(differences, current_level, OBJECTS_DIFFER,
                            message=f'Reference and test objects are not equal.')
        return differences

    return differences


def has_matplotlib():
    """
    True if matplotlib is installed.
    """
    try:
        import matplotlib
        return True
    except ImportError:
        return False


def get_tfp(file_name=None):
    """
    get test file path
    Args:
        file_name: optional the name of the file to get the path to

    Returns:
        full path if a file is given or the path of the testfile directory is file_name is None
    """
    tfp = os.path.split(__file__)[0]
    if file_name is None:
        return tfp
    else:
        return os.path.join(tfp, file_name)


def get_sp(struc):
    """
    The path to a structure in the testfiles/structures folder.

    Args:
        struc (str): the name of the structure.

    Returns:
        str: the absolute path to the coord file.
    """
    return os.path.join(get_tfp(), 'structures', struc)


def get_control_integration(filename):
    """
    The path to a reference control file in the testfiles/integration/control folder.

    Args:
        filename (str): the name of the reference control file.

    Returns:
        str: the absolute path to the coord file.
    """
    return os.path.join(get_tfp(), 'integration', 'control', filename)


class ItestConfig:
    """
    Helper class to store the configuration parameters used inside the integration
    test function. The actual values are expected to be set by pytest.
    """
    define_timeout = 10
    generate_ref = False
    dryrun = False
    dryrun_fpath = 'dryrun_itest.json'
    tol = 1e-4
    delete_tmp_dir = True
    dryrun_use_ref_control = False


ignored_itest_parsed_keys = ["start_time", "end_time", "wall_time", "cpu_time", "construction_timings",
                             "host", "version", "build", "eigenvectors", "reduced_masses", "electric_dipole",
                             "magnetic_dipole", "electric_quadrupole", "@version"]


def run_itest(executables, define_options, coord_filename, control_reference_filename, file_classes,
              arguments=None, datagroups_options=None):
    """
    Runs the integration tests. First define is executed and the produced control file is compared
    with the reference. If successful the required turbomole executables are run and the numerical values
    of the outputs will be checked with the reference.

    Args:
        executables: string or list of strings identifying the list of programs that should be executed.
        define_options (dict): the options passed to DefineRunner.
        coord_filename (str): the filename of the coord used. Will be taken from the testfiles/structures folder.
        control_reference_filename (str): the name of the reference control file used to check the correctness of
            the execution of DefineRunner. Will be taken from the testfiles/integration/control folder.
        file_classes: a list of classes subclassing BaseData (or even a single one if only one program required)
            they will be used to generate outputs from stdout and compared with the reference.
        arguments: string or list of strings with the arguments to be passed to the each executable.
        datagroups_options (dict): a dict of the form {"datagroup_name": "datagroup_value"}. Can contain
            additional datagroups that will be set with the cdg utility before running the calculation.

    Returns:
        bool: True if the test passed successfully
    """

    if not isinstance(executables, (list, tuple)):
        executables = [executables]

    if not isinstance(file_classes, (list, tuple)):
        file_classes = [file_classes]

    if arguments is None:
        arguments = [None] * len(executables)
    elif not isinstance(arguments, (list, tuple)):
        arguments = [arguments]

    opt_define_timeout = ItestConfig.define_timeout
    opt_generate_ref = ItestConfig.generate_ref
    opt_dryrun = ItestConfig.dryrun
    dryrun_fpath = ItestConfig.dryrun_fpath
    opt_tol = ItestConfig.tol

    dryrun_differences = []

    with temp_dir(ItestConfig.delete_tmp_dir) as tmp_dir:
        # get the coord file (for the structure defined in the string)
        shutil.copyfile(get_sp(coord_filename), 'coord')

        dr = DefineRunner(define_options, timeout=opt_define_timeout)
        define_out = dr.run_full()

        # define should complete normally
        if not define_out:
            raise ItestError("Define did not complete normally.")

        if datagroups_options:
            c = Control.from_file()
            for k, v in datagroups_options.items():
                c.cdg(k, v)
            c.to_file()

        if opt_generate_ref:  # pragma: no cover
            shutil.copy2("control", get_control_integration(control_reference_filename))

        ref_control = Control.from_file(get_control_integration(control_reference_filename))

        current_control = Control.from_file("control")
        compare_control = current_control.compare(ref_control, tol=opt_tol)
        if ItestConfig.dryrun_use_ref_control:
            shutil.copy("control", "control_generated")
            shutil.copy(get_control_integration(control_reference_filename), "control")
        # print the output of Control.compare if the compare fails
        if opt_dryrun and compare_control is not None:
            dryrun_differences.append(('control', compare_control))
        else:
            assert compare_control is None, compare_control

        for iexec, (executable, exec_args, out_parser) in enumerate(zip(executables, arguments, file_classes)):
            cmd = [executable]
            if exec_args:
                cmd += shlex.split(exec_args)
            process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE, encoding='utf-8')
            program_std_out, program_std_err = process.communicate()

            try:
                ret_code = process.wait()

                # The riper executable does not echo "ended normally" in the std err file...
                if ret_code or (executable != "riper" and "ended normally" not in program_std_err):
                    raise ItestError("Executable {} has failed with return code {}".format(executable, ret_code))

                if out_parser:
                    # for jobex the outputs are not in the stdout but in the job.last file
                    if executable == "jobex":
                        out = out_parser.from_file("job.last").as_dict()
                    else:
                        out = out_parser.from_string(program_std_out).as_dict()
                    out_ref_path = os.path.join(get_tfp(), "integration", "logs_json",
                                                "{}_{}.json".format(control_reference_filename, executable))
                    if opt_generate_ref:
                        dumpfn(out, out_ref_path)
                    out_ref = loadfn(out_ref_path).as_dict()
                    if opt_dryrun:
                        diffs = compare_differences(out, out_ref, atol=opt_tol)
                        if diffs:
                            dryrun_differences.append((f'{executable} ({iexec})', diffs))
                    else:
                        assert_almost_equal(out, out_ref, atol=opt_tol, ignored_values=ignored_itest_parsed_keys)

                c = Control.from_file("control")
                e = c.energy
                if e is not None:
                    e_ref_path = os.path.join(get_tfp(), "integration", "energy",
                                              "{}_{}.json".format(control_reference_filename, executable))
                    if opt_generate_ref:
                        dumpfn(e, e_ref_path)
                    e_ref = loadfn(e_ref_path)
                    if opt_dryrun:
                        diffs = compare_differences(e.scf, e_ref.scf, atol=opt_tol)
                        if diffs:
                            dryrun_differences.append((f'{executable} ({iexec}) Energy.scf', diffs))
                        diffs = compare_differences(e.total, e_ref.total, atol=opt_tol)
                        if diffs:
                            dryrun_differences.append((f'{executable} ({iexec}) Energy.total', diffs))
                    else:
                        np.testing.assert_allclose(e.scf, e_ref.scf, atol=opt_tol)
                        np.testing.assert_allclose(e.total, e_ref.total, atol=opt_tol)

                g = c.gradient
                if g is not None:
                    g_ref_path = os.path.join(get_tfp(), "integration", "gradient",
                                              "{}_{}.json".format(control_reference_filename, executable))
                    if opt_generate_ref:
                        dumpfn(g, g_ref_path)
                    g_ref = loadfn(g_ref_path)
                    if opt_dryrun:
                        diffs = compare_differences(g.gradients, g_ref.gradients, atol=opt_tol)
                        if diffs:
                            dryrun_differences.append((f'{executable} ({iexec}) Gradient', diffs))
                    else:
                        np.testing.assert_allclose(g.gradients, g_ref.gradients, atol=opt_tol)

                # check that the output from eiger and our parser give the same results
                states = States.from_file()
                eiger_runner = EigerRunner()
                eiger_runner.run()
                eiger_out = eiger_runner.get_eiger_output()
                eiger_comp = eiger_out.compare_states(states)
                if opt_dryrun:
                    if eiger_comp is not None:
                        dryrun_differences.append((f'{executable} ({iexec}) Eiger comparison', diffs))
                else:
                    assert eiger_comp is None

                with open("{}_stdout".format(executable), "w") as f:
                    f.write(program_std_out)
            except:
                # if an exception is raised write down the output file for debugging then reraise
                with open("{}_stdout".format(executable), "w") as f:
                    f.write(program_std_out)
                with open("{}_stderr".format(executable), "w") as f:
                    f.write(program_std_err)

                raise

        if opt_dryrun:
            if dryrun_differences:
                frames = inspect.getouterframes(inspect.currentframe())
                test_frame = frames[1]
                fname_itest = test_frame.filename.split('/turbomoleio/')[-1].strip()
                funct_itest = test_frame.function
                line_itest = test_frame.lineno
                if not os.path.exists(dryrun_fpath):
                    alldiffs = []
                else:
                    with open(dryrun_fpath, 'r') as f:
                        alldiffs = json.load(f)
                diff_identifier = (fname_itest, funct_itest, coord_filename, control_reference_filename, line_itest)
                if diff_identifier in [tuple(did) for did, d in alldiffs]:
                    raise RuntimeError('Difference already in list. This might be that you are running the dry-run '
                                       'mode a second time. You should delete the dry-run file ("dryrun_itest.json") '
                                       'or specify a different one with the --dryrun-fpath=DRYRUN_FPATH option '
                                       'in pytest.')
                alldiffs.append((diff_identifier, dryrun_differences))
                with open(dryrun_fpath, 'w') as f:
                    json.dump(alldiffs, f, indent=2)

        return True


def generate_control_for_test(test_definition):
    """
    Regenerates control from the definition of the test.

    Args:
        test_definition: Dictionary with the definition of the test

    Returns:
        None
    """
    if test_definition['define']['template']:
        dp = get_define_template(test_definition['define']['template'])
    else:
        dp = {}
    if test_definition['define']['parameters']:
        dp.update(test_definition['define']['parameters'])
    datagroups = test_definition.get('datagroups')
    add_datagroups = None
    change_datagroups = None
    if datagroups:
        add_datagroups = test_definition.get('datagroups', {}).get('add', None)
        change_datagroups = test_definition.get('datagroups', {}).get('change', None)
    dr = DefineRunner(dp)
    run_full_ok = dr.run_full()
    if not run_full_ok:
        raise DefineError('Generation of control file with the DefineRunner failed.')
    if add_datagroups:
        for dg, val in add_datagroups.items():
            adg(dg, val)
    if change_datagroups:
        for dg, val in change_datagroups.items():
            cdg(dg, val)


def generate_reference_out_parser_files(log_fpath, outdir=None):
    """
    Helper function to generate the reference files for the test of the files and parser objects
    for a given log file.

    Args:
        log_fpath (str): Path to the logfile.
        outdir (str): Directory to output reference File and Parser_methods objects.
    """
    if outdir is None:
        outdir = '.'
    fname = os.path.split(log_fpath)[1]
    if fname == 'job.last':
        tm_exec = 'jobex'
    else:
        tm_exec = os.path.splitext(fname)[0]

    if tm_exec == 'jobex':
        out = JobexOutput.from_file(log_fpath)
    elif tm_exec == 'escf':
        out = exec_to_out_obj['escf'].from_file(log_fpath)
        out_escf_only = EscfOnlyOutput.from_file(log_fpath)
        dumpfn(out_escf_only, os.path.join(outdir, 'ref_escf_output.json'), indent=2)
    else:
        out = exec_to_out_obj[tm_exec].from_file(log_fpath)
    dumpfn(out, os.path.join(outdir, 'ref_output.json'), indent=2)

    parser = Parser.from_file(log_fpath)
    parsed_data = {}
    for m in PARSER_METHODS:
        data = getattr(parser, m)
        parsed_data[m] = data

    dumpfn(parsed_data, os.path.join(outdir, 'ref_parser.json'), indent=2)
    shutil.copy(log_fpath, outdir)


def generate_reference_output(test_definition,
                              ):
    """
    Executes the list of Turbomole commands in the definition of the test
    Args:
        test_definition: Dictionary with the definition of the test

    Returns:
        None
    """
    cmds = test_definition['commands']
    for cmd in cmds:
        cmd_split = cmd.split()
        tm_exec = cmd_split[0]
        options = ''
        if len(cmd_split) > 1:
            options = ' '.join(cmd_split[1:])
            options = f' {options}'
        os.system(f'{tm_exec}{options} > {tm_exec}.log 2> {tm_exec}.err')
