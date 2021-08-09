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
Script to generate reference output files for a new TM version and check differences.
"""


import argparse
import os
import sys
from turbomoleio.testfiles.utils import TM_VERSIONS
from monty.serialization import loadfn
from turbomoleio.testfiles.utils import TESTDIR
from turbomoleio.core.utils import get_tm_version
from turbomoleio.testfiles.utils import generate_control_for_test
from turbomoleio.testfiles.utils import generate_reference_output
from turbomoleio.testfiles.utils import generate_reference_out_parser_files
from turbomoleio.testfiles.utils import PARSER_METHODS
from turbomoleio.input.utils import get_define_template
from turbomoleio.input.define import DefineRunner
from turbomoleio.core.control import cpc
from turbomoleio.core.control import Control
from turbomoleio.output.files import exec_to_out_obj
from turbomoleio.output.files import EscfOnlyOutput
from turbomoleio.output.files import JobexOutput
from turbomoleio.output.parser import Parser
from turbomoleio.testfiles.utils import compare_differences
from monty.os import cd, makedirs_p
from monty.serialization import dumpfn
import shutil


gen_dir = os.path.join(TESTDIR, 'outputs', 'generation')
OUTPUTS_BASENAMES = loadfn(os.path.join(gen_dir, 'tests_config.yaml'))['testlist']
exec_to_out_obj = dict(exec_to_out_obj)
exec_to_out_obj['jobex'] = JobexOutput


def get_args(parser):
    """Get arguments for script.

    Args:
        parser: ArgumentParser object.

    Returns:
        Namespace with all the arguments from the parser.
    """
    parser.add_argument("--list", help="Print list of all tests.",
                        action="store_true", default=False)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--test", help="Run specific test. "
                                      "Can either be the name of a given executable (e.g. \"aoforce\"), "
                                      "in which case all tests related to aoforce are performed, "
                                      "or the name of a specific test (e.g. \"h2o_std\"), "
                                      "in which case all the tests with this name are performed, "
                                      "or the name of the executable and the name of the test "
                                      "(e.g. \"aoforce h2_numforce\"), "
                                      "in which case the specific test is performed.", nargs='+')
    parser.add_argument("--dryrun", help="Perform a dry run of the tests (default is to generate the files). "
                                         "A diff file is generated based on the last TM version.",
                        action="store_true", default=False)
    parser.add_argument("--dryrun_fname", help="Name of the file with the differences.",
                        type=str, default="differences.json")
    parser.add_argument("--compare_to", help="Version directory to compare to when executing in dryrun mode.",
                        type=str, default=None)
    parser.add_argument("--keep_rundirs", help="Whether to keep the run directories used for the generation "
                                               "(default is: do not keep the run directories).",
                        action="store_true", default=False)
    parser.add_argument("--print_diffs", help="Whether to print the differences to stdout.",
                        action="store_true", default=False)
    parser.add_argument("--generate_control", help="Whether to regenerate the control file "
                                                   "(default is: do not regenerate the control file).",
                        action="store_true", default=False)
    parser.add_argument("--only_control", help="Whether to only regenerate the control file and not the output files "
                                               "(default is to regenerate both the control and the output files).",
                        action="store_true", default=False)
    parser.add_argument("--no_gen_control_update",
                        help="Do not update and backup coord, control and "
                             "other related files in the generation folder.",
                        action="store_true", default=False)
    parser.add_argument("--check_update_json_files",
                        help="Checks and updates all reference parser and output json files. This is needed "
                             "when additional data is parsed (e.g. new parsing method in the parser, more data "
                             "in a Data or File object, ...).",
                        action="store_true", default=False)
    parser.add_argument("--force", help="Overwrite existing files if any (default is: "
                                        "do not overwrite existing files).",
                        action="store_true", default=False)
    parser.add_argument("--version_dir", help="Specify version directory name for generated files.",
                        type=str, default=None)
    parser.add_argument("--atol", help="Absolute tolerance for comparison.",
                        type=float, default=0.0)
    parser.add_argument("--rtol", help="Relative tolerance for comparison.",
                        type=float, default=1e-7)
    args = parser.parse_args()
    if args.only_control and not args.generate_control:
        parser.error(f'"only_control" parameter can only be True if "generate_control" is True.')
    if args.compare_to is not None and not args.dryrun:
        parser.error(f'"compare_to" can only be used in dry run mode.')
    return args


def print_tests():
    """Print a list of all the tests."""
    for tm_exec, test_names in OUTPUTS_BASENAMES.items():
        print(tm_exec, test_names)


def get_tests_list(test, parser):
    """Get the list of tests to be performed.

    Args:
        test: String identifying the test or None. If None, all tests are returned.
        parser: ArgumentParser object. Used to signal a parser error.

    Returns:
        list: List of the tests to be performed.
    """
    tests_list = []
    if test:
        if len(test) == 1:
            name = test[0]
            if name in OUTPUTS_BASENAMES:
                tests_list = [(name, test) for test in OUTPUTS_BASENAMES[name]]
            else:
                tests_list = [
                    (tm_exec, test)
                    for tm_exec, tests in OUTPUTS_BASENAMES.items()
                    for test in tests if name == test
                ]
        elif len(test) == 2:
            tm_exec, tname = test
            if tm_exec not in OUTPUTS_BASENAMES:
                parser.error(f'TM exec "{tm_exec}" not found.')
            if tname not in OUTPUTS_BASENAMES[tm_exec]:
                parser.error(f'Test "{tm_exec} {tname}" not found.')
            tests_list = [(tm_exec, tname)]
        else:
            parser.error(
                f'argument --test accepts either 1 or 2 values ({len(test)} given).'
            )
    else:
        tests_list = [
            (tm_exec, test) for tm_exec, tests in OUTPUTS_BASENAMES.items()
            for test in tests
        ]

    if not tests_list:
        parser.error('No test found, check script options (--test) and list of tests (--list)')
    return tests_list


def get_version_dir(version_dir, force, compare_to):
    """Get the version directories: current version used and reference version.

    Args:
        version_dir: Explicitly defined directory for the tests. Default directory is set up if it is None.
        force: Whether to force overwriting of existing test directories.
        compare_to: To which version directory the current generated tests should be compared. Default is None, i.e.
            compare to the previous version.

    Returns:
        tuple: Directory for the tests with the current version and directory of the reference version.
    """
    tm_version = get_tm_version()
    print(f'Turbomole version {tm_version} detected')
    if not version_dir:
        version_dir = f'TM_v{tm_version}'
    vdir_path = os.path.join(TESTDIR, 'outputs', version_dir)
    if not force and os.path.exists(vdir_path):
        print('Directory exists. If generation of existing tests is performed, code will exit. '
              'Use --force to overwrite.')
    print(f'Version directory will be <{version_dir}>.\n'
          f'New reference output files will be generated in <testfiles/outputs/{version_dir}>.')
    if compare_to is None:
        ref_vdir_path = os.path.join(TESTDIR, 'outputs', TM_VERSIONS[-1])
    else:
        ref_vdir_path = os.path.join(TESTDIR, 'outputs', compare_to)
    if not os.path.exists(ref_vdir_path):
        print('Reference version directory does not exist')
        exit()
    return vdir_path, ref_vdir_path


def get_paths_and_deftest(force, version_dir_path, ref_version_dir_path, tm_exec, test_name):
    """Set up the directories and get the deftest config for the test tm_exec/test_name.

    Args:
        force: Whether to force overwriting of existing test directories.
        version_dir_path: Absolute path of the current version directory.
        ref_version_dir_path: Absolute path of the reference version directory.
        tm_exec: String identifying the Turbomole executable being tested.
        test_name: String identifying the name of the test.

    Returns:
        tuple: Absolute path of the test directory (where new reference files should be copied),
            of the run directory (where Turbomole is run before output files are compared and then
            copied to the test directory), of the reference test directory (where reference files
            are stored), of the generation test directory (where information about generation is
            stored), and configuration of the test as dictionary.
    """
    test_dir_path = os.path.join(version_dir_path, tm_exec, test_name)
    if not force and os.path.exists(test_dir_path):
        print('Test directory exists. Use --force to overwrite.')
        exit()
    test_run_dir_path = os.path.join(test_dir_path, 'run')
    ref_test_dir_path = os.path.join(ref_version_dir_path, tm_exec, test_name)
    gen_test_dir_path = os.path.join(gen_dir, tm_exec, test_name)
    deftest_fpath = os.path.join(gen_test_dir_path, "test.yaml")
    if not os.path.isfile(deftest_fpath):
        raise RuntimeError(f'No test.yaml file for test {tm_exec}/{test_name}')
    deftest = loadfn(deftest_fpath)
    return test_dir_path, test_run_dir_path, ref_test_dir_path, gen_test_dir_path, deftest


def generate_control(print_diffs, deftest, gen_test_dir, test_run_dir, all_diffs):
    """
    Generate the control file and related files and compare the generated control file with the reference
    control file from the generation directory.

    Args:
        print_diffs: Whether to print the differences between the generated and reference control files.
        deftest: Configuration of the test specifying how to run define (i.e. which define template/parameters to use).
        gen_test_dir: Generation test directory (where information about generation is stored).
        test_run_dir: Run test directory (where Turbomole is run before output files are compared and then
            copied to the test directory).
        all_diffs: Dictionary containing the differences for this test. Will be updated if differences are found
            between the new and the reference control files.
    """
    # Some tests use a fixed control file
    if deftest.get('fixed_control', False):
        generate_mos(deftest, test_run_dir, gen_test_dir)
    elif deftest.get('define', None) is None:
        raise ValueError('No define template and/or parameters provided for reference test generation.')
    else:
        generate_control_for_test(test_definition=deftest)
        ref_control_fpath = os.path.join(
            gen_test_dir, deftest['control']
            if 'control' in deftest else 'control'
        )
        if os.path.exists(ref_control_fpath):
            ref_control = Control.from_file(ref_control_fpath)
        else:
            ref_control = Control.empty()
        test_control = Control.from_file(os.path.join(test_run_dir, 'control'))
        control_diffs = test_control.compare(ref_control, return_all_diffs=True)
        if control_diffs:
            if print_diffs:
                msg = 'There are differences in the control files generated:\n'
                for idiff, diff in enumerate(control_diffs, start=1):
                    msg += f'#{idiff} {diff}\n'
                print(msg)
            all_diffs['control'] = control_diffs


def get_coord(test_definition, rundir, gendir):
    """
    Copy the coord file to the run directory of the test

    Args:
        test_definition: Configuration of the test
        rundir: Directory where the test is run.
        gendir: Directory where the generation information is stored.
    """
    if 'coord' in test_definition:
        coord_fpath = os.path.join(gendir, test_definition['coord'])
    else:
        coord_fpath = os.path.join(gendir, 'coord')

    shutil.copy(coord_fpath, os.path.join(rundir, 'coord'))


def cpc_backup_previous_control(gen_test_dir, compare_to, no_gen_control_update):
    """
    Back up the coord, control and basis files of the previous Turbomole reference folder and copy the new
    ones to the generation folder.

    Note:
        This should be called from the corresponding test run directory.

    Args:
        gen_test_dir: Generation test directory where to copy the new coord, control and basis files.
        compare_to: Reference version for which a backup is needed. If None, the previous version directory is used.
        no_gen_control_update: Whether to actually copy the new coord, control and basis files to the generation
            directory.
    """
    # Keep a copy of the coord, control and basis files of the previous Turbomole reference folder
    if not no_gen_control_update:

        if compare_to is None:
            vdir = TM_VERSIONS[-1]
        else:
            vdir = compare_to
        prev_gen_control_dir = os.path.join(gen_test_dir, vdir)
        if os.path.exists(os.path.join(gen_test_dir, 'control')):
            # This is to back up the control, coord and basis files used for the previous version
            # Not performed for a completely new test
            if os.path.exists(prev_gen_control_dir):
                print(f' ... control, coord and basis files for {vdir} already backed up in generation folder, '
                      f'will not overwrite!')
            else:
                makedirs_p(prev_gen_control_dir)
                cpc(prev_gen_control_dir, control_dir=gen_test_dir)
        # Copy the coord, control and basis files of the current Turbomole version to the
        # generation directory
        for fname in ['control', 'coord', 'basis', 'auxbasis']:
            if os.path.exists(fname):
                shutil.copy(fname, gen_test_dir)


def generate_mos(deftest, test_run_dir, gen_test_dir):
    """
    Generate initial MO files (mos or alpha/beta).

    This uses the existing control and related files in the generation directory of the test.
    First a "cpc" (copy control) is performed from gen_test_dir to test_run_dir, the run
    directory for the test.

    Args:
        deftest: Configuration of the test.
        test_run_dir: Directory where the test is run.
        gen_test_dir: Directory where the generation information is stored. The control and related files
            are also stored in this directory.

    Returns:
        None
    """
    if deftest['define'] and deftest['define']['template']:
        dp = get_define_template(deftest['define']['template'])
    else:
        dp = {}
    if deftest['define'] and deftest['define']['parameters']:
        dp.update(deftest['define']['parameters'])
    cpc(test_run_dir, force_overwrite=False, control_dir=gen_test_dir)
    dr = DefineRunner(dp)
    dr.run_generate_mo_files()


def get_log_fpath(test_run_dir, tm_exec):
    """Get the absolute path of the log file.

    Args:
        test_run_dir: Absolute path to the run test directory.
        tm_exec: Turbomole executable being tested.

    Returns:
        Absolute path to the log file of the run.
    """
    if tm_exec == 'jobex':
        return os.path.join(test_run_dir, 'job.last')
    else:
        return os.path.join(test_run_dir, f'{tm_exec}.log')


def compare_to_reference_files(log_fpath, ref_test_dirpath, out_cls, rtol, atol, print_diffs, all_diffs,
                               tm_exec, ref_fname='ref_output.json'):
    """
    Generate the output object (and its serialized version) and compare it to the references one.

    Args:
        log_fpath: Absolute path to the generated log file.
        ref_test_dirpath: Absolute path of the reference test directory.
        out_cls: Output object used to parse the log file.
        rtol: relative tolerance.
        atol: absolute tolerance.
        print_diffs: Whether to print the differences found.
        all_diffs: Dictionary containing the differences for this test. Will be updated if differences are found
            between the new and the reference log files.
        tm_exec: Turbomole executable being tested.
        ref_fname: Filename of the json-serialized reference output object.

    Returns:
        tuple: Output object and serialized Output object.
    """
    ref_fpath = os.path.join(ref_test_dirpath, ref_fname)
    if os.path.exists(ref_fpath):
        # Comparison below is performed on the dictionary directly, not on the object, hence the cls=None
        ref_out = loadfn(ref_fpath, cls=None)
    else:
        # When creating a new test, the reference object does not exist yet
        ref_out = None

    gen_out = out_cls.from_file(log_fpath)
    gen_out_dict = gen_out.as_dict()

    out_diffs = compare_differences(gen_out_dict, ref_out, rtol=rtol, atol=atol)
    if out_diffs:
        if print_diffs:
            print(f'There are differences in the parsed {outurbomoleio/testfiles/utils.pyt_cls.__name__} output file objects:')
            for idiff, diff in enumerate(out_diffs, start=1):
                print(f'#{idiff} {diff[0]}\n  {diff[1]}')
        all_diffs[f'{tm_exec}_{out_cls.__name__}'] = out_diffs

    return gen_out, gen_out_dict


def compare_to_reference_parser(log_fpath, ref_test_dirpath, rtol, atol, print_diffs, all_diffs, tm_exec):
    """
    Generate all results of parser methods and compare to the reference parser file.

    Args:
        log_fpath: Absolute path to the generated log file.
        ref_test_dirpath: Absolute path of the reference test directory.
        rtol: relative tolerance.
        atol: absolute tolerance.
        print_diffs: Whether to print the differences found.
        all_diffs: Dictionary containing the differences for this test. Will be updated if differences are found
            between the new and the reference parser files.
        tm_exec: Turbomole executable being tested.

    Returns:
        dict: Dictionary mapping each parser method with it's resulting parsed data.
    """
    ref_fpath = os.path.join(ref_test_dirpath, 'ref_parser.json')
    if os.path.exists(ref_fpath):
        # Comparison below is performed on the dictionary directly, not on the object, hence the cls=None
        ref_parser_methods = loadfn(ref_fpath, cls=None)
    else:
        # When creating a new test, the reference object does not exist yet
        ref_parser_methods = None
    tmio_parser = Parser.from_file(log_fpath)
    parsed_data = {}
    for m in PARSER_METHODS:
        data = getattr(tmio_parser, m)
        parsed_data[m] = data

    # Compare the parsed data with the reference
    parser_diffs = compare_differences(parsed_data, ref_parser_methods, rtol=rtol, atol=atol)
    if parser_diffs:
        if print_diffs:
            print(f'There are differences in the parsed methods:')
            for idiff, diff in enumerate(parser_diffs, start=1):
                print(f'#{idiff} {diff[0]}\n  {diff[1]}')
        all_diffs[f'{tm_exec}_Parser'] = parser_diffs

    return parsed_data


def get_outputs_and_parser_objects_and_compare(log_fpath, ref_test_dirpath,
                                               rtol, atol, print_diffs, all_diffs, tm_exec):
    """
    Get the serialized

    Args:
        log_fpath: Absolute path to the generated log file.
        ref_test_dirpath: Absolute path of the reference test directory.
        rtol: relative tolerance.
        atol: absolute tolerance.
        print_diffs: Whether to print the differences found.
        all_diffs: Dictionary containing the differences for this test.
            Will be updated if differences are found
            between the new and the reference parser and/or output files.
        tm_exec: Turbomole executable being tested.

    Returns:
        dict: Dictionary mapping the reference outputs and parser filenames
            to the serialized outputs and parser objects.
    """
    out_cls = exec_to_out_obj[tm_exec]
    output, output_dict = compare_to_reference_files(log_fpath, ref_test_dirpath, out_cls, rtol, atol,
                                                     print_diffs, all_diffs,
                                                     tm_exec, ref_fname='ref_output.json')
    outputs_parser_dicts = {'ref_output.json': output_dict}
    if tm_exec == 'escf':
        escf_output, escf_output_dict = compare_to_reference_files(
            log_fpath, ref_test_dirpath, EscfOnlyOutput, rtol, atol,
            print_diffs, all_diffs,
            tm_exec, ref_fname='ref_escf_output.json'
        )
        outputs_parser_dicts['ref_escf_output.json'] = escf_output_dict

    outputs_parser_dicts['ref_parser.json'] = compare_to_reference_parser(
        log_fpath, ref_test_dirpath, rtol, atol,
        print_diffs, all_diffs, tm_exec
    )

    return outputs_parser_dicts


def update_json_files(dryrun, rtol, atol, print_diffs):
    """
    Updates all json files of all versions when the parser and/or the data/file objects are modified or extended.

    Args:
        dryrun: Whether to perform a dryrun of the generation and dump the differences between the previous reference
            files with the serialized parser methods results and file objects.
        rtol: relative tolerance.
        atol: absolute tolerance.
        print_diffs: Whether to print the differences found.
    """
    print("Updating json files")
    for tm_version in TM_VERSIONS:
        vdir_path = os.path.join(TESTDIR, 'outputs', tm_version)
        vdir_testlist = loadfn(os.path.join(vdir_path, 'tests_config.yaml'))['testlist']
        version_all_diffs = {}
        for tm_exec, testnames in vdir_testlist.items():
            for test_name in testnames:
                all_diffs = {}
                test_dirpath = os.path.join(vdir_path, tm_exec, test_name)
                # Get the new the serialized output and parser objects and compare to the old ones
                log_fpath = get_log_fpath(test_dirpath, tm_exec)
                outputs_parser_dicts = get_outputs_and_parser_objects_and_compare(
                    log_fpath, test_dirpath,
                    rtol, atol, print_diffs, all_diffs, tm_exec
                )

                # Dump the new json-serialized reference output and parser objects
                if not dryrun:
                    for fname, ref_dict in outputs_parser_dicts.items():
                        dumpfn(ref_dict, os.path.join(test_dirpath, fname), indent=2)
                # Add the differences found in this test to the full list (dictionary) of differences in this version
                if tm_exec not in version_all_diffs:
                    version_all_diffs[tm_exec] = {}
                if test_name in version_all_diffs[tm_exec]:
                    raise RuntimeError('Test already in the list of differences')
                version_all_diffs[tm_exec][test_name] = all_diffs

        # Dump the file containing the differences found in the tests for this version
        if dryrun:
            dumpfn(version_all_diffs, os.path.join(vdir_path, 'parsing_update_differences.json'), indent=2)


def main():
    """Main function"""
    parser = argparse.ArgumentParser()
    args = get_args(parser)

    if args.list:
        return print_tests()

    tests_list = get_tests_list(args.test, parser)
    dryrun = args.dryrun

    if args.check_update_json_files:
        return update_json_files(dryrun, args.rtol, args.atol, args.print_diffs)

    version_dirpath, ref_version_dirpath = get_version_dir(args.version_dir, args.force, args.compare_to)

    all_tests_diffs = {}

    for tm_exec, test_name in tests_list:
        # Set up test directories and get deftest configuration
        test_dirpath, test_run_dirpath, ref_test_dirpath, gen_test_dirpath, deftest = get_paths_and_deftest(
            args.force, version_dirpath, ref_version_dirpath, tm_exec, test_name
        )
        if deftest.get('disable', False):
            print(f'!Generation of outputs for test {tm_exec}/{test_name} is disabled')
            continue

        print(f'Generation of outputs for test {tm_exec}/{test_name}')

        # Create run directory for test and copy coord file
        makedirs_p(test_run_dirpath)
        get_coord(test_definition=deftest, rundir=test_run_dirpath, gendir=gen_test_dirpath)

        with cd(test_run_dirpath):
            all_diffs = {}
            # Generate the control file
            if args.generate_control:
                generate_control(args.print_diffs, deftest, gen_test_dirpath, test_run_dirpath, all_diffs)
            # Copy the control file and generate the initial mos
            else:
                generate_mos(deftest, test_run_dirpath, gen_test_dirpath)

            # Run the Turbomole executables
            if not args.only_control:
                if not dryrun:
                    # Backup the previous coord, control and basis files and copy the new ones to the generation test
                    # directory
                    cpc_backup_previous_control(gen_test_dirpath, args.compare_to, args.no_gen_control_update)

                # Run Turbomole and get the new log files.
                generate_reference_output(test_definition=deftest)

                # Compare the generated log files with the references using the serialized output and parser
                # objects
                log_fpath = get_log_fpath(test_run_dirpath, tm_exec)
                outputs_parser_dicts = get_outputs_and_parser_objects_and_compare(
                    log_fpath, ref_test_dirpath,
                    args.rtol, args.atol, args.print_diffs, all_diffs, tm_exec
                )

                # Dump the new json-serialized reference output and parser objects
                if not dryrun:
                    for fname, ref_dict in outputs_parser_dicts.items():
                        dumpfn(ref_dict, os.path.join(test_dirpath, fname), indent=2)
                    shutil.copy(log_fpath, test_dirpath)

        # Remove the run directories
        if not args.keep_rundirs:
            shutil.rmtree(test_run_dirpath)

        # Add the differences found in this test to the full list (dictionary) of differences
        if tm_exec not in all_tests_diffs:
            all_tests_diffs[tm_exec] = {}
        if test_name in all_tests_diffs[tm_exec]:
            raise RuntimeError('Test already in the list of differences')
        all_tests_diffs[tm_exec][test_name] = all_diffs

    # Dump the file containing the differences found in the tests
    if dryrun:
        dumpfn(all_tests_diffs, os.path.join(version_dirpath, args.dryrun_fname), indent=2)


if __name__ == "__main__":
    sys.exit(main())
