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
from turbomoleio.testfiles.utils import OUTPUTS_BASENAMES
from turbomoleio.testfiles.utils import TM_VERSIONS
from monty.serialization import loadfn
from turbomoleio.testfiles.utils import TESTDIR
from turbomoleio.core.utils import get_tm_version
from turbomoleio.testfiles.utils import generate_control_for_test
from turbomoleio.testfiles.utils import generate_reference_output
from turbomoleio.input.utils import get_define_template
from turbomoleio.input.define import DefineRunner
from turbomoleio.core.control import cpc
from turbomoleio.core.control import Control
from turbomoleio.output.files import exec_to_out_obj
from turbomoleio.output.files import JobexOutput
from turbomoleio.testfiles.utils import compare_differences
from monty.os import cd, makedirs_p
import shutil


gen_dir = os.path.join(TESTDIR, 'outputs', 'generation')


def get_args(parser):
    """Get arguments for script."""
    parser.add_argument("--list", help="Print list of all tests.",
                        action="store_true", default=False)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--all", help="Run generation for all tests.",
                       action="store_true", default=True)
    group.add_argument("--test", help="Run specific test.", nargs='+')
    parser.add_argument("--dryrun", help="Perform a dry run of the tests."
                                         "A diff file is generated based on the last TM version",
                        action="store_true", default=True)
    parser.add_argument("--compare_to", help="Version directory to compare to when executing in dryrun mode.",
                        type=str, default=None)
    parser.add_argument("--generate_control", help="Whether to regenerate the control file.",
                        action="store_true", default=False)
    parser.add_argument("--only_control", help="Whether to regenerate the control file.",
                        action="store_true", default=False)
    parser.add_argument("--force", help="Overwrite existing files if any.",
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
    return args


def print_tests():
    """Print a list of all the tests."""
    for tm_exec, test_names in OUTPUTS_BASENAMES.items():
        print(tm_exec, test_names)


def get_tests_list(args, parser):
    """Get the list of tests to be performed."""
    tests_list = []
    if args.test:
        args.all = False
        if len(args.test) == 1:
            name = args.test[0]
            if name in OUTPUTS_BASENAMES:
                tests_list = [(name, test) for test in OUTPUTS_BASENAMES[name]]
            else:
                tests_list = [
                    (tm_exec, test)
                    for tm_exec, tests in OUTPUTS_BASENAMES.items()
                    for test in tests if name == test
                ]
        elif len(args.test) == 2:
            tm_exec, tname = args.test
            if tm_exec not in OUTPUTS_BASENAMES:
                parser.error(f'TM exec "{tm_exec}" not found.')
            if tname not in OUTPUTS_BASENAMES[tm_exec]:
                parser.error(f'Test "{tm_exec} {tname}" not found.')
            tests_list = [(tm_exec, tname)]
        else:
            parser.error(
                f'argument --test accepts either 1 or 2 values ({len(args.test)} given).'
            )
    elif args.all:
        tests_list = [
            (tm_exec, test) for tm_exec, tests in OUTPUTS_BASENAMES.items()
            for test in tests
        ]
    else:
        raise RuntimeError('No test given')
    if not tests_list:
        parser.error('No test found, check script options (--all and --test) and list of tests (--list)')
    return tests_list


def get_version_dir(args):
    """Get the version directory"""
    tm_version = get_tm_version()
    print(f'Turbomole version {tm_version} detected')
    version_dir = args.version_dir
    if not version_dir:
        version_dir = f'TM_v{tm_version}'
    vdir = os.path.join(TESTDIR, 'outputs', version_dir)
    if not args.force and os.path.exists(vdir):
        print('Directory exists. Use --force to overwrite.')
        exit()
    print(f'Version directory will be <{version_dir}>.\n'
          f'New reference output files will be generated in <testfiles/outputs/{version_dir}>.')
    if args.compare_to is None:
        ref_version_dir = os.path.join(TESTDIR, 'outputs', TM_VERSIONS[-1])
    else:
        ref_version_dir = os.path.join(TESTDIR, 'outputs', args.compare_to)
    if not os.path.exists(ref_version_dir):
        print('Reference version directory does not exist')
        exit()
    return vdir, ref_version_dir


def get_coord(test_definition, rundir, gendir):
    if 'coord' in test_definition:
        coord_fpath = os.path.join(gendir, test_definition['coord'])
    else:
        coord_fpath = os.path.join(gendir, 'coord')

    shutil.copy(coord_fpath, os.path.join(rundir, 'coord'))


exec_to_out_obj = dict(exec_to_out_obj)
exec_to_out_obj['jobex'] = JobexOutput


def main():
    """Main function"""
    parser = argparse.ArgumentParser()
    args = get_args(parser)

    if args.list:
        return print_tests()

    tests_list = get_tests_list(args, parser)
    dryrun = args.dryrun
    version_dir, ref_version_dir = get_version_dir(args)

    for tm_exec, test_name in tests_list:
        test_dir = os.path.join(version_dir, tm_exec, test_name)
        test_run_dir = os.path.join(test_dir, 'run')
        ref_test_dir = os.path.join(ref_version_dir, tm_exec, test_name)
        gen_test_dir = os.path.join(gen_dir, tm_exec, test_name)
        deftest_fpath = os.path.join(gen_test_dir, "test.yaml")
        if not os.path.isfile(deftest_fpath):
            raise RuntimeError(f'No test.yaml file for test {tm_exec}/{test_name}')
        deftest = loadfn(deftest_fpath)
        if deftest.get('disable', False):
            print(f'Generation of outputs for test {tm_exec}/{test_name} is disabled')
            continue

        # Create run directory for test and copy coord file
        makedirs_p(test_run_dir)
        get_coord(test_definition=deftest, rundir=test_run_dir, gendir=gen_test_dir)

        with cd(test_run_dir):
            if args.generate_control:
                if deftest['define'] is None:
                    raise ValueError('Cannot generate control file from define parameters (control file is fixed).')
                generate_control_for_test(test_definition=deftest)
                if dryrun:
                    ref_control_fpath = os.path.join(
                        gen_test_dir, deftest['control']
                        if 'control' in deftest else 'control'
                    )
                    ref_control = Control.from_file(ref_control_fpath)
                    test_control = Control.from_file(os.path.join(test_run_dir, 'control'))
                    control_diffs = test_control.compare(ref_control, return_all_diffs=True)
                    if control_diffs:
                        print('There are differences in the control files generated:')
                        for idiff, diff in enumerate(control_diffs, start=1):
                            print(f'#{idiff} {diff}')
            else:
                if deftest['define'] and deftest['define']['template']:
                    dp = get_define_template('ridft')
                else:
                    dp = {}
                if deftest['define'] and deftest['define']['parameters']:
                    dp.update(deftest['define']['parameters'])
                print('Copying control file and related files')
                cpc(test_run_dir, force_overwrite=False, control_dir=gen_test_dir)
                dr = DefineRunner(dp)
                dr.run_generate_mo_files()

            if not args.only_control:
                generate_reference_output(test_definition=deftest)

                if dryrun:
                    ref_fpath = os.path.join(ref_test_dir, 'ref_output.json')
                    ref_out = loadfn(ref_fpath, cls=None)
                    out_cls = exec_to_out_obj[tm_exec]
                    if tm_exec == 'jobex':
                        log_fpath = os.path.join(test_run_dir, 'job.last')
                    else:
                        log_fpath = os.path.join(test_run_dir, f'{tm_exec}.log')
                    gen_out = out_cls.from_file(log_fpath).as_dict()

                    diffs = compare_differences(gen_out, ref_out, rtol=args.rtol, atol=args.atol)
                    if diffs:
                        print('There are differences in the parsed output file objects:')
                        for idiff, diff in enumerate(diffs, start=1):
                            print(f'#{idiff} {diff[0]}\n  {diff[1]}')


if __name__ == "__main__":
    sys.exit(main())



# from turbomoleio.core.utils import get_tm_version
# import os
# import sys
# from turbomoleio.testfiles.utils import OUTPUTS_BASENAMES
# from turbomoleio.testfiles.utils import generate_reference_output
# from monty.serialization import loadfn
# from turbomoleio.output.files import exec_to_out_obj
# from turbomoleio.output.files import JobexOutput
# from turbomoleio.testfiles.utils import TM_VERSIONS
# from turbomoleio.testfiles.utils import compare_differences
# # from turbomoleio.output.tests.test_files import generate_files as generate_files_refs
# # from turbomoleio.output.tests.test_parser import generate_files as generate_parser_refs
# import shutil
# import argparse
#
#
# # parser = argparse.ArgumentParser()
# # parser.parse_args()
# # parser.add_argument("--all", help="Run generation for all tests.",
# #                     action="store_true", default=True)
# # parser.add_argument("--dryrun", help="Perform a dry run of the tests."
# #                                      "A diff file is generated based on the last TM version",
# #                     action="store_true", default=True)
#
#
# tm_version = get_tm_version()
# print(f'Turbomole version {tm_version} detected')
# version_dir = input(f'Enter name of generated version directory (default: <TM_v{tm_version}>): ')
# if not version_dir:
#     version_dir = f'TM_v{tm_version}'
#
# print(f'Version directory will be <{version_dir}>.\n'
#       f'New reference output files will be generated in <testfiles/outputs/{version_dir}>.')
#
# outputs_path = os.path.join('outputs', version_dir)
# # if os.path.exists(outputs_path):
# #     print(f'The {outputs_path} already exists.\n'
# #           f'Either remove or change name of generated version directory.')
# #     exit()
#
# while True:
#     test = input('Which output files would you like to generate ?\n'
#                  ' - Enter <a> to generate all\n'
#                  ' - Enter <l> to list the output files\n'
#                  ' - Enter index of the output file to generate\n'
#                  ' - Enter <q> to quit\n')
#     if test == 'l':
#         print('Output files:')
#         for ioutput, (dirname, filebase) in enumerate(OUTPUTS_FILES_LIST):
#             print(f'#{ioutput} {dirname} : {filebase}')
#         continue
#     if test == 'q':
#         exit()
#     if test == 'a':
#         files = OUTPUTS_FILES_LIST
#         break
#     else:
#         try:
#             idx = int(test)
#         except ValueError:
#             print('Wrong key')
#             continue
#         files = [OUTPUTS_FILES_LIST[idx]]
#         break
#
#
# test = input(f'If you want to compare to previous outputs,\n'
#              f' - Enter <l> to compare with last outputs ({TM_VERSIONS[-1]})\n'
#              f' - Enter <n> if you do not want to compare\n')
#
# compare = test != 'n'
#
#
# for dirname, filebase in files:
#     refgendir = os.path.join('outputs', 'generation', dirname)
#     gendir = os.path.join(outputs_path, dirname)
#     rundir = os.path.join(gendir, f'run_{filebase}')
#     coord_fpath = os.path.join(refgendir, f'{filebase}.coord')
#     if not os.path.exists(coord_fpath):
#         continue
#     dp_dg_fpath = os.path.join(refgendir, f'{filebase}_dp_dg.json')
#     if not os.path.exists(dp_dg_fpath):
#         continue
#     tm_execs_fpath = os.path.join(refgendir, f'{filebase}_tm_execs.json')
#     if not os.path.exists(tm_execs_fpath):
#         continue
#     control_fpath = os.path.join(refgendir, f'{filebase}.control')
#     if os.path.exists(control_fpath):
#         control_file = control_fpath
#     else:
#         control_file = None
#     dp_dg = loadfn(dp_dg_fpath)
#     dp = dp_dg['dp']
#     adg = dp_dg['adg']
#     cdg = dp_dg['cdg']
#     tm_execs = loadfn(tm_execs_fpath)
#     generate_reference_output(coord_file=coord_fpath, dp=dp,
#                               tm_execs=tm_execs,
#                               rundir=rundir,
#                               control_file=control_file,
#                               add_datagroups=adg, change_datagroups=cdg,
#                               )
#
#     if compare:
#         comparedir = os.path.join('outputs', TM_VERSIONS[-1], dirname)
#         ref_fpath = os.path.join(comparedir, f'{filebase}_outfile.json')
#         ref_out = loadfn(ref_fpath, cls=None)
#         out_cls = exec_to_out_obj[tm_execs[-1]]
#         log_fpath = os.path.join(rundir, f'{tm_execs[-1]}.log')
#         gen_out = out_cls.from_file(log_fpath).as_dict()
#         diffs = compare_differences(gen_out, ref_out)
#
#         for diff, diffdescr in diffs:
#             print(diff)
#             print(diffdescr)
#
#         new_log_fpath = os.path.join(gendir, f'{filebase}.log')
#         copy = input(f'Enter <y> to copy <{log_fpath}> to <{new_log_fpath}> : ')
#         if copy == 'y':
#             shutil.copy(log_fpath, new_log_fpath)
#     # print(file)


