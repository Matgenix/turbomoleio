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

import os
import pytest
import tempfile
import shutil

from turbomoleio.core.datagroups import cleanup_string
from turbomoleio.core.datagroups import DataGroups
from turbomoleio.core.datagroups import remove_comments
from turbomoleio.core.datagroups import remove_dg_from_list
from turbomoleio.core.datagroups import split_string_to_dg_list, compare_datagroup_string
from turbomoleio.testfiles.utils import temp_dir, touch_file


def test_cleanup():
    """Testing cleanup function."""
    # Testing : No $ sign in the string
    sample_string = 'nodollarsign'
    with pytest.raises(expected_exception=ValueError,
                       match=r'^No dollar \("\$"\) sign in the string.$'):
        cleanup_string(string=sample_string,
                       cleanup_types=None)

    # Testing : Invalid type of cleanup
    sample_string = '$title'
    with pytest.raises(expected_exception=ValueError,
                       match=r'^Cleanup of type "UNKNOWN" is not a valid '
                             r'type.$'):
        cleanup_string(string=sample_string,
                       cleanup_types=['UNKNOWN'])

    # Testing : Characters preceding a $ sign
    sample_string = '$title\na$dummy'
    with pytest.raises(expected_exception=ValueError,
                       match=r'^Some character\(s\) are preceding a '
                             r'dollar \("\$"\) sign.$'):
        cleanup_string(string=sample_string,
                       cleanup_types=None)

    # Testing : inline comments are acceptable
    sample_string = '$title\na#this is a comment line\n$end'
    new_string = cleanup_string(string=sample_string,
                                cleanup_types=None)
    assert new_string == sample_string

    # Testing : Removing of lines before the first dollar sign
    sample_string = 'these two first lines\n' \
                    'should be removed\n' \
                    '$first dollar sign\n' \
                    '$second dollar sign\n' \
                    '$end'
    new_string = cleanup_string(string=sample_string,
                                cleanup_types=['BEFORE_FIRST_DOLLAR'])
    assert new_string == '$first dollar sign\n$second dollar sign\n$end'

    # Testing : Removing blank lines
    sample_string = '1\n2\n3\n' \
                    '$dummy abcd\n  \n' \
                    '$key2 var2\n' \
                    '    data\n\n\n\n\n' \
                    '$end\n'
    new_string = cleanup_string(string=sample_string,
                                cleanup_types=['BLANK_LINES'])
    assert new_string == '1\n2\n3\n$dummy abcd\n$key2 var2\n    data\n$end\n'

    # Testing : Removing leading spaces preceding $ signs
    sample_string = '1\n2\n3\n' \
                    '   \t$dummy abcd\n  \n' \
                    '$key2 var2\n' \
                    '    data\n\n\n\n\n' \
                    '$end\n'
    new_string = cleanup_string(string=sample_string,
                                cleanup_types=['LEADING_SPACES_DOLLAR'])
    assert new_string == '1\n2\n3\n$dummy abcd\n  \n' \
                         '$key2 var2\n    data\n\n\n\n\n$end\n'

    # Testing : Removing leading spaces preceding # signs
    sample_string = '1\n2\n3\n' \
                    '   \t #this is a comment\n  \n' \
                    '$key2 var2\n' \
                    '    data\n\n\n\n\n' \
                    '  # this is a comment\n\n\n\n\n' \
                    '$end\n'
    new_string = cleanup_string(string=sample_string,
                                cleanup_types=['LEADING_SPACES_HASH'])
    assert new_string == '1\n2\n3\n#this is a comment\n  \n$key2 var2\n' \
                         '    data\n\n\n\n\n' \
                         '# this is a comment\n\n\n\n\n$end\n'

    # Testing : No "$end" in the string
    sample_string = '$title\n' \
                    '$coord  file=coord\n'
    with pytest.raises(expected_exception=ValueError,
                       match=r'^No "\$end" in the string.$'):
        cleanup_string(string=sample_string)

    # Testing: $ in comment before first data group
    # this tests that cleanup_string does not raise an exception even if
    # there is a $ with some text in front of it, because it is in a comment.
    sample_string = '# a comment with a $ inside\n'\
                    '$now the first datagroup\n'\
                    '$end'
    new_string = cleanup_string(string=sample_string,
                                cleanup_types=None)
    assert new_string == sample_string


def test_remove_comments():
    """Testing the removing of comments."""
    # Testing : Removing lines starting with a # sign
    sample_string = 'first line #this inline comment ' \
                    'should not be removed here\n' \
                    '#this line should be removed\n' \
                    '$first dollar sign\n'
    new_string = remove_comments(string=sample_string,
                                 comment_types=['HASH_START'])
    assert new_string == 'first line #this inline comment should not be ' \
                         'removed here\n$first dollar sign\n'

    # Testing : Invalid type of comment
    with pytest.raises(expected_exception=ValueError,
                       match=r'^Comment of type "UNKNOWN" is not a '
                             r'valid type.$'):
        remove_comments(string=sample_string,
                        comment_types=['UNKNOWN'])

    # Testing : Removing $dummy data groups
    sample_string = '1\n2\n3\n' \
                    '$key1\n' \
                    '$dummy abcd\n\n' \
                    '$dummy \n' \
                    '$key2 var2\n' \
                    '    data\n' \
                    '$dummy abcd\n' \
                    'end of file'
    new_string = remove_comments(string=sample_string,
                                 comment_types=['DUMMY'])
    assert new_string == '1\n2\n3\n$key1\n$key2 var2\n    data\n'

    # Testing : Removing lines after $end
    sample_string = '1\n2\n3\n' \
                    '$dummy abcd\n\n' \
                    '$key2 var2\n' \
                    '    data\n' \
                    '$end\n' \
                    'line to be removed\n' \
                    'line to be removed\n'
    new_string = remove_comments(string=sample_string,
                                 comment_types=['AFTER_END'])
    assert new_string == '1\n2\n3\n$dummy abcd\n\n$key2 var2\n    data\n$end\n'

    # Testing : Always a new line after "$end"
    sample_string = '1\n2\n3\n' \
                    '$dummy abcd\n\n' \
                    '$key2 var2\n' \
                    '    data\n' \
                    '$end'
    new_string = remove_comments(string=sample_string,
                                 comment_types=['AFTER_END'])
    assert new_string == '1\n2\n3\n$dummy abcd\n\n$key2 var2\n    data\n$end\n'

    # Testing : Removing "inline" comments
    sample_string = '1\n2\n3\n' \
                    '#This hash starting line should not be removed\n' \
                    '$dummy abcd #This is an inline comment ($dummy is also ' \
                    'a comment by the way but we don\'t remove it here)\n\n' \
                    '$key2 var2\n' \
                    '    data  #  Another inline # comment\n' \
                    '$end'
    new_string = remove_comments(string=sample_string,
                                 comment_types=['HASH_INLINE'])
    assert new_string == '1\n2\n3\n' \
                         '#This hash starting line should not be removed\n' \
                         '$dummy abcd \n\n' \
                         '$key2 var2\n' \
                         '    data  \n' \
                         '$end'


def test_split_string_to_dg_list():
    """Testing splitting of a string to a data group list."""
    sample_string = '$title\n' \
                    '$operating system unix\n' \
                    '$coord    file=coord\n' \
                    '$scfdamp   start=0.600  step=0.050  min=0.100\n' \
                    '$scfintunit\n' \
                    ' unit=30       size=0        file=twoint\n' \
                    '$end\n'
    dg_list = split_string_to_dg_list(sample_string)
    assert len(dg_list) == 6
    assert dg_list == ['$title\n',
                       '$operating system unix\n',
                       '$coord    file=coord\n',
                       '$scfdamp   start=0.600  step=0.050  min=0.100\n',
                       '$scfintunit\n'
                       ' unit=30       size=0        file=twoint\n',
                       '$end\n'
                       ]

    sample_string = '\n\nbefore first dollar\n' \
                    '$title\n' \
                    '$operating system unix\n' \
                    '#This is a comment\n' \
                    '$symmetry c1\n' \
                    '$redundant    file=coord\n' \
                    '$user-defined bonds    file=coord\n' \
                    '$coord    file=coord\n' \
                    '\n\n' \
                    '$dummy\n' \
                    ' this is a comment\n' \
                    '$optimize\n' \
                    ' internal   on\n' \
                    ' redundant  on\n' \
                    '#This is a comment\n' \
                    ' cartesian  off\n' \
                    ' global     off\n' \
                    ' basis      off\n' \
                    '$atoms\n' \
                    'c  1                                                                           \\\n' \
                    '   basis =c SVP\n' \
                    'h  2-5                                                                         \\\n' \
                    '   basis =h SVP\n' \
                    '$basis    file=basis\n' \
                    '$end\n'
    dg_list = split_string_to_dg_list(sample_string)
    assert len(dg_list) == 10
    assert dg_list == ['$title\n',
                       '$operating system unix\n',
                       '$symmetry c1\n',
                       '$redundant    file=coord\n',
                       '$user-defined bonds    file=coord\n',
                       '$coord    file=coord\n',

                       '$optimize\n'
                       ' internal   on\n'
                       ' redundant  on\n'
                       ' cartesian  off\n'
                       ' global     off\n'
                       ' basis      off\n',
                       '$atoms\n'
                       'c  1                                                                           \\\n'
                       '   basis =c SVP\n'
                       'h  2-5                                                                         \\\n'
                       '   basis =h SVP\n',
                       '$basis    file=basis\n',
                       '$end\n']

    sample_string = '$title\n' \
                    '#This is a comment\n' \
                    '$operating system unix\n' \
                    '$coord    file=coord\n' \
                    '\n\n' \
                    '$dummy\n' \
                    ' this is a comment\n' \
                    '$end\n'
    dg_list = split_string_to_dg_list(sample_string,
                                      cleanup_types=[],
                                      remove_comment_types=[])
    assert len(dg_list) == 5
    assert dg_list == ['$title\n'
                       '#This is a comment\n',
                       '$operating system unix\n',
                       '$coord    file=coord\n\n\n',
                       '$dummy\n'
                       ' this is a comment\n',
                       '$end\n'
                       ]


def test_remove_dg_from_list():
    """Testing the removal of one data group from a list."""
    # Testing strict removal
    sample_string = '$title\n' \
                    '$operating\tsystem unix\n' \
                    '$coord    file=coord\n' \
                    '$scfdamp   start=0.600  step=0.050  min=0.100\n' \
                    '$scfintunit\n' \
                    ' unit=30       size=0        file=twoint\n' \
                    '$end\n'
    dg_list = split_string_to_dg_list(sample_string)
    new_dg_list = remove_dg_from_list(dg_to_remove='operating',
                                      dg_list=dg_list,
                                      strict=True)
    assert '$operating\tsystem unix\n' not in new_dg_list

    # Testing not strict removal
    new_dg_list = remove_dg_from_list(dg_to_remove='$opera',
                                      dg_list=dg_list,
                                      strict=False)
    assert '$operating\tsystem unix\n' not in new_dg_list

    new_dg_list = remove_dg_from_list(dg_to_remove='scf',
                                      dg_list=dg_list,
                                      strict=False)
    assert '$operating\tsystem unix\n' in new_dg_list


def test_compare_datagroup_string():
    dg1 = "$forceupdate\n   ahlrichs numgeo=0  mingeo=3 maxgeo=4 modus=<g|dq> dynamic fail=0.3\n   threig=0.005  reseig=0.005  thrbig=3.0  scale=1.00  damping=0.0\n"
    dg1_b = "$forceupdate\n   ahlrichs numgeo=0  mingeo=3 maxgeo=4 modus=<g|dq> dynamic fail=0.3\n   threig=0.005  reseig=0.004  thrbig=3.0  scale=1.00  damping=0.0\n"
    dg1_c = "$forceupdate\n   xxxahlrichs numgeo=0  mingeo=3 maxgeo=4 modus=<g|dq> dynamic fail=0.3\n   threig=0.005  reseig=0.005  thrbig=3.0  scale=1.00  damping=0.0\n"

    assert compare_datagroup_string(dg1, dg1)
    assert not compare_datagroup_string(dg1, dg1_b)
    assert not compare_datagroup_string(dg1, dg1_b, tol= 1e-4)
    assert compare_datagroup_string(dg1, dg1_b, tol= 1e-2)
    assert not compare_datagroup_string(dg1, dg1_c)
    assert not compare_datagroup_string(dg1, dg1_c, tol= 1e-2)

    dg2 = "$thize     0.10000000E-04"
    dg2_b = "$thize 0.10000000E-04"
    dg2_c = "$thize     0.10000000E-05"

    assert compare_datagroup_string(dg2, dg2)
    assert compare_datagroup_string(dg2, dg2_b)
    assert compare_datagroup_string(dg2, dg2_c, tol=1e-4)
    assert not compare_datagroup_string(dg2, dg2_c, tol=1e-6)

    dg3 = "$last step     dscf\n"
    dg3_b = "$last step   \n dscf\n"
    dg3_c = "$last step     dscf a\n"
    assert not compare_datagroup_string(dg3, dg3_b)
    assert not compare_datagroup_string(dg3, dg3_c)


def test_datagroups():
    """Testing the DataGroups object."""
    # Testing : no input
    with pytest.raises(expected_exception=ValueError,
                       match=r'^Both "string" and "dg_list" are None.$'):
        DataGroups()

    # Testing initialization
    dg = DataGroups(string=' $title\n$end')
    assert dg.initial_string == ' $title\n$end'
    assert str(dg) == '$title\n$end\n'
    assert dg.dg_list == ['$title\n', '$end\n']
    assert dg.ndg == 1

    dg = DataGroups(dg_list=['$title title\n',
                             '$coord  file=coord\n',
                             '$end\n'])
    assert dg.initial_string == '$title title\n$coord  file=coord\n$end\n'

    assert str(dg) == '$title title\n$coord  file=coord\n$end\n'

    dg = DataGroups(string=' $title\n$end',
                    dg_list=['$title title\n',
                             '$coord  file=coord\n',
                             '$end\n'])
    assert dg.initial_string == ' $title\n$end'
    assert str(dg) == '$title title\n$coord  file=coord\n$end\n'

    # Testing add data group
    dg.add_data_group('scfiterlimit', '30')
    assert dg.dg_list[-2] == '$scfiterlimit 30\n'
    assert dg.initial_string == ' $title\n$end'
    assert str(dg) == '$title title\n$coord  file=coord\n' \
                      '$scfiterlimit 30\n$end\n'

    # Testing show data group
    coord_data_block = dg.show_data_group('coord')
    assert coord_data_block == '  file=coord\n'
    assert dg.show_data_group("non_existing_dg") is None
    assert dg.show_data_group("non_existing_dg", default="fake_value") == "fake_value"

    # Testing kill data group
    dg.kill_data_group('coord')
    assert dg.dg_list == ['$title title\n',
                          '$scfiterlimit 30\n',
                          '$end\n']

    # Testing add data group with "-" and space
    dg = DataGroups.empty()
    dg.add_data_group('user-defined bonds', '    file=coord')
    assert dg.dg_list[-2] == '$user-defined bonds    file=coord\n'

    # Testing add data group with existing data group
    dg.add_data_group('scfiterlimit', '30')
    with pytest.raises(expected_exception=RuntimeError,
                       match=r'^Data group "\$scfiterlimit" already exists '
                             r'in this DataGroups object.$'):
        dg.add_data_group('scfiterlimit', '60')

    # Testing add data group with invalid name
    dg = DataGroups.empty()
    invalid_name_regex = r'^Data group should start with a letter and be ' \
                         r'followed by alphanumeric characters, a space, "-", "_", "\(" or "\)".$'
    with pytest.raises(expected_exception=ValueError,
                       match=invalid_name_regex):
        dg.add_data_group('0datagroup', '')

    with pytest.raises(expected_exception=ValueError,
                       match=invalid_name_regex):
        dg.add_data_group('data.group', '')

    # Testing show data group with strict=False
    dg = DataGroups(dg_list=['$title title\n',
                             '$coord  file=coord\n',
                             '$scfiterlimit 30\n',
                             '$scfconv 6\n',
                             '$scfdamp   start=0.700  '
                             'step=0.050  min=0.050\n',
                             '$end\n'])
    scfiter_data_block = dg.show_data_group('scfiter', strict=False)
    assert scfiter_data_block == '30\n'
    scfda_data_block = dg.show_data_group('scfda', strict=False)
    assert scfda_data_block == '  start=0.700  step=0.050  min=0.050\n'

    # Testing show data group with strict=False with multiple matches
    with pytest.raises(expected_exception=RuntimeError,
                       match=r'^Found multiple occurrences of data group '
                             r'"\$scf".$'):
        dg.show_data_group('scf', strict=False)

    # Testing change data group
    scfiterlimit_data_block = dg.show_data_group('scfiterlimit', strict=True)
    assert scfiterlimit_data_block == ' 30\n'
    dg.change_data_group('scfiterlimit', '40')
    scfiterlimit_data_block = dg.show_data_group('scfiterlimit', strict=True)
    assert scfiterlimit_data_block == ' 40\n'
    assert str(dg).endswith('step=0.050  min=0.050\n$scfiterlimit 40\n$end\n')
    # when the value is None change_data_group should remove the data group
    dg.change_data_group('scfiterlimit', None)
    assert "scfiterlimit" not in str(dg)

    # Testing as_dict/from_dict
    dg = DataGroups(dg_list=['$title title\n',
                             '$coord  file=coord\n',
                             '$scfiterlimit 30\n',
                             '$scfconv 6\n',
                             '$scfdamp   start=0.700  '
                             'step=0.050  min=0.050\n',
                             '$end\n'])
    dg_dict = dg.as_dict()
    dg_from_dict = DataGroups.from_dict(dg_dict)
    assert dg.dg_list == dg_from_dict.dg_list
    assert str(dg) == str(dg_from_dict)
    assert dg.initial_string == dg_from_dict.initial_string

    dg = DataGroups(dg_list=['$title title\n',
                             '$coord  file=coord\n',
                             '$scfiterlimit 30\n',
                             '$scfconv 6\n',
                             '$scfdamp   start=0.700  '
                             'step=0.050  min=0.050\n',
                             '$end\n'],
                    string='$title\n$end\n')
    dg_dict = dg.as_dict()
    dg_from_dict = DataGroups.from_dict(dg_dict)
    assert dg.dg_list == dg_from_dict.dg_list
    assert str(dg) == str(dg_from_dict)
    assert dg.initial_string == dg_from_dict.initial_string

    # Test to_file and from_file
    with tempfile.TemporaryDirectory() as tmpdir:
        fname = os.path.join(tmpdir, 'control_test')
        dg.to_file(filename=fname)
        dg_from_file = DataGroups.from_file(filename=fname)
        assert dg.dg_list == dg_from_file.dg_list


@pytest.mark.parametrize('control_filename', ['control_test-Control'])
def test_empty(control):
    c = DataGroups.empty()
    assert "$end" == str(c).strip()


def test_modify_data_group_options(testdir):
    cc = DataGroups.from_file(os.path.join(testdir, 'control', "control_test-energy"))
    cc.modify_data_group_options("dft", {"gridsize": "gridsize m5", "ntheta": "ntheta 30"})

    dg = cc.show_data_group("dft")
    mod_lines = set(s.strip() for s in dg.splitlines())
    assert mod_lines == {"gridsize m5", "ntheta 30", "functional b-p", ""}

    cc = DataGroups.from_file(os.path.join(testdir, 'control', "control_test-energy"))
    cc.modify_data_group_options("dft", {"gridsize": None, "ntheta": None})

    dg = cc.show_data_group("dft")
    mod_lines = set(s.strip() for s in dg.splitlines())
    assert mod_lines == {"functional b-p", ""}

    cc = DataGroups.empty()
    cc.modify_data_group_options("dft", {"gridsize": None, "ntheta": "ntheta 30"})
    dg = cc.show_data_group("dft")
    mod_lines = set(s.strip() for s in dg.splitlines())
    assert mod_lines == {"ntheta 30", ""}


def test_show_data_group_option(testdir):
    cc = DataGroups.from_file(os.path.join(testdir, 'control', "control_test-energy"))

    assert cc.show_data_group_option("dft", "gridsize") == "   m3"
    assert cc.show_data_group_option("fake_dg", "fake_option") is None
    assert cc.show_data_group_option("dft", "fake_option") is None
    assert cc.show_data_group_option("dft", "fake_option", default="test") == "test"


def test_sdg_subfiles(testdir, delete_tmp_dir):
    with temp_dir(delete_tmp_dir):
        for fname in ['control', 'coord', 'energy', 'mos', 'basis']:
            shutil.copy2(os.path.join(testdir, 'control',
                                      '{}_test-subfiles'.format(fname)),
                         fname)
        cc = DataGroups.from_file("control")
        energy_dg = cc.sdg(data_group='energy', show_from_subfile=True,
                           raise_if_multiple_subfiles=False,
                           raise_if_missing_subfile=False,
                           raise_if_regular_and_subfile=False)
        assert energy_dg == '      SCF               ' \
                            'SCFKIN            SCFPOT\n' \
                            '     1   -40.48328750375' \
                            '    39.86406749278   -80.34735499653\n'
        udef_bonds_dg = cc.sdg(data_group='user-defined bonds',
                               show_from_subfile=True,
                               raise_if_multiple_subfiles=False,
                               raise_if_missing_subfile=False,
                               raise_if_regular_and_subfile=False)
        assert udef_bonds_dg == '\n'
        udef_bonds_dg = cc.sdg(data_group='user-defined bonds',
                               show_from_subfile=False,
                               raise_if_multiple_subfiles=False,
                               raise_if_missing_subfile=False,
                               raise_if_regular_and_subfile=False)
        assert udef_bonds_dg == '    file=coord\n'
        redundant_dg = cc.sdg(data_group='redundant',
                              show_from_subfile=True,
                              raise_if_multiple_subfiles=False,
                              raise_if_missing_subfile=False,
                              raise_if_regular_and_subfile=False)
        assert redundant_dg.startswith('\n     number_of_atoms'
                                       '             5\n'
                                       '     degrees_of_freedom'
                                       '          9\n'
                                       '     internal_coordinates'
                                       '        9\n'
                                       '     frozen_coordinates'
                                       '          0\n'
                                       '   1 k  1.0000000000000 stre'
                                       '    1    2'
                                       '           val=   2.06326\n')
        assert redundant_dg.endswith('         7\n           8'
                                     '           0.699953222    8    0\n'
                                     '         8\n           9'
                                     '           0.699953222    9    0\n'
                                     '         9\n')
        basis_dg = cc.sdg(data_group='basis',
                          show_from_subfile=True,
                          raise_if_multiple_subfiles=False,
                          raise_if_missing_subfile=False,
                          raise_if_regular_and_subfile=False)
        assert basis_dg.startswith('\n*\nc SVP\n*\n')
        assert basis_dg.endswith('   1  p\n'
                                 ' 0.80000000000       1.0000000000\n*\n')
        mos_dg = cc.sdg(data_group='scfmo',
                        show_from_subfile=True,
                        raise_if_multiple_subfiles=False,
                        raise_if_missing_subfile=False,
                        raise_if_regular_and_subfile=False)
        assert mos_dg.startswith('    scfconv=6   format(4d20.14)\n'
                                 '     1  a1     '
                                 'eigenvalue=-.98973943603420D+01'
                                 '   nsaos=6\n')
        assert mos_dg.endswith('-.11000018463806D+01-.60830569311504D+00'
                               '0.13711445998604D+010.13421520432673D+01\n'
                               '0.53804368223874D+00-.91129252301913D+00'
                               '0.11101532709297D+01\n')
        scfintunit_dg = cc.sdg(data_group='scfintunit',
                               show_from_subfile=True,
                               raise_if_multiple_subfiles=False,
                               raise_if_missing_subfile=False,
                               raise_if_regular_and_subfile=False)
        assert scfintunit_dg == '\n unit=30       size=35' \
                                '       file=twoint1\n' \
                                ' unit=31       size=35' \
                                '       file=twoint2\n'
        with pytest.raises(RuntimeError,
                           match='Multiple "file=FILENAME" directives.'):
            cc.sdg(data_group='scfintunit',
                   show_from_subfile=True,
                   raise_if_multiple_subfiles=True,
                   raise_if_missing_subfile=False,
                   raise_if_regular_and_subfile=False)

        with pytest.raises(RuntimeError,
                           match=r'Both a reference to a file '
                                 r'\("file=FILENAME"\) and regular data '
                                 r'blocks are in the data group.'):
            cc.sdg(data_group='user-defined',
                   show_from_subfile=True,
                   raise_if_multiple_subfiles=False,
                   raise_if_missing_subfile=False,
                   raise_if_regular_and_subfile=True)

        with pytest.raises(RuntimeError,
                           match='File "gradient" for data group "grad" '
                                 'is missing.'):
            cc.sdg(data_group='grad',
                   show_from_subfile=True,
                   raise_if_multiple_subfiles=False,
                   raise_if_missing_subfile=True,
                   raise_if_regular_and_subfile=False)

        grad_db = cc.sdg(data_group='grad',
                         show_from_subfile=True,
                         raise_if_multiple_subfiles=False,
                         raise_if_missing_subfile=False,
                         raise_if_regular_and_subfile=False)
        assert grad_db == '    file=gradient\n'

        # create a gradient empty file and check that the behavior is the same
        # as when the file is missing
        touch_file("gradient")
        with pytest.raises(RuntimeError,
                           match='File "gradient" for data group "grad" '
                                 'is missing.'):
            cc.sdg(data_group='grad',
                   show_from_subfile=True,
                   raise_if_multiple_subfiles=False,
                   raise_if_missing_subfile=True,
                   raise_if_regular_and_subfile=False)

        grad_db = cc.sdg(data_group='grad',
                         show_from_subfile=True,
                         raise_if_multiple_subfiles=False,
                         raise_if_missing_subfile=False,
                         raise_if_regular_and_subfile=False)
        assert grad_db == '    file=gradient\n'


def test_show_subfile_fname(testdir):
    cc = DataGroups.from_file(os.path.join(testdir, 'control', "control_test-energy"))
    assert cc.show_subfile_fname("redundant") == "coord"

    with pytest.raises(RuntimeError):
        cc.show_subfile_fname("scfintunit", raise_if_regular_and_subfile=True)


def test_compare(testdir):
    cc1 = DataGroups.from_file(os.path.join(testdir, 'control', "control_test-energy"))
    cc2 = DataGroups.from_file(os.path.join(testdir, 'control', "control_test-energy"))
    assert cc1.compare(cc2) is None

    dg_missing = cc2.dg_list.pop(0)
    assert dg_missing in cc1.compare(cc2)

    cc2 = DataGroups.from_file(os.path.join(testdir, 'control', "control_test-energy"))
    cc2.cdg("symmetry", "c1")
    assert "symmetry" in cc1.compare(cc2)

    cc2 = DataGroups.from_file(os.path.join(testdir, 'control', "control_test-energy"))
    cc2.cdg("thize", "0.001")
    assert "thize" in cc1.compare(cc2)
    assert "thize" in cc1.compare(cc2, tol=1e-6)
    assert cc1.compare(cc2, tol=1e-2) is None
