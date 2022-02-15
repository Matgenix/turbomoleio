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

"""Utility module for the "data group" model of TurboMole.

This module contains utility functions to parse, manipulate and write files
formatted using the "data group" convention of TurboMole. Examples of files
that are written in this fashion are the "control" file, the "coord" file, the
"basis" file and many other files written by TurboMole.

Files using the "data group" model of TurboMole are human readable ASCII files
in which data is formally split into consistent "groups". Each "data group"
consists of a keyword starting with a dollar ("$") sign (e.g. "$symmetry",
"$open shell", "$optimize" ...) and its associated data. The associated data
can be empty or just one value/word but it can also be very large (e.g.
molecular orbitals : "$scfmo" data group contained in the "mos" file).

A complete description of the data group model of TurboMole can be found in
the user's manual, as well as a description of all the keywords available.
"""

import re
import os
from collections import namedtuple

from monty.json import MSONable


def cleanup_string(string, cleanup_types=None):
    """Returns a cleaned data groups string.

    This function allows to clean up a data groups string. It allows to :
    - remove everything before the first dollar ("$") sign,
    - remove blank lines,
    - remove spaces and tabs before "#" comments,
    - remove spaces and tabs before "$" signs.

    Args:
        string (str): A string that should be in the "data group" format of
            TurboMole. Typically, this string is read from file (e.g. control,
            coord, basis, ...).
        cleanup_types (list of str, optional): List of the cleanup types that
            should be performed. Valid types are "BEFORE_FIRST_DOLLAR",
            "BLANK_LINES", "LEADING_SPACES_DOLLAR" and "LEADING_SPACES_HASH".
            Defaults to None in which case all four types of cleanup will be
            performed.

    Returns:
        str: A cleaned up string.

    Raises:
        ValueError: If one of the cleanup_types is not valid, if there is no
            dollar ("$") sign, if one of the dollar ("$") or hash ("#") signs
            is preceded by characters other than a space or a tab in a line or
            if there is no "$end" at the end of the file.
    """
    if '$' not in string:
        raise ValueError('No dollar ("$") sign in the string.')
    if cleanup_types is None:
        cleanup_types = ["BLANK_LINES", "LEADING_SPACES_DOLLAR",
                         "LEADING_SPACES_HASH", "BEFORE_FIRST_DOLLAR"]
    for cleanup_type in cleanup_types:
        # Removing lines before the first dollar sign (except if this dollar
        # sign is in a hash comment)
        if cleanup_type == 'BEFORE_FIRST_DOLLAR':
            newlines = []
            after_first_dollar = False
            for line in string.split('\n'):
                if after_first_dollar or line.strip().startswith('#'):
                    newlines.append(line)
                elif line.strip().startswith('$'):
                    after_first_dollar = True
                    newlines.append(line.lstrip())
            string = '\n'.join(newlines)
        # Removing all blank lines
        elif cleanup_type == 'BLANK_LINES':
            string = re.sub(r'(\n[ \t]*){2,}', '\n', string,
                            flags=re.MULTILINE)
        # Removing all spaces and tabs before dollar signs
        elif cleanup_type == 'LEADING_SPACES_DOLLAR':
            string = re.sub(r'^[ \t]*\$', '$', string,
                            flags=re.MULTILINE)
        # Removing all spaces before hash signs
        elif cleanup_type == 'LEADING_SPACES_HASH':
            string = re.sub(r'^[ \t]*#', '#', string,
                            flags=re.MULTILINE)
        else:
            raise ValueError('Cleanup of type "{}" '
                             'is not a valid type.'.format(cleanup_type))

    for match in re.findall(pattern=r'^[ \t\S]+\$',
                            string=string,
                            flags=re.MULTILINE):
        if not match.strip().startswith('#'):
            raise ValueError('Some character(s) are preceding a '
                             'dollar ("$") sign.')

    if '$end' not in string:
        raise ValueError('No "$end" in the string.')
    return string


def remove_comments(string, comment_types=None):
    """Returns a data group string with comment lines removed.

    Comments or additional information in TurboMole is typically specified in
    three different ways : a line starting with a hash ("#") symbol, a data
    group with the keyword "$dummy" (can span multiple lines), anything
    after the "$end" keyword or anything after a hash ("#") symbol (so-called
    "inline" comments). This function removes all types of comments in
    a data group formatted string. Note also that if a newline is added after
    the "$end" keyword if not present.

    Args:
        string (str): A string that should be in the "data group" format of
            TurboMole. Typically, this string is read from file (e.g. control,
            coord, basis, ...).
        comment_types (:obj:`list`, optional): List of the types of comments
            that should be removed. Valid types are "HASH_START", "DUMMY",
            "AFTER_END" and "HASH_INLINE". Defaults to None in which case all
            four types of comments will be removed.

    Returns:
        str: A string without the comment lines.

    Raises:
        ValueError: If one of the comment_types is not valid.
    """
    if comment_types is None:
        comment_types = ["HASH_START", "DUMMY", "AFTER_END", "HASH_INLINE"]
    for comment_type in comment_types:
        # Removing lines starting with a hash symbol
        if comment_type == 'HASH_START':
            string = re.sub(r'^#.*\n?', '', string, flags=re.MULTILINE)
        # Removing portions of the file corresponding to a $dummy data group
        elif comment_type == 'DUMMY':
            old = None
            # Loop is needed because replaced characters are not used in
            # consecutive regex searches
            while old != string:
                old = string
                string = re.sub(r'^\$dummy[^$]*\$', '$', string,
                                flags=re.MULTILINE)
            # Dummy data group at the end of the file
            string = re.sub(r'^\$dummy[\s\S]*\Z', '', string,
                            flags=re.MULTILINE)
        # Removing lines after $end
        elif comment_type == 'AFTER_END':
            string = re.sub(r'^\$end[\s\S]*\Z', '$end\n', string,
                            flags=re.MULTILINE)
        # Removing inline comments i.e. #-starting comments in the middle of
        # a line (e.g. : "$symmetry c1 # Symmetry of the system")
        # Note that lines starting with a hash symbol will also be replaced
        # by a blank line
        elif comment_type == 'HASH_INLINE':
            # Only matches if there is something other than white spaces
            # before the hash ("#") symbol. Only the first capturing group
            # is reinserted with the "\1", while the lookahead ("(?=$)") is
            # checking for the end of the line but leaving it as it is, be it
            # a newline or an end of file.
            string = re.sub(r'(^\s*\S.*?)(#.*)(?=$)', r'\1', string, flags=re.MULTILINE)
        else:
            raise ValueError('Comment of type "{}" '
                             'is not a valid type.'.format(comment_type))
    return string


def split_string_to_dg_list(string, cleanup_types=None,
                            remove_comment_types=None):
    r"""Splits a string into a list of data groups.

    Args:
        string (str): A string that should be in the "data group" format of
            TurboMole. Typically, this string is read from file (e.g. control,
            coord, basis, ...).
        cleanup_types (list of str, optional): List of the cleanup types that
            should be performed before splitting the string into a list of
            data groups. Valid types are "BEFORE_FIRST_DOLLAR",
            "HASH_START", "DUMMY", "AFTER_END" and "BLANK_LINES". Defaults to
            None in which case all five types of comments will be removed.
        remove_comment_types (list of str, optional): List of the types of
            comments that should be removed before splitting the string into a
            list of data groups. Valid types are "BEFORE_FIRST_DOLLAR",
            "HASH_START", "DUMMY", "AFTER_END" and "BLANK_LINES". Defaults to
            None in which case all five types of comments will be removed.

    Returns:
        list: The list of data groups in the string.

    Raises:
        ValueError: If there is no dollar ("$") sign in the string, i.e.
            the string is not or incorrectly formatted with respect to the
            data group format of TurboMole.

    Examples:
        >>> split_string_to_dg_list("$title\n$end\n")
        ['$title\n', '$end\n']
        >>> split_string_to_dg_list("no dollar sign")
        []
        >>> split_string_to_dg_list("\n\n$title\n"
        ...                         "$optimize\n"
        ...                         " internal   on\n"
        ...                         " redundant  on\n"
        ...                         " cartesian  off\n"
        ...                         " global     off\n"
        ...                         " basis      off\n")
        ['$title\n', '$optimize\n internal   on\n redundant  on\n cartesian  off\n global     off\n basis      off\n']
    """
    string = cleanup_string(string=string, cleanup_types=cleanup_types)
    string = remove_comments(string=string, comment_types=remove_comment_types)
    return ['$' + dg for dg in string.split('$')[1:]]


def remove_dg_from_list(dg_to_remove, dg_list, strict=True):
    """Returns new list with a data group removed from the initial list.

    This function will remove dg_to_remove from a list of data groups. If the
    data group is duplicated, all occurrences are removed. The standard usage
    is to remove the data group if it is an exact match (i.e. if the data
    group to be removed is followed by a space, a tab, a new line or a return).
    One can also remove all data groups starting with `dg_to_remove` by using
    strict=False.

    Args:
        dg_to_remove (str): Data group to be removed from the list. Can be
            given with or without the dollar ("$") sign.
        dg_list (list): List of strings. Each string should be a data group,
            i.e. should start with a dollar ("$") sign and should not contain
            any other dollar sign.
        strict (bool): If False, all the data groups starting with
            `dg_to_remove` will be removed from the list. Otherwise, the data
            group is removed only if it is an exact match of `dg_to_remove`.

    Returns:
        list: The list of data groups with `dg_to_remove` removed.
    """
    if dg_to_remove[0] == "$":
        dg_to_remove = dg_to_remove[1:]
    new_dg_list = []
    if strict:
        pattern = r'^\$' + dg_to_remove + r'\s'
    else:
        pattern = r'^\$' + dg_to_remove
    for dg in dg_list:
        if not re.match(pattern=pattern, string=dg):
            new_dg_list.append(dg)
    return new_dg_list


def compare_datagroup_string(dg1, dg2, tol=None):
    """
    Compares to datagroup strings as split from the DataGroups object.
    Return True if the two match.
    The two strings should have the same number of lines. Each line should have the
    same number of chunks. Chunks will be separated using spaces and "=". Each chunk
    should match exactly as a string, except if they are numbers. In that, if a tol
    is specified they will be converted to float and their difference should be lower
    than the tolerance.

    Args:
        dg1 (str): the first datagroup string to be compared.
        dg2 (str): the second datagroup string to be compared.
        tol (float): the tolerance allowed when comparing numbers. If None
            even the number should match as strings.

    Returns:
        bool: True if the string match.
    """

    lines1 = dg1.splitlines()
    lines2 = dg2.splitlines()

    if len(lines1) != len(lines2):
        return False

    for l1, l2 in zip(lines1, lines2):
        strings1 = re.split(r"\s+|=", l1)
        strings2 = re.split(r"\s+|=", l2)
        if len(strings1) != len(strings2):
            return False
        for s1, s2 in zip(strings1, strings2):
            if tol is not None:
                try:
                    f1 = float(s1.replace("D", "E"))
                    f2 = float(s2.replace("D", "E"))
                    if abs(f1 - f2) > tol:
                        return False
                except ValueError:
                    if s1 != s2:
                        return False
            elif s1 != s2:
                return False

    return True


class DataGroups(MSONable):
    """Generic class for data group formatted files and strings.

    DataGroups is a generic class for parsing, manipulating and generating
    strings/files using the data group format of TurboMole.
    """

    def __init__(self, string=None, dg_list=None):
        """Initializes a `DataGroups` class.

        The DataGroups object defines a set of key-value-like pairs. A data
        group is a key starting with a dollar ("$") sign. The value
        corresponding to the data group is called the data block. When
        initializing `DataGroups` from a string (e.g. reading a file), comments
        and/or unnecessary lines/spaces/... may be present in that string. In
        that case, the initial string is kept as a reference. Note that if
        changes are applied to the DataGroups object, the initial string is
        not updated as comments may interfere with the modifications. The
        actual string corresponding to the data group list (`dg_list`) is
        always in line with the data group list.

        Args:
            string (str): A string in the data group format.
            dg_list (list of str): A list of the data groups.
        """
        if dg_list is None:
            if string is None:
                raise ValueError('Both "string" and "dg_list" are None.')
            self.initial_string = string
            self.dg_list = split_string_to_dg_list(string=string,
                                                   cleanup_types=None,
                                                   remove_comment_types=None)
        else:
            self.dg_list = dg_list
            if string is None:
                self.initial_string = ''.join(dg_list)
            else:
                self.initial_string = string

    def kill_data_group(self, data_group, strict=True):
        """Removes `data_group` from this `DataGroups` object.

        Args:
            data_group (str): Data group to be removed.
            strict (bool): If True `data_group` should be an exact match.
                If False, any Data group starting with `data_group` will be
                removed.
        """
        self.dg_list = remove_dg_from_list(dg_to_remove=data_group,
                                           dg_list=self.dg_list,
                                           strict=strict)

    kdg = kill_data_group

    def add_data_group(self, data_group, data_block):
        """Adds `data_group`->`data_block` to this `DataGroups` object.

        Args:
            data_group (str): Data group (key) to be added. The dollar ("$")
                sign will be automatically added if not present.
            data_block (str): Data block corresponding to the data group.

        Raises:
            RuntimeError: if data_group already exists.
        """
        if data_group[0] != '$':
            data_group = '$' + data_group
        if self.show_data_group(data_group=data_group,
                                strict=True) is not None:
            raise RuntimeError('Data group "{}" already exists in this '
                               'DataGroups object.'.format(data_group))
        if not re.fullmatch(pattern=r'^\$[a-zA-Z][ \-\w\(\)]*', string=data_group):
            raise ValueError('Data group should start with a letter and be '
                             'followed by alphanumeric characters, a space, '
                             '"-", "_", "(" or ")".')
        if not data_block.endswith('\n'):
            data_block = data_block + '\n'
        if data_block.startswith((' ', '\t', '\n')):
            self.dg_list.insert(-1, '{}{}'.format(data_group, data_block))
        else:
            self.dg_list.insert(-1, '{} {}'.format(data_group, data_block))

    adg = add_data_group

    def show_data_group(self, data_group, strict=True, default=None,
                        show_from_subfile=True,
                        raise_if_multiple_subfiles=False,
                        raise_if_missing_subfile=False,
                        raise_if_regular_and_subfile=False):
        """Shows `data_group` from this `DataGroups` object.

        Args:

            data_group (str): Data group (key) to be added. The dollar ("$")
                sign will be automatically added if not present.
            strict (bool): Whether `data_group` should be an exact match or
                if data groups starting with `data_group` are allowed.
            default (str): the default value that will be returned if the
                data_group is not present in the list. Default is None.
            show_from_subfile (bool): If True, will show `data_group` from
                within the "subfile" if the data block contains a
                "file=FILENAME". If False the "file=FILENAME" block is simply
                returned. This supposes that file exists in the current
                directory.
            raise_if_multiple_subfiles (bool): Whether to raise an error if
                multiple "file=" directives are present in this data group. If
                False, just returns the standard control data block.
            raise_if_missing_subfile (bool): Whether to raise an error if the
                subfile does not exist in the current directory or exists but
                is empty. If False, just returns the standard data block.
            raise_if_regular_and_subfile (bool): Whether to raise an error if
                data group contains both a reference to a file with a
                "file=FILENAME" and regular data block options. If
                False, just returns the standard data block.

        Returns:
            str: the value of the selected datagroup. None if no `data_group`
                is found.

        Raises:
            RuntimeError: if multiple occurrences of `data_group` are found
                (and raise_if_multiple_subfiles is True), if the subfile is
                missing or empty (and raise_if_missing_subfile is True) or if
                regular data block options coexist with a reference to a subfile
                 (and raise_if_regular_and_subfile is True).
        """

        if data_group[0] == "$":
            data_group = data_group[1:]
        if strict:
            pattern = r'^\$' + data_group + r'[=\s]'
        else:
            pattern = r'^\$' + data_group

        matches = []
        for dg in self.dg_list:
            if re.match(pattern=pattern, string=dg):
                matches.append(dg)
        if len(matches) == 0:
            return default
        dollar_data_group = '$' + data_group
        if len(matches) > 1:
            raise RuntimeError('Found multiple occurrences of data group '
                               '"{}".'.format(dollar_data_group))
        if strict:
            data_block = matches[0].replace(dollar_data_group, '')
        else:
            match = re.sub(r'^\$' + data_group + r'[-\w]*', '', matches[0])
            data_block = re.sub(r'^ ', '', match)

        if not show_from_subfile:
            return data_block

        subfiles_count = data_block.count('file=')
        if subfiles_count == 0:
            return data_block
        elif subfiles_count == 1:
            filename = self._get_subfile_fname(data_block=data_block,
                                               raise_if_regular_and_subfile=raise_if_regular_and_subfile)
            if os.path.exists(filename):
                with open(filename, 'r') as f:
                    subfile_string = f.read()

                if subfile_string.strip():
                    subfile = DataGroups(subfile_string)
                    return subfile.show_data_group(data_group=data_group)

            if raise_if_missing_subfile:
                raise RuntimeError('File "{}" for data group "{}" is '
                                   'missing or empty.'.format(filename, data_group))
            return data_block
        elif raise_if_multiple_subfiles:
            raise RuntimeError('Multiple "file=FILENAME" directives.')

        return data_block

    sdg = show_data_group

    def change_data_group(self, data_group, data_block):
        """
        Changes the value of the `data_group` to the new `data_block` in this
        `DataGroups` object. If the key is not present will be added.

        It will first apply the kill_data_group and then add_data_group. Note
        that in the current implementation, the data group will be moved to
        the end of the data groups object (just before $end).

        If `data_block` is None it is equivalent to kill_data_group. To set
        a data group with no explicit value use an empty string for `data_block`
        (e.g. to add $uhf the input should be data_group="uhf" and data_block="").

        Args:
            data_group (str): Data group (key) to be changed. The dollar ("$")
                sign will be automatically added if not present.
            data_block (str): Data block corresponding to the data group.
                If None it will simply kill_data_group for the specified
                data_group value.
        """

        self.kill_data_group(data_group=data_group, strict=True)
        if data_block is not None:
            self.add_data_group(data_group=data_group, data_block=data_block)

    cdg = change_data_group

    def modify_data_group_options(self, data_group, options):
        """
        Given a data group that allows several options on separate line
        (e.g. $dft), updates the values of the options according to the
        dictionary provided.
        The option dictionary should have the form:

        .. code-block:: text

            {"option_name1": "option_name1 option_value",
             "option_name2": "option_name2=option_value"}

        The key will be used to identify the line to be modified and that line
        would be entirely replaced by the value.
        Since different options may be defined in different ways, no
        attempt is made here to identify the suitable format for the option.
        It is responsibility of the caller to specify the line for the option
        in the correct format.

        If the datagroup is not present will be created with the specified options.

        If the entire data group should be modified it would be safer to use
        change_data_group.

        Args:
            data_group (str): Data group (key) to be changed. The dollar ("$")
                sign will be automatically added if not present.
            options (dict): The options that should be added or updated. If the
                value of one key is None the option will be removed if present.
        """

        dg = self.show_data_group(data_group, strict=True)

        # start from a line  with a space if dg is not present (should not be empty)
        if dg is None:
            dg = " "

        # copy the dict as it will modified
        options = dict(options)
        new_lines = []
        for l in dg.splitlines():
            # check if one of the lines matches the options. If yes replace
            # the line and remove it from the optiond dict, otherwise keep
            # the oriiginal line.
            for opt in list(options.keys()):
                if opt in l:
                    opt_kv = options.pop(opt)
                    # if None skip the line
                    if opt_kv is not None:
                        new_lines.append("   " + opt_kv)
                    break
            else:
                new_lines.append(l)

        # add all the options that were not present in the original set
        # except for None. Sorted so that the results are deterministic
        new_values = set(options.values())
        new_values.discard(None)
        for opt_kv in sorted(new_values):
            new_lines.append("   " + opt_kv)

        self.change_data_group(data_group, "\n".join(new_lines))

    mdgo = modify_data_group_options

    def show_data_group_option(self, data_group, option, default=None):
        """
        Given a data group that allows several options on separate line
        (e.g. $dft), returns the value of the option provided.

        Since each value may be defined in a different way, the returned
        value will include anything after the option specified. It is up
        to the caller determine if there are symbols like "=" that should
        be removed, depending on the data group and option that is queried.

        Args:
            data_group (str): Data group (key) to be shown. The dollar ("$")
                sign will be automatically added if not present.
            option (str): The options that should be returned.
            default (str): the default value that will be returned if the
                data_group or the option are not present in the list.

        Returns:
            str: The value of the option.
        """

        dg = self.show_data_group(data_group, strict=True, show_from_subfile=True,
                                  raise_if_missing_subfile=False)

        if dg is None:
            return default

        compiled_re = re.compile(r"\s*{}(.*)$".format(option), re.DOTALL)
        for l in dg.splitlines():
            match = compiled_re.match(l)
            if match:
                return match.group(1)

        return default

    sdgo = show_data_group_option

    @property
    def number_of_data_groups(self):
        return len(self.dg_list) - 1

    ndg = number_of_data_groups

    @classmethod
    def empty(cls):
        """Creates an empty DataGroups object.

        An empty DataGroups object only contains the "$end" data group.
        """
        return cls(dg_list=['$end\n'])

    def as_dict(self):
        """Returns a dictionary representing the DataGroups object."""
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "string": self.initial_string,
                "dg_list": self.dg_list}

    @classmethod
    def from_dict(cls, d):
        """Creates DataGroups object from a dictionary."""
        return cls(string=d['string'], dg_list=d['dg_list'])

    def __str__(self):
        """Returns string representation of this DataGroups object.

        The string representation of the DataGroups object is supposed to be
        written to a file and used as an input to TurboMole executables.
        """
        return ''.join(self.dg_list)

    def to_file(self, filename):
        """Writes this DataGroups object to a file.
        If the file exists it will be overwritten.

        Args:
            filename (str): Name of the file to which this DataGroups object
                should be written.
        """
        with open(filename, 'w') as f:
            f.write(self.__str__())

    @classmethod
    def from_file(cls, filename):
        """Creates DataGroups object reading from a given file.

        Args:
            filename (str): Name of the file from which this DataGroups object
                should be read.
        """
        with open(filename, 'r') as f:
            string = f.read()
        return cls(string=string)

    @staticmethod
    def _get_subfile_fname(data_block, raise_if_regular_and_subfile=False):
        """
        Extract the name of the file from a datagroup string.

        Args:
            data_block (str): the string with the content of the data group.
            raise_if_regular_and_subfile (bool): if True will raise an error if the
                string contains both a file=FILENAME and other data. If False will
                return the filename anyway.

        Returns:
            (str): the name of the file, None if no match for "file=filename" is found.

        Raises:
            RuntimeError: if raise_if_regular_and_subfile is True and file=FILENAME
                is not the only content of the string.
        """
        pattern = r'file=([^\s]+)'
        match = re.search(pattern=pattern, string=data_block)

        if not match:
            return None

        data_block_test = data_block.replace(match.group(0), '')
        data_block_test = data_block_test.strip()
        if data_block_test != '' and raise_if_regular_and_subfile:
            raise RuntimeError('Both a reference to a file '
                               '("file=FILENAME") and regular data '
                               'blocks are in the data group.')
        return match.group(1)

    def show_subfile_fname(self, data_group, raise_if_regular_and_subfile=False):
        """
        Extract the name of the file from a datagroup.

        Args:
            data_group (str): the string with the content of the data group.
            raise_if_regular_and_subfile (bool): if True will raise an error if the
                string contains both a file=FILENAME and other data. If False will
                return the filename anyway.

        Returns:
            (str): the name of the file, None if no match for "file=filename" is found.

        Raises:
            RuntimeError: if raise_if_regular_and_subfile is True and file=FILENAME
                is not the only content of the string.
        """
        return self._get_subfile_fname(self.sdg(data_group, show_from_subfile=False),
                                       raise_if_regular_and_subfile=raise_if_regular_and_subfile)

    def compare(self, datagroups, tol=None, ignored_dg=None, return_all_diffs=False):
        """
        Compares the current object with another DataGroups. Aside from the datagroups
        that should be ignored, (define in the ignored_dg argument), the two should contain
        the same amount of elements in dg_list and each one from the current object
        should have a matching elements in the one that is passed as argument.
        To match two strings should have the same number of lines. Each line should have the
        same number of chunks. Each chunk should match exactly as a string, except if
        they are numbers. In that, if a tol is specified they will be converted to float
        and their difference should be lower than the tolerance.

        If the two objects match None will be returned, otherwise a message describing
        the property that caused the match to fail.

        Note that, given the variety of notations supported by Turbomole, this function is to
        be intended more for testing purposes rather than as a complete tool to guarantee the
        equivalence of two control file.

        Args:
            datagroups (DataGroups): a Datagroups that should match with the current instance.
            tol (float): the tolerance allowed when comparing numbers. If None
                even the number should match as strings.
            ignored_dg (list): a list of datagroups that should be ignored for the comparison.
                It is expected to be the name of the datagroup, thus starting with a "$".
                If the "$" symbol is not present will be added before trying to match the
                datagroup to skip.
            return_all_diffs (bool): If True, a list of all differences is returned. If False
                (default), a string describing the first source of difference found is returned.

        Returns:
            None if the results match or a string describing a source of difference otherwise
                (if return_all_diffs is False) or a list of all differences (if return_all_diffs
                is True).
        """

        dg_other_list = list(datagroups.dg_list)
        diffs = []

        if ignored_dg:
            # add the initial $ to each datagroup if not already present
            ignored_dg = ["$" + dg if not dg.startswith("$") else dg for dg in ignored_dg]

        for dg1 in self.dg_list:
            # skip the check if dg1 matches the list of datagroups that should be ignored
            if ignored_dg and any(dg1.startswith(dg) for dg in ignored_dg):
                continue
            for i, dg2 in enumerate(dg_other_list):
                if compare_datagroup_string(dg1, dg2, tol):
                    dg_other_list.pop(i)
                    break
            else:
                msg = "Datagroup does not match to any of the references: {}".format(dg1)
                if return_all_diffs:
                    diffs.append(msg)
                else:
                    return msg

        # check that all the remaining datagroups in the second control belong to the ignore list
        for dg2 in dg_other_list:
            if not ignored_dg or not any(dg2.startswith(dg) for dg in ignored_dg):
                msg = "Datagroup in the reference does not match to any of the current " \
                      "control: {}".format(dg2)
                if return_all_diffs:
                    diffs.append(msg)
                else:
                    return msg

        if return_all_diffs:
            return diffs if len(diffs) > 0 else None
        return None
