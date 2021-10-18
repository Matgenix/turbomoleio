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

from unittest import mock
import turbomoleio.core.utils
from monty.tempfile import ScratchDir
import shutil
import os


class TestGetVersion:

    def test_V72(self, testdir):
        def copy_define_log():
            shutil.copy(
                os.path.join(testdir, "define", "define-qq.log_V7.2"),
                "define.log",
            )

        with mock.patch("turbomoleio.core.utils.define_quit", side_effect=copy_define_log):
            with ScratchDir('.'):
                tm_version = turbomoleio.core.utils.get_tm_version()
                assert tm_version == '7.2'

    def test_V731(self, testdir):
        def copy_define_log():
            shutil.copy(
                os.path.join(testdir, "define", "define-qq.log_V7.3.1"),
                "define.log",
            )

        with mock.patch("turbomoleio.core.utils.define_quit", side_effect=copy_define_log):
            with ScratchDir('.'):
                tm_version = turbomoleio.core.utils.get_tm_version()
                assert tm_version == '7.3.1'

    def test_V741(self, testdir):
        def copy_define_log():
            shutil.copy(
                os.path.join(testdir, "define", "define-qq.log_V7.4.1"),
                "define.log",
            )

        with mock.patch("turbomoleio.core.utils.define_quit", side_effect=copy_define_log):
            with ScratchDir('.'):
                tm_version = turbomoleio.core.utils.get_tm_version()
                assert tm_version == '7.4.1'

    def test_V751(self, testdir):
        def copy_define_log():
            shutil.copy(
                os.path.join(testdir, "define", "define-qq.log_V7.5.1"),
                "define.log",
            )

        with mock.patch("turbomoleio.core.utils.define_quit", side_effect=copy_define_log):
            with ScratchDir('.'):
                tm_version = turbomoleio.core.utils.get_tm_version()
                assert tm_version == '7.5.1'
