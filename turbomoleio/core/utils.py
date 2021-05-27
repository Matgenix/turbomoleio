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
Module containing miscellaneous utilities
"""

from monty.tempfile import ScratchDir
import re
import pexpect.popen_spawn
import signal


def define_quit():
    """Runs define but directly quits by sending "qq" interactively.
    """
    with open('define.log', 'wb') as logf:
        spawn = pexpect.popen_spawn.PopenSpawn('define', timeout=2, logfile=logf)
        try:
            spawn.expect(["D E F I N E"])
            spawn.sendline('qq')
        finally:
            if spawn.proc.poll() is None:
                spawn.kill(sig=signal.SIGKILL)


def get_tm_version():
    """Get the turbomole version currently in use.

    This basically runs a define and exits at the very beginning and extracts
    the version from the header of the output.
    """
    with ScratchDir('.'):
        define_quit()
        with open('define.log', 'r') as f:
            data = f.read()
    r = r'TURBOMOLE.*V([\d\.]+\d)\s+.*'
    match_version = re.search(r, data)
    return match_version.group(1)
