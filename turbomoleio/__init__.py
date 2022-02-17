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

"""
Turbomoleio is a package acting as an interface with the TURBOMOLE package.

It allows to generate inputs and parse outputs easily.
This is the root package.
"""

# Useful aliases for commonly used objects and modules.
# Allows from turbomoleio import <class> for quick usage.
from turbomoleio.core.control import Control  # noqa
from turbomoleio.core.datagroups import DataGroups  # noqa
from turbomoleio.core.molecule import MoleculeSystem  # noqa
from turbomoleio.input.define import DefineRunner  # noqa
from turbomoleio.input.utils import get_define_template  # noqa
from turbomoleio.output.files import (  # noqa
    AoforceOutput,
    EgradOutput,
    EscfOnlyOutput,
    EscfOutput,
    GradOutput,
    JobexOutput,
    RelaxOutput,
    ScfOutput,
    StatptOutput,
)
from turbomoleio.output.parser import Parser  # noqa
from turbomoleio.output.states import States  # noqa

from .__version__ import __version__  # noqa
