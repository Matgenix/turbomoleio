..
    The turbomoleio package, a python interface to Turbomole
    for preparing inputs, parsing outputs and other related tools.

    Copyright (C) 2018-2021 BASF SE, Matgenix SRL.

    This file is part of turbomoleio.

    Turbomoleio is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Turbomoleio is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with turbomoleio (see ~turbomoleio/COPYING). If not,
    see <https://www.gnu.org/licenses/>.

=========
Changelog
=========

Release 1.0.0 (Mar 11, 2021)
============================

First official release (GPL v3 license) of turbomoleio, the python interface to
Turbomole. Turbomoleio is expected to work with Turbomole version 7.3 and higher.
Some parts may work with previous versions as well.

* Automatic definition and modification of inputs
    * running define from a template
    * datagroup-like files
    * control file
* Automatic parsing of outputs (log files) from:
    * dscf
    * ridft
    * escf
    * grad
    * egrad
    * relax
    * statpt
    * aoforce
    * jobex
