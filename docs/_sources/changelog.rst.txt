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

Release 1.1.1 (Feb 20, 2022)
============================

Features added
--------------

* Handling of periodic systems for the riper module.
    * Definition of a PeriodicSystem object to handle 1, 2 or 3-dimensionally periodic
      structures.
    * Adaptation of the Gradient object to extract the lattice and lattice gradients
      in addition to the atomic coordinates and gradients.
    * Note that this feature is currently experimental, meaning that it is a work in
      progress and the API might change without prior notice.
* Parsing of periodic calculations performed with riper [preliminary].

Release 1.1.0 (Sep 20, 2021)
============================

Features added
--------------

* Parsing of MP2 calculations performed with mpgrad and rimp2/ricc2 [preliminary].

Bug fixing
----------

* Corrected parsing of convergence information in statpt and relax programs.

Testing and development
-----------------------

* Reorganization of test files for the output parsing
    * A generation folder has been created in the testfiles folder. This folder
      contains all the available tests and the required information to
      regenerate reference output files and objects for a new Turbomole versions.
    * Reference Turbomole output files and serialized data objects are stored in
      specific version-related folders.
* New development script to regenerate reference outputs, data and file objects for the
  testing of the output parsing.

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
