..
    The turbomoleio package, a python interface to Turbomole
    for preparing inputs, parsing outputs and other related tools.

    Copyright (C) 2018-2022 BASF SE, Matgenix SRL.

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

.. _developer_versioning:

==========
Versioning
==========

turbomoleio adopts a slightly modified `semantic versioning <https://semver.org/>`_ to identify
the different versions of the packages when released. In particular, given a version number
MAJOR.MINOR.PATCH, increment the:

#. MAJOR version when you make incompatible API changes,
#. MINOR version corresponds to a specific version of Turbomole (e.g. 1.0.x is for Turbomole
version 7.3, 1.1.y is for Turbomole version 7.4, ...), and
#. PATCH version when you add new functionalities in a backward-compatible manner or make
backwards-compatible bug fixes.

To ease the access to the version number, this is stored as a string in a single
point in the python code, i.e. in the :mod:`turbomoleio.__version__` module.
This file is accessed by the ``turbomoleio/__init__.py`` module, by the ``setup.py``
file and by the ``doc/conf.py`` file used to configure the generation of the documentation.
When creating a new version of turbomoleio you should thus only update the version in
``turbomoleio/__version__.py``.

======================
Backward-compatibility
======================

The latest version of turbomoleio supports the latest version of Turbomole. The latest version
thus supports the generation of input files using the `define` executable, the reading, parsing and
modification of datagroup-like files, as well as the parsing of the output files (log files)
of the latest version of Turbomole. In addition, the latest version of turbomoleio (corresponding
to the latest version of Turbomole) also supports the parsing of the output files generated with
older versions of Turbomole (starting at version 7.3).

For example, turbomoleio 1.2.0 supports Turbomole version 7.5 for all the features and supports
parsing of outputs for Turbomole version 7.3 and 7.4.
