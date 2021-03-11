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

==========
Versioning
==========

turbomoleio adopts a `semantic versioning <https://semver.org/>`_ to identify the different
versions of the packages when released. In particular, given a version number
MAJOR.MINOR.PATCH, increment the:

#. MAJOR version when you make incompatible API changes,
#. MINOR version when you add functionality in a backwards-compatible manner, and
#. PATCH version when you make backwards-compatible bug fixes.

To ease the access to the version number, this is stored as a string in a single
point in the python code, i.e. in the :mod:`turbomoleio.__version__` module.
This file is accessed by the ``turbomoleio/__init__.py`` module, by the ``setup.py``
file and by the ``doc/conf.py`` file used to configure the generation of the documentation.
When creating a new version of turbomoleio you should thus only update the version in
``turbomoleio/__version__.py``.
