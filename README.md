Installation
============

Dependencies
------------

For testing the following modules should be installed:
pytest
pytest-cov

Usage
=====

Developing
==========

tests are run via the runtests script in the root directory, via py.test

* `runtests.sh -u` runs only the units tests (but requires TM as well to be installed)
* `runtests.sh -i` runs only the integration tests
* `runtests.sh` runs all the test

Use `coverage html -d coverage_html` afterwards to generate html reports of the coverage.

License
=======

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
