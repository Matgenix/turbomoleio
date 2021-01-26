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
