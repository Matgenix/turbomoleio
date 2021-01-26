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
