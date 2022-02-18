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

.. _developer_contributing:

============
Contributing
============

With a few exceptions, that will be described in the following sections, most of the
objects and features present in turbomoleio do not require specific explanations.
A developer should nonetheless comply with the general coding guidelines provided
here.

The most basic aspect concerns the coding in general. turbomoleio is written in **python 3**
and should support versions greater or equal to 3.8. It is important to have an
understanding of the basic principles of object orientation programming, like
inheritance and abstraction, since these are often used in the code.
The coding follows mainly the conventions defined in the standard
`PEP 8 <https://www.python.org/dev/peps/pep-0008/>`_ style guide. Although this is not
strictly enforced and deviations are acceptable, it would be advisable to stick to it
when possible.

Providing detailed **in-code documentation** in accordance with the standards is also essential
to help future developers and to automatically generate the API documentation. In turbomoleio
the `Google style <https://github.com/google/styleguide/blob/gh-pages/pyguide.md#38-comments-and-docstrings>`_
is used for docstrings and comments. New developments should also be described in the user
and developer guides whose sources can be found in the ``turbomoleio/docs_rst`` folder. They are
typeset using the `reStructuredText <http://docutils.sourceforge.net/rst.html>`_
and converted to html format using `sphinx <http://www.sphinx-doc.org>`_. The html files are located
in the ``turbomoleio/docs`` folder and can be generated from the source rst files using the
`invoke <https://www.pyinvoke.org/>`_ `make-doc` task defined in
`~turbomoleio/tasks.py`.:

    invoke make-doc

Implementing unit tests and, if suitable, integration tests is also a requirement for the
acceptance of new developments in turbomoleio. See the :ref:`developer_testing` section for more details.

A final remark: since other packages rely on turbomoleio, changes to the current API are strongly
discouraged. In case a backward incompatible change would lead to clear advantages these should be
clearly highlighted and, possibly, the changes should be agreed with maintainers of packages that
rely on turbomoleio.
