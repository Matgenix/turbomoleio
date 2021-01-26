.. _developer_contributing:

============
Contributing
============

With a few exceptions, that will be described in the following sections, most of the
objects and features present in turbomoleio do not require specific explanations.
A developer should nonetheless comply with the general coding guidelines provided
here.

The most basic aspect concerns the coding in general. turbomoleio is written in **python 3**
and should support versions greater or equal to 3.6. It is important to have an
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
and developer guides whose sources can be found in the ``turbomoleio/docs`` folder. They are
typeset using the `reStructuredText <http://docutils.sourceforge.net/rst.html>`_
and converted to html format using `sphinx <http://www.sphinx-doc.org>`_.

Implementing unit tests and, if suitable, integration tests is also a requirement for the
acceptance of new developments in turbomoleio. See the :ref:`developer_testing` section for more details.

A final remark: since other packages rely on turbomoleio, changes to the current API are strongly
discouraged. In case a backward incompatible change would lead to clear advantages these should be
clearly highlighted and, possibly, the changes should be agreed with maintainers of packages that
rely on turbomoleio.
