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

.. _developer_parse_logs:

===================
Output logs parsing
===================

Most of the details about how the parsing of log output files is explained in the
:ref:`parsing_outputs` section of the user guide. Here we will provide more information
about the different objects involved and what are the constraints that should be
followed to add the parsing of new quantities or of new types of outputs to the
current implementation.

Let us first have a more in-depth inspection of the :class:`turbomoleio.output.parser.Parser`
object. The ``Parser`` takes the string of a TURBOMOLE log file as an input. It is then used to extract
the relevant data information required for each Data or File object. The Parser is made of several
parsing methods, each of which is in charge of parsing a specific portion of the log file.
Some of the parsing methods are very specific to one type of calculation (e.g. only relevant for an escf
output, or for a statpt output), while others are common to several of them (e.g.
:meth:`turbomoleio.output.parser.Parser.basis`, :meth:`turbomoleio.output.parser.Parser.header`).
Each parsing method returns a dictionary with the data parsed or ``None`` if the section to be parsed
could not be found in the string. This dictionary is meant to be used by the Data and File objects
during their instantiation.

The parsing methods are implemented as *lazy* properties, meaning that they will be generated only
once, using the ``lazy_property`` context manager available in the
`monty package <http://guide.materialsvirtuallab.org/monty/monty.os.html#monty.os.cd>`_.
This is advantageous since some sections contain mixed information of the same property and it
may be called more than once while building a File object. The data is thus stored temporarily
in the ``Parser`` instance, but this should not be a problem, because the information
extracted is relatively small and the parser instance is meant to be disposed at the end of the
generation of the object.

In general, when possible, the parsing methods first narrow down the section of the text that
contains the information that should be extracted and then work on this to extract the exact data
needed. This is usually because the outputs of TURBOMOLE do not have an organized structure
and this allows to have a target specific lines of the output more easily.

The information parsed with the ``Parser`` are meant to be used by the Data objects to create
an instance. In particular all the Data and File objects are subclasses of the
:class:`turbomoleio.output.data.BaseData` abstract class. This has a single abstract method,
``from_parser``, that should implemented by subclasses and should be able to create a new instance
of the Data or File object from a ``Parser`` instance. The ``BaseData`` class provides a
``from_string`` and a ``from_file`` method that rely on ``from_parser`` to create an instance of
the subclass from the string of the output or from the file.

The Data objects are meant to contain data that can be grouped together by meaning or by type.
In the ``from_parser`` method it can use the information coming from a single property of the
``Parser`` or from more than one. On the other hand the File objects in the ``from_parser``
preferably do not access directly the properties in the ``Parser`` but instead create instances
of Data object and store mainly Data instances as attributes.

From a design point of view, this whole approach has been chosen since it allows to:

* share parsing functions for different types of calculations,
* quickly instantiate Data objects based on a single portion of the output, e.g. if the user
  only wants information about the energy, he can obtain it running ``ScfEnergiesData.from_file("ridft.log")``
  without the waste of parsing all the rest of the file.
* update the parsing of a specific section for new versions of TURBOMOLE without impacting the
  other parsing methods. This is possible since the parsing methods for the different sections
  are isolated.

.. _developer_parse_new_quantity:

Parsing a new quantity
======================

If you want to start parsing a new quantity for one of the types of log files that are already
handled you should proceed from the bottom up. First check if the quantity that you want to parse
is part of a section that is already parsed by some method of the ``Parser``. If that is the case
you can probably modify that method, otherwise you should create a new property, mark it with
the ``lazy_property`` context manager and implement the parsing inside it. A good approach is
to use regular expressions both to narrow down the section that you want to focus on and to
extract the exact values that you want, but for small sections analyzing the lines one by one
is also acceptable.

.. note::
  When writing the parser always consider that the output of TURBOMOLE can change considerably
  depending on the different options provided (e.g. an entire section or single values might
  appear/disappear). Always check that your parser is working under different conditions.

Once the ``Parser`` has been modified you should either update one (or more) of the Data objects
or write a new one, depending on the type of information that you extracted. In the latter case
subclass the ``BaseData`` class and get your newly parsed information in the ``from_parser`` method.

Finally, if a new Data object has been created, modify the appropriate File objects that should
contain it.


.. warning::

  When adding new attributes to Data and File objects always set a ``None`` as a default
  in the ``__init__``. In this way data that have been json serialized with previous versions
  of the code will still be deserialized correctly (only your newly added property will
  be missing). If instead you change or remove one of the existing attributes the old
  data will not be deserialized anymore. Backward compatibility changes as these should
  be though carefully and agreed upon by the community of users.


Parsing a new type of log
=========================

When parsing the output provided by an executable that is not supported by the current
version of turbomoleio you should proceed in a similar manner as when
:ref:`developer_parse_new_quantity`. In this case you will probably need to add several
new methods to the ``Parser``, addressing the different quantities and information that
you want to extract. Also check which ones of the existing properties are working for
your kind of output, since some of the are likely to be compatible (e.g. the ``header``
attribute will very likely be equivalent in your case as well).

After creating the parser methods, your should encapsulate that information inside
Data objects and create a new File object where you will store all the extracted data.
The same recommendations given in the previous sections hold here as well.

Lastly, if this is suitable, add your object to the :data:`turbomoleio.outputs.files.exec_to_out_obj`
dictionary. This will be used as a reference to decide which File object to use when
running a specific executable. In particular it will be used in the unit tests, as
explained below.


Tests
=====

For the main discussion concerning the testing in turbomoleio you should refer to the
:ref:`developer_testing` section of this developer guide. However, given the particular
nature of the unit tests implemented for the log parsing we will provide some more
explanations here.

The tests for the ``Parser`` object are performed running all the methods implemented
on a series a TURBOMOLE output files. The generated dictionaries are then compared
with references stored in the ``turbomoleio/testfiles/outputs`` folder as JSON files.
A tolerance is allowed, given potentially small differences that can happen while converting
strings to floats, but in general the numerical value should be exactly equivalent to
those parsed. Note that for some combinations of files and methods the output will simply
be ``None``.

In case you want to add a new output file to be tested you should add it in the ``testfiles``
folder and also add it to the list in :data:`turbomoleio.output.tests.test_parser.files_list`.
If instead you are adding one or more new methods to the ``Parser``, remember to add their name
to the list in :data:`turbomoleio.output.tests.test_parser.parser_methods`.

In case you need to generate the reference JSON file again, maybe because you have modified
one of the existing ``Parser``'s method or because you added a new one, you can use the
:func:`turbomoleio.output.tests.test_parser.generate_files` helper function. This will
generate a new JSON and **overwrite the old one** for all the files and methods that have
been given in input. So you should be extremely careful when running it, since the generated
files will become the new reference. If a bug is introduced in the ``Parser``, the reference
files will be generated with a bugged version and the tests will partially loose their use.

A similar approach has been chosen for the testing of the Data and Files objects. Since
all the Data objects are contained in at least one File object, the tests concerning the parsing
will be performed only at the level of the File objects, since repeating them for the Data as
well would just be redundant.

The structure of the tests is similar to those for the ``Parser``. The test output files are
parsed using the ``from_file`` method of the corresponding File object. The object is converted
to a dictionary and compared with the reference stored in the corresponding JSON file.

As before, if you want to add a new output file to be parsed, you should add it in the
``testfiles`` folder and in the list in :data:`turbomoleio.output.tests.test_files.files_list`.
In addition if you want to add a new type of File object you should either add it to the
:data:`turbomoleio.outputs.files.exec_to_out_obj` or make it available in the
:func:`turbomoleio.output.tests.test_files.cls_dict_path` fixture (follow the example
in the case of ``EscfOnlyOutput`` there).

A function to generate the reference JSON files similar to the one described above is available:
:func:`turbomoleio.output.tests.test_files.generate_files`. The same warnings of dealing
with it carefully should be kept in mind here as well.
