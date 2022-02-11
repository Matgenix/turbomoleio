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

==============
The coord file
==============

In TURBOMOLE the coord file mainly contains the coordinates and the types of the atoms
of the molecule that should be simulated. These information can be stored in
`pymatgen <http://pymatgen.org/>`_ ``Molecule`` object. However there are a set
of data, stored in different data groups, that provide further information
about the geometry of the molecule, more precisely on the dynamic during a geometry
optimization. For this reason in turbomoleio the ``coord`` file is represented
through a :class:`turbomoleio.core.molecule.MoleculeSystem` object. This has methods
to read and generate a ``coord`` file that can be read by TURBOMOLE. Other file formats
can be written taking advantage of the conversion from ``Molecule``, but the output will be
limited to the geometry of the molecule.

.. note::

    TURBOMOLE is also capable of performing periodic calculations (1, 2 or 3-dimensional).
    The :class:`turbomoleio.core.molecule.MoleculeSystem` object is only used for isolated
    systems, i.e. molecules. For periodic systems, the
    :class:`turbomoleio.core.periodic.PeriodicSystem` should be used.
    See :ref:`Periodic systems<periodic_systems>` below for more information.

Internally the ``MoleculeSystem`` stores the geometry of the molecule using a
pymatgen ``Molecule``. This allows to exploit all the functionalities available in
the pymatgen package. Other attributes are used to determine the values of other data
groups. The simplest one is ``frozen_indices`` that determines which atoms should have the
cartesian coordinates frozen. Translated in the ``coord`` file format the atoms with an ``f``
at the end of the coordinate line. For example with ``frozen_indices={0,2}`` the coord file
will look like::

    $coord
       -0.00000000000000      0.00000000000000      0.00000000000000  n  f
       -1.17616323936717      2.03717448857876      0.00000000000000  o
       -1.17616323936717     -2.03717448857876      0.00000000000000  o  f
        2.35232647873438      0.00000000000000      0.00000000000000  o
    $end

.. warning::

    For all the properties set in ``MoleculeSystem`` that refers to the indices of
    the atoms in the ``Molecule`` the convention is to use **0-based indices**.
    This is true for ``frozen_indices``, ``int_def`` and ``user_defined_bonds``.
    The values used then in the ``coord`` file are 1-based instead. The conversion
    is done internally and you should always use 0-based indices in the python code.

An important and more complicated attribute is ``int_def`` that, if present, is used to
fill the ``$intdef`` data group in the ``coord`` file. The attribute should be a list of
internal coordinates given as subclasses of :class:`turbomoleio.core.molecule.InternalDefinition`.
Each subclass define a type of internal variable as defined from TURBOMOLE. For example the
:class:`turbomoleio.core.molecule.Distance` represent an internal variable of type ``stre``.
For the internal definitions that are most commonly used (``stre``, ``bend`` and ``tors``) helper
methods are provided in ``MoleculeSystem`` to append additional internal definition to an existing
instance. For example the following python code:

.. code-block:: python

    m = MoleculeSystem.from_file("coord")
    m.add_distance(atom1=0, atom2=1, value=1.5)
    m.add_bond_angle(atom1=1, atom2=2, vertex=3, status="k")


will result in a coord file with this structure::

    $coord
        -0.00000000000000 0.00000000000000 0.00000000000000 n
        -1.17616323936717 2.03717448857876 0.00000000000000 o
        -1.17616323936717 -2.03717448857876 0.00000000000000 o
        2.35232647873438 0.00000000000000 0.00000000000000 o
    $intdef
        1 f 1.0 stre 1 2 val=1.5
        2 k 1.0 bend 2 3 4
    $user-defined bonds
        1-2, 2-4, 3-4
    $end


We stress again that the indices provided by the are 0-based as opposed in the one produced in the
file, that are 1-based.

.. warning::

    Following pymatgen conventions the ``Molecule`` object contained in ``MoleculeSystem`` has coordinates
    in angstrom, while the ``coord`` file contains coordinates in Bohr. The values provided for the
    internal definitions in the helper methods and in the ``InternalDefinition`` subclasses are written
    in ``$intdef`` as they are provided. Their units are the same as those defined in TURBOMOLE for
    that specific type of internal coordinate.

Note that the helper methods have also automatically set the ``$user-defined bonds``. This brings us
to the last attribute that can be set in ``MoleculeSystem``: ``user_defined_bonds``. This can contain
a set of tuples of the form ``(index1, symbol, index2)``, where the atom indices are 0-based
and the symbol can be ``-`` or ``|``.

These attributes could be used to set the complete list of internal coordinates required by TURBOMOLE,
but in general this would be more useful to define a few of them and call ``define`` afterwards with
the ``ired`` option to let it generate a complete list of redundant internal coordinates. Given its
potentially complicated structure, the ``$redundant`` data group is never parsed nor stored when
reading a ``coord`` file that contains it.

.. _periodic_systems:

Periodic systems
================

.. warning::

    This is an experimental feature in turbomoleio. As a work-in-progress, the API might still
    change without prior notice.

For periodic systems, the ``coord`` file and the periodicity is represented
through a :class:`turbomoleio.core.periodic.PeriodicSystem` object. This has methods
to read and generate a ``coord`` file that can be read by TURBOMOLE, as well as the corresponding
periodic attributes.
