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

=============
Define runner
=============

The description about the main features of :class:`turbomoleio.input.define.DefineRunner`
has been given in the corresponding section of the user guide: :ref:`running_define`.
Here we will briefly discuss a few implementation aspects that might be useful in case
you want to add some more option to the execution of ``define``.

.. note::

  Remember that if the option that you want to add is a trivial one,
  i.e. that could be simply added with a ``change_data_group``, you should probably consider
  to set it in this way, instead of passing it through ``define``, which will be much
  more complicated and error prone to implement, while bringing very little advantage, if any.
  For particular data groups an alternative possibility is to add a helper method
  in the ``Control`` object. See for example :meth:`turbomoleio.core.control.Control.set_disp` or
  :meth:`turbomoleio.core.control.Control.add_cosmo`.

If you wish to add a new call to one of the options in the ``define``'s menu you should first
explore the sequence of operations that are performed by ``DefineRunner``. The best way
would be to start from the ``run_full`` method and follow the procedure. You can thus
identify the most suitable point where your option should be called and start to plan how
to call it. The best option is probably to create a new protected method (i.e. starting
with a ``_``) and that will be called from one of the other methods at the right moment.

This new option will probably need an additional input parameter. You should make sure that
if your option is turned off your new portion of code is skipped.

.. note::

  If you introduce a new option in the list of parameters passed to ``DefineRunner``
  always remember to carefully document it in the user guide.

Inside your method you can send one or more commands to ``define`` in order to set the
options that you desire. There are a few points that you should keep in mind:

* **Avoid accessing directly to the ``self.define`` instance of ``spawn``**. This is
  acceptable if you need to do some advanced parsing of the returned values, but
  in any case you should always pass through the ``DefineRunner._sendline`` and
  ``DefineRunner._expect`` when sending commands and expecting the replies.
  This automatically keeps track of all the operations in the ``history`` and
  ``sent_commands`` attributes.
* After sending a command **always make an expect**, even a relatively trivial one
  with only one option, to ensure that the command has been received.
* If there is any potential known failure that can happen consider to intercept it
  (for example by "expecting" a specific string that would signal the problem)
  and raise a ``DefineError`` or one of its subclasses, if suitable.
* Make sure that your new option is not disrupting the flow of any of the ``run_``
  methods, since the methods are shared by all of them.

In any case it is always a good idea to try to mimic the behavior of other methods
already implemented.

If instead you need implement a completely new functionality, that runs define
only performing a specific set of operations in order to modify an already existing
set of input files, you should take the ``run_update_internal_coords`` and
``run_generate_mo_files`` as examples and try to implement the full set of operations
inside a new method. Try to reuse the available internal methods if possible.


Tests
=====

For the main discussion concerning the testing in turbomoleio you should refer to the
:ref:`developer_testing` section of this developer guide. However, given the particular
nature of the tests implemented for the various functionalities of ``DefineRunner``
we will provide some more explanations here.

The unit tests for ``DefineRunner`` rely heavily on mocking, since no call to ``define``
can be done. The ad-hoc ``dr_data`` fixture is available, providing all that is needed
to mock all the values returned by the ``expect`` command.

When adding a new feature you should introduce unit tests for the different options
and paths that the code could follow. In addition to this, given the interconnectedness
of the different methods, it is possible that some tests will start to fail after the
introduction of your new feature. In that case it is your responsibility to fix them
and make sure that all succeed before trying to add your new feature to the main
distribution of turbomoleio.

Concerning the integration tests, you should introduce one or more tests
where your new option is used, so that it will also be checked against a real execution
of ``define`` and that is producing the correct inputs for some TURBOMOLE calculation.

Validation
==========

The validation implemented in the :mod:`turbomoleio.input.utils` module has already been
described in the user guide. Being an experimental feature you are welcome to contribute to
it and improve its level of validation. Also, if you introduce a new option for
``DefineRunner`` you should at least add the name and the type to the
:data:`turbomoleio.input.utils.schema_define_params` dictionary, that contains the
validation schema according to the cerberus conventions.

The `cerberus <https://docs.python-cerberus.org/en/stable/>`_ is a powerful tool that
validates a dictionary based on a schema. The best place to learn about all of its options
is is directly in the cerberus documentation.

Since in ``DefineRunner`` different methods can be called, some of which targeting only a
specific subset of the ``define`` menus, you should also make your new attribute as ``nullable``.
If possible try to include the list of ``allowed`` values (or try to add it for the parameters
already defined) and add the required dependecies, if any.
