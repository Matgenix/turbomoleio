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

.. _running_define:

==============
Running define
==============

A crucial point of the execution of TURBOMOLE is generating the input through the ``define`` executable.
Since ``define`` is executed interactively this step is challenging when aiming at automating
the calculations with TURBOMOLE.
This problem is solved with the :class:`turbomoleio.input.define.DefineRunner` object by automating
the interaction with ``define`` with the `pexpect <https://pexpect.readthedocs.io>`_ module, that spawns
child applications and responds to expected patterns in their output. This can be used inside
some complex workflow, but also to speed up the manual execution of ``define``.

``DefineRunner``, based on the parameters provided in input as a dictionary, runs ``define`` and
navigates through the menus, setting the different options. It also stores all the history of the
submitted commands in the ``sent_commands`` attribute and the comments about the commands and the expected
replies from the code in the ``history`` attribute.

In case of failure an exception will be raised, that can be :class:`turbomoleio.input.define.DefineError`
or one of its subclasses, depending on the type of error encountered. In particular, it should be noted
that for each submitted command code will expect some specific string coming from ``define``, if this
does not happen within the number of seconds specified by the ``timeout`` argument a
:class:`turbomoleio.input.define.DefineExpectError` is raised. This can happen mainly for two reasons,
the most common is that a combination of the files already available and/or of the options selected
brought the execution of ``define`` in a situation that was not foreseen. However this can also happen
in case of really low responsiveness from the system, taking too much time to provide the output.
Usually 60-120 seconds is a good compromise, but in case you encounter this failure and you cannot
figure out the correct reason, you might want to check if there is some system problem by increasing the
value of the ``timeout``.

Given an instance of ``DefineRunner`` there are three methods that can be called and will result in
different kind of executions of ``define``:

* ``run_full()`` will generate the whole set of inputs running all the menus in ``define``. It expects
  to only have the ``coord`` file in the working directory. If a ``control`` file is already present
  this will likely lead to some ``DefineError``.
* ``run_update_internal_coords()`` runs ``define`` just to update the internal coordinates of a set
  of inputs previously generated. In this case the ``control`` file should already be present in the
  working directory. In this case the execution is interrupted after going through the molecular geometry
  menu, so ``define`` ends *abnormally* (by sending the ``qq`` command).
* ``run_generate_mo_files()`` is used to regenerate the molecular orbitals files in a folder where
  some inputs are already present. As for the previous case the ``control`` file should be already present
  and ``define`` will be terminated *abnormally* after the molecular orbital menu.

In order to help new users to start using ``DefineRunner`` a small set of predefined input parameters
for the most common type of calculations have been prepared. These are stored as YAML files in the
``turbomoleio/input/templates`` folder, but can also be easily accessed with the
:func:`turbomoleio.input.utils.get_define_template` function. A typical execution of ``DefineRunner``
in a folder already containing a ``coord`` file can be done in this way

.. code-block:: python

    from turbomoleio.input.define import DefineRunner
    from turbomoleio.input.utils import get_define_template

    dp = get_define_template("ridft")

    dr = DefineRunner(parameters=dp)
    dr.run_full()

This will prepapare the inputs for an ``ridft`` calculation. Of course you can tune the ``parameters``
according to the type of calculation that you want to perform. A good way of proceeding would be
to prepare the dictionary with options that you need and store it as a YAML or JSON file, so
that it can be easily retrieved afterwards to generate other calculations.

In general, even though it does not cover every single scenario available in ``define``,
``DefineRunner`` supports a large variety of options and can produce the inputs for different types
of calculations. Anyway it is focused on providing options for those inputs that can only be set through
define. Simple data groups can be set after the execution of ``DefineRunner`` by the user with the
standard functions available to modify a data group file (e.g. with the ``cdg`` function.
See :ref:`mod_data_group` and :ref:`control_object`).

The DefineRunner parameters
===========================

In this section we provide a full list of all the parameters that can be passed to the ``parameter``
argument when creating an instance of ``DefineRunner``.

* ``title`` (:py:class:`str`): title of the job passed to ``define``.
* ``metric`` (:py:class:`int`): sets the ``$metric`` keywork in the ``control`` file before running ``define``.
* ``copymo`` (:py:class:`str`): path to a directory containing the ``mos``, ``alpha`` and ``beta``
  files that will be copied in the current working directory at the end of ``define``. The control
  file should be present in the folder as well, since it will be used to extract the value of the symmetry,
  overriding the ``sym`` and ``desy`` values.
* ``sym`` (:py:class:`str`): the value will be passed to force a specific symmetry with the ``sy``
  command in the molecular geometry menu.
* ``sym_eps`` (:py:class:`float`): if present will be added to the ``sym`` option as a tolerance for
  the ``sy`` command in the molecular geometry menu.
* ``desy`` (:py:class:`bool`): the system will determine the symmetry of the molecule in the molecular
  geometry menu. Only used if ``sym`` is not defined.
* ``desy_eps`` (:py:class:`float`): if present will be added as a tolerance for the ``desy`` command
  in the molecular geometry menu.
* ``ired`` (:py:class:`bool`): generates the internal coordinates with the ``ired`` command.
* ``usemo`` (:py:class:`str`): path to a ``control`` file or to a directory containing it. The file
  will be passed to ``define`` with the ``use`` command.
* ``ex_method`` (:py:class:`str`): method used to calculate the excited states. Available options:
  ``rpa``, ``cis``, ``dynpol``, ``polly``.
* ``ex_multi`` (:py:class:`str`): multiplicity of the excited states. Available options: ``singlet``, ``triplet``.
  This will be only applied for closed shell calculations to distinguish between ``rpas/ciss`` and ``rpat/cist``.
  The type of calculation will be determined according to the options available for the calculations based on the
  outcome of the EHT. The value of ``ex_multi`` will be ignored if the calculation is an UHT or if ``ex_method``
  is ``dynpol`` or ``polly``.
* ``ex_all_states`` (:py:class:`int`): the number of excited states for all the irreps present in the
  system according to ``define``. Since ``define`` does not accept values larger than the number
  of states available, the code will check all the available states and set the number to the minimum
  between ``ex_all_states`` and the number of actual available states for each specific irrep.
* ``ex_irrep_states`` (:py:class:`dict`): a dictionary of the form ``{"irrep": num_exc_states}``, with
  the key representing the irrep and the value the number of excited states (e.g. ``{"a1": 10}``).
  It will override the values of ``ex_all_states`` for the specific states mentioned in this dictionary.
* ``ex_mp2`` (:py:class:`dict`): list of excited states for the *mp2* calculations. Should be a dict of
  the form ``{"irrep": [multiplicity_int, num_exc_states]}``, with the key representing the irrep and the
  value a list with an :py:class:`int` representing the multiplicity (1=singlet, 2=doublet, 3=triplet)
  and one one representing the number of excited states (e.g. ``{"a1": [1, 10]}``). This will be provided
  to the ``exci`` command with  ``irrep=a1 multiplicity=1 nexc=10``.
* ``ex_frequency`` (:py:class:`float`): value of the frequency for the calculation of dynamic polarisabilities.
  Default 589 nm if not defined and excited states are required.
* ``ex_frequency_unit`` (:py:class:`str`): units of the frequency for the calculation of dynamic
  polarisabilities. Default nm.
* ``ex_exopt`` (:py:class:`int`): explicitly enforces treatment of the n-th state. Set directly in the
  ``control`` file with the ``$exopt`` keyword using ``cdg`` while ``define`` is being executed.
* ``method`` (:py:class:`str`): calculation method. Available options ``dft``, ``hf``, ``mp2``, ``adc(2)``,
  ``ccsd(t)`` (or ``ccsdt``).
* ``mp2energy`` (:py:class:`bool`): if ``True`` the calculations will be limited to the energy
  (i.e. in the ``ricc2`` menu only the ``method`` will be provided. Otherwise the ``geoopt method``
  option will be given to ``define``.
* ``basis`` (:py:class:`str`): the basis that will be used for all the elements in the system
  (e.g. ``b all def2-SV(P)``).
* ``basis_atom`` (:py:class:`dict`): a dictionary of the type ``{"atom": "basis"}`` defining the
  basis for specic atoms, where the key of the dictionary can be any string accepted by ``define``
  (e.g. ``1,2,4-6``, ``c``). If ``basis`` also defined, it will be used to set the basis for all the
  atoms and  then ``basis_atom`` will override specific atoms.
* ``charge`` (:py:class:`int`): charge defined in the extended Hueckel guess.
* ``unpaired_electrons`` (:py:class:`int`): number of unpaired electrons for UHF.
* ``rijk`` (:py:class:`bool`): activates rijk calculation.
* ``ri`` (:py:class:`bool`): activates ``ri`` calculation, only if ``rijk`` is disabled.
* ``marij`` (:py:class:`bool`): activates ``marij`` calculation, only if ``rijk`` is disabled and
  ``ri`` is enabled.
* ``functional`` (:py:class:`str`): functional for DFT. All available values in TURBOMOLE.
* ``gridsize`` (:py:class:`str`): size of the grid in DFT calculation.
* ``maxcor`` (:py:class:`float`): memory set in *mp2* calculations.
* ``use_f12`` (:py:class:`bool`): enables *f12* calculation for methods ``mp2``, ``adc(2)`` and ``ccsd(t)``.
* ``use_f12*`` (:py:class:`bool`): enables *f12\** calculation for method ``ccsd(t)``. Adds the line
  ``ccsdapprox  ccsd(f12*)`` to the ``$rir12`` data group. Requires ``use_f12``.
* ``maxiter`` (:py:class:`int`): maximum number of iterations for *mp2* calculations.
* ``scfiterlimit`` (:py:class:`int`): maximum number of scf iterations.
* ``scfconv`` (:py:class:`int`): accuracy of scf energy. A number in the range 4-9.
* ``coord_file`` (:py:class:`str`): path to the ``coord`` file. By default uses the ``coord`` file in
  the working folder.
* ``disp`` (:py:class:`str`): activates dispersion correction according to the provided value. Accepted
  values are ``DFT-D1``, ``DFT-D2``, ``DFT-D3``, ``DFT-D3 BJ``. N.B. the values will be set directly on
  the control file, not using ``define``.

In addition the following keywords are related to cosmo and are set directly in the ``control`` file after
``define`` has completed. The functionality is activated setting the ``use_cosmo`` in the parameters.
This adds the ``$cosmo`` data group to the control file with the following options plus the
``$cosmo_out = out.cosmo`` data group. The additional options will modify the values inside the
``$cosmo`` data group.

* ``use_cosmo`` (:py:class:`bool`): If True enables the calculation with cosmo.
* ``epsilon`` (:py:class:`float`): permittivity used for scaling of the screening charges.
* ``nppa`` (:py:class:`int`): number of basis grid points per atom.
* ``nspa`` (:py:class:`int`): number of segments per atom.
* ``disex`` (:py:class:`float`): distance threshold for A matrix elements (Angstrom).
* ``rsolv`` (:py:class:`float`): distance to outer solvent sphere for cavity construction (Angstrom).
* ``routf`` (:py:class:`float`): factor for outer cavity construction in the outlying charge correction.
* ``cavity`` (:py:class:`str`): acceptable values are "open" (leave untidy seams between atoms) and
  "closed" (pave intersection seams with segments).
* ``use_old_amat`` (:py:class:`bool`): if True adds the ``use_old_amat`` to the ``$cosmo`` data
  group, i.e. uses A matrix setup of TURBOMOLE 5.7.

Validation
----------

An **experimental feature** is available to validate the correctness of the dictionary that you want
to use as ``paramaters`` for ``DefineRunner``. It is based on the
`cerberus <https://docs.python-cerberus.org/en/stable/>`_ package and is implemented in the
:mod:`turbomoleio.input.utils` module.

The simplest way to use it is:

.. code-block:: python

    from turbomoleio.input.utils import get_define_template, validate_parameters

    dp = get_define_template("dscf")
    dp["use_cosmo"] = True
    validate_parameters(dp)

This will return ``True`` if the passed dictionary is correct according to the defined schema
and ``False`` otherwise. At the moment the function will provide the following kind of validations:

* Check that all the keys of the dictionary are acceptable ones. This should prevent
  inserting typos in the keys.
* Check that the values are of the correct type, according to the definitions above.
* Validate some dependencies among the different options. For example ``use_f12*`` can be
  ``True`` only if ``use_f12`` is also ``True``.
* Validate the values of some options. For example that ``method`` is among one the allowed values:
  ``dft``, ``hf``, ``mp2``, ``adc(2)``, ``ccsd(t)``, ``ccsdt``

None of the parameters is required and thus an empty dictionary will be considered as valid.

Alternatively you can directly access :data:`turbomoleio.input.utils.define_parameters_validator`,
that is an instance of a cerberus ``Validator``, so that you can take full advantages of its features.

Being experimental this feature has not been extensively tested, and you should use it with a bit
of care. If the validation of one of your dictionaries fails but you are absolutely certain that
it is correct you can probably ignore the failure in the validation. In addition for some of the
parameters a validation of the possible values is missing. Just to mention one, the validation
only checks that the value of ``functional`` is a string, but no check is performed to verify
that the value is the name of an existing functional.
