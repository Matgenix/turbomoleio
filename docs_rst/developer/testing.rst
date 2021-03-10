.. _developer_testing:

=======
Testing
=======

Tests are an essential part of the development process. Their aim is to check that the
different portions of the code meet their design and behave as intended. Providing an exhaustive
guide about how to write tests in python is beyond the scope of this document, so in this section we
will instead explain the details of the tests in turbomoleio: how they are organized,
what are the tools available and how you should proceed to introduce tests for your new developments.

In turbomoleio we divide the tests in two groups. The definition could be considered somewhat
arbitrary with respect to the standard ones, but in the following we use the following naming convention:

* **unit tests**: aiming at testing each single component on its own (no TURBOMOLE executable
  runs inside these tests).
* **integration tests**: testing the integration among different components and in particular
  with the TURBOMOLE executables.

When implementing a new feature you should always implement a series of tests before
adding it to turbomoleio. It is also important to consider adding new tests when a bug is identified
and fixed. This will prevent that such a bug will be reintroduced in later developments.

.. _developer_test_suite:

The test suite
==============

The testing suite is based on **pytest** and makes wide use of its functionalities, like for example
fixtures, marking and custom command line arguments. The pytest package should thus be installed
and the tests should be executed using the ``pytest`` command. For a reference on all its features you
can refer to the `pytest documentation <https://docs.pytest.org>`_. To get the coverage of the
tests you also need to install the **pytest-cov** package.

Summing up, aside from the standard dependencies in BASFflows, the easiest way to prepare your environment
to run the tests is to run

.. code-block:: bash

    pip install pytest pytest-cov mongomock

To run the integration tests of course a standard installation of Turbomole should be available, so
that the tests can call the different executables.

Once this is done you can run each test separately going to the different folders and executing
the ``pytest`` command or use the ``runtest.sh`` script present in the root folder of the project.
Running ``runtest.sh -u`` will run all the unit tests, while ``runtest.sh -i`` will run all
the integration tests. In both cases it will produce a detailed report of the level of coverage.
If you wish, you can then run

.. code-block:: bash

    coverage html -d coverage_html

to generate an html version of the coverage report. The configurations for the coverage are given
in the ``.coveragerc_unit`` and ``.coveragerc_integration`` files, present in the root folder of
the project.


Unit tests
==========

As mentioned above, we call unit tests those that test each single component on its own.
Each module folder in turbomoleio contains a ``tests`` subfolder where the unit tests for
all the components in that module are implemented.
All the functions, objects and their methods have some unit test that will verify that
it behaves as expected. This is done, for example, by checking:

* the returned values of functions and methods,
* the attributes set in objects,
* the side effects of a function (e.g. creation or modification of files),
* that the expected exceptions are raised under specific conditions,

These tests will be limited to a specific component and thus may require to
`mock <https://docs.python.org/3/library/unittest.mock.html>`_ functions or objects to determine
a specific outcome of the test (see `this guide <https://realpython.com/python-mock-library/>`_
for some examples). Mocks are also required for
all the calls towards external applications and in particular to the TURBOMOLE executables.
Since in the turbomoleio code the calls to TURBOMOLE executables is done for the ``DefineRunner``
class, specific mocks will be done for its testing.

For all the cases where inputs and outputs should be read or modified a set of prepared files are
used. These files, mainly Turbomole inputs and outputs, are stored in the
``turbomoleio/testfiles`` folder.
Using these test files and with the additional help of the mocking it is possible to direct and
verify in detail what happens inside the component under test. For this reason the unit tests can achieve
a very high coverage rate, of the code, where by coverage we intend the fraction of the lines of
code that are executed by at least one test. As of version 1.0.0 the global coverage of the unit
tests alone is 99% and developers should aim at keeping this value as high as possible.


Integration tests
=================

Integration tests aim to verify the integration between the different components and
in particular with the Turbomole executables. For this reason they are all collected in the
``turbomoleio/integrations_tests`` folder and, to allow to execute them separately, they
are all *marked* by setting the ``integration`` pytest mark for each module:

.. code-block:: python

    pytestmark = pytest.mark.integration

Since, aside from in ``DefineRunner``, in the turbomole source code there are no explicit
calls to TURBOMOLE executable, the integration tests will also check if the I/O features implemented
keep working in the same way when the files are directly used or produced by TURBOMOLE.
To be more explicit one can check if a file that has been produced or modified with
turbomoleio keeps being read correctly in TURBOMOLE and if an output file produced
by an execution of TURBOMOLE is still parsed correctly. For some functionalities anyway
there is not really a meaning in performing an integration test. In these cases the
portion of the code can be excluded from the coverage analysis.

The integration tests are all performed using the :func:`turbomoleio.testfiles.utils.run_itest`
helper function. This first runs ``define`` using ``DefineRunner`` based on the parameters given
in input, the produced control file is compared with a reference stored previously. If successful
the required turbomole executables are run and the outputs are extracted using the objects
implemented here. The numerical values of these outputs are checked with a respect to
references stored in JSON files. Note that the comparison is performed using the
:func:`turbomoleio.testfiles.utils.assert_almost_equal`. For this comparison some values
are outright excluded from the check. In some cases a comparison will be obviously meaningless
(e.g. the date of execution, the elapsed time), while in some other cases the values,
while physically and chemically equivalent, may be reported in slightly different ways in
the TURBOMOLE output and it would be extremely impractical to make a meaningful
comparison for them (e.g. the  order of the acoustic modes from aoforce
can vary on different machines, even for the same version of TURBOMOLE. While checking the
list of eigenvalues is trivial, verifying that the two lists of eigenvectors, represented
in different order, are indeed equivalent up to a numerical tolerance would be cumbersome).

This procedure is repeated for different kinds of inputs and for the different
TURBOMOLE executables currently supported by the turbomoleio objects.

Aside from preventing the introduction of errors when modifying existing parts of the code,
integration tests can thus be used to verify that the implemented I/O functionalities
will keep working for **new versions of TURBOMOLE**. It is inevitable that when moving
from one TURBOMOLE version to another at least some of the integration tests will fail. This
might be for some changes in the inputs, in the underlying implementations or even for some
options being removed or renamed altogether. Another possibility is that the output given
in the ``stdout`` has changed its format and that the parser now fails to correctly
extract the results of the calculation. turbomoleio has been initially developed to
target TURBOMOLE version 7.3 and the tests are designed for that specific version only.
As new versions of TURBOMOLE are released the tests can be updated to match the new input and outputs
definitions, but care should be taken by checking exactly which part of the test is leading
to a failure.

In light of all this, two custom options have been added to pytest in turbomoleio. The first is a way
to tune the value of the tolerance when comparing numerical values. This can be changed by
running the tests with the ``--itest-tol``. For example running

.. code-block:: bash

    pytest --itest-tol=0.01

performs a much looser comparison (note that this is an absolute tolearance).
The second option is ``--generate-itest-ref``. When running
the tests with this option, instead of using the reference files to check the
correctness of the outputs produced by the tests, the outputs will instead be used to generate
a new version of the JSON reference files and overwrite the previous ones.

.. warning::

    The ``--generate-itest-ref`` should be used with extreme care. The original files
    provided with turbomoleio have been tested across different machines
    to verify that the tests pass on different environments. This option should
    be used only if the output objects parser are changed or if the reference version of TURBOMOLE
    is updated. In any case it would be wise not to run the whole set of tests with this option,
    but instead target only the specific test for which the reference file should be updated.
    For example:

    .. code-block:: bash

        pytest --generate-itest-ref test_dscf.py -k test_run_dscf_hf

Writing integration tests that trigger all the possible cases in the code is basically impossible.
Some errors and checks in the source code are only triggered in presence of exceptional errors
(e.g. wrongly formatted output files) and this cannot be easily reproduced systematically.
In addition, while there is no strict limit for the execution time of a single test, these should be
short enough to be executed in a reasonable amount of time with small computational resources.
Lastly for some functions and objects there is little gain in testing their behavior in direct connection
with an actual TURBOMOLE calculation. For all these reasons the aim for the coverage of the integration
tests is lower compared to the unit tests. As of version 1.0.0 the integration tests alone cover 89% of
the turbomoleio source code.


Writing new tests
=================

New contributions to turbomoleio should always come with their set of unit and, if suitable, integration tests.
The tests that you add should be organized in the same way as the other tests already available,
as described in the previous sections.

The test functions and files should be named starting with ``test_`` to be correctly discovered by
pytest (this can be customized, but sticking to the standard is the easiest option). If you have some
experience with the standard python testing utils, notice that in turbomoleio you should *not*
subclass ``unittest.TestCase``, since this is incompatible with some of the functionalities in
pytest. All the options provided by ``TestCase`` can be easily replaced by pytest options and
fixtures.

As for the other developments, you can use the tests already present as a reference for implementing
your own. Several resources are available online discussing the testing best practices and the
pytest documentation is a good starting point as well. Here we will limit to the description of
the options and utils specific to turbomoleio.

Helper functions and utilities specific for testing can be found in the :mod:`turbomoleio.testfiles.utils`
module. In particular this contains a context manager that is used in most of the tests present
in the test suite: :func:`turbomoleio.testfiles.utils.temp_dir`.
Whenever you write a test that needs to access files from or write files to the file system this should be done
in a temporary directory specific for the test, otherwise the files will be written all over the project
folders. This context manager can be easily used to create the temporary directory, optionally change directory
to that one, and changing back to the initial directory when leaving the context manager, optionally
deleting the temporary directory..

This function can be used in connection with the ``delete_tmp_dir`` fixture, that should be passed
to its ``delete`` argument. The ``delete_tmp_dir``, a boolean value ``True`` by default, can be switched
to ``False`` calling pytest with the ``--keep-tmpdir`` option. When this happens the path to the
temporary folder will be printed (with the ``print`` function) by ``temp_dir``. This means that you will
have this in the output for the failing tests and this will give you the way to inspect the files that were
used and produced by the test and better understand why a test has failed.

``delete_tmp_dir``, as well as the other general fixtures, is implemented in the pytest standard ``conftest.py``
file in the root folder of the project, with the others being mainly options to ease the access to the
test files. These files are stored in the ``turbomoleio/testfiles`` folder. You can add there the files that
might be needed for your unit and integration tests. Modifying existing files is possible, even though discouraged.
If it cannot be avoided (e.g. for the update to a new TURBOMOLE version) you should at least
check that all the tests relying on the files that you plan to modify will keep working as expected.
