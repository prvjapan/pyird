Test code for developers
==========================

PyIRD includes test code to ensure the maintainability of the package, using the `pytest <https://docs.pytest.org/en/stable/>`_ testing framework.

How to run tests
^^^^^^^^^^^^^^^^^
To run the tests, make sure that `pytest` is installed in your Python environment.
Then, run the following commands from the root directory of the repository:

.. code-block:: bash

   pytest tests/unittests
   pytest tests/integration

Or, 

.. code-block:: bash

   cd tests/unittests
   pytest


Tests in ``tests/unittests`` are automatically run by GitHub Actions on pull requests to develop or master.

Unit tests
""""""""""""

``tests/unittests``: Tests for individual functions and classes. 
Developers encauraged to create new tests when adding new features.

Integration tests
"""""""""""""""""""

``tests/integration``: Tests for the behavior of multiple functions and classes working together.