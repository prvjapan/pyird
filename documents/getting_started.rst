Getting Started
==================================

Installation
---------------------

Install PyIRD using `pip`:

.. code-block:: 

   pip install pyird

or by cloning the GitHub repository:

.. code-block:: 

   git clone https://github.com/prvjapan/pyird.git
   cd pyird
   pip install .

Verify Installation
^^^^^^^^^^^^^^^^^^^^

To verify that PyIRD is correctly installed, run:

.. code-block::

   python -m pyird --version

or simply import it in Python:

.. code-block:: python

   import pyird
   print(pyird.__version__)

Running the Tutorial
---------------------

You can follow the tutorial provided in :doc:`./tutorials/IRD_stream` to learn how to
process and analyze IRD data step-by-step.

.. note::
   The tutorial script ``IRD_stream.py`` is located in the ``examples/python/`` folder.

Before running the tutorial, make sure that the sample data
(from Zenodo) is extracted in the correct directory.

Sample Data
-------------

The raw data used in the tutorial (IRD_stream.py) can be downloaded from the `Zenodo repository <https://zenodo.org/records/14614004>`_.
The total file size after extraction is approximately 7 GB.

Input Data Overview
^^^^^^^^^^^^^^^^^^^^^^

For details about supported input data formats, see :doc:`./input_data`.

Next Steps
^^^^^^^^^^^^^^

- :doc:`./input_data` — Explanation of input data structure
- :doc:`./tutorials/IRD_stream` — Step-by-step usage guide