.. pyird documentation master file, created by
   sphinx-quickstart on Mon Jan 11 14:38:51 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: ./figures/pyird_logo.png
   :scale: 30%

PyIRD
==================================

`PyIRD` is a Python-based pipeline for reducing spectroscopic data obtained with `IRD <https://ird.mtk.nao.ac.jp/IRDpub/index_tmp.html>`_ and `REACH <https://secondearths.sakura.ne.jp/reach/>`_ on the Subaru Telescope. 
It is designed to process raw images into one-dimensional spectra in a semi-automatic manner. 

.. toctree:: 
   :maxdepth: 2
   :caption: Introduction:

   getting_started.rst
   input_data.rst

.. toctree::
   :maxdepth: 2
   :caption: Tutorial:

   tutorials/IRD_stream.rst
   tutorials/plot_spec.rst
   tutorials/other_outputs.rst

.. toctree::
   :maxdepth: 2
   :caption: Userguide:

   userguide/backgrounds.rst
   userguide/aptrace_flat_star.rst
   userguide/other_settings.rst
   developers.rst


.. toctree::
   :maxdepth: 1
   :caption: API:

   pyird/pyird.rst




License & Attribution
---------------------

Copyright 2021-2024, 
Contributors:

- `Yui Kasagi <https://yuikasagi.com>`_ (maintainer), 
- Hajime Kawahara (co-maintainer), 
- Ziying Gu, 
- Teruyuki Hirano, 
- Takayuki Kotani, 
- Masayuki Kuzuhara, 
- Kento Masuda

PYIRD is free software made available under the MIT License. For details
see the `LICENSE <https://github.com/prvjapan/pyird/blob/develop/LICENSE>`_.
   
Community Guidelines
---------------------

For developers, please read `CONTRIBUTING.md <https://github.com/prvjapan/pyird/blob/develop/CONTRIBUTING.md>`_ for details how to contribute to the project.
If you encounter a bug or have suggestions, feel free to open an `Issue <https://github.com/prvjapan/pyird/issues>`_.