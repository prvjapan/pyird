Useful Options
==============

- :ref:`handle-fitsid`
- :ref:`read-hotpixmask`
- :ref:`output-format`
- :ref:`reach-reduction`

.. _handle-fitsid:

Increment/Decrement FITS IDs
----------------------------

Frame numbers in IRD data are assigned as even numbers for the Y/J band
and odd numbers for the H band.

Manual Increment/Decrement
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    import pathlib
    from pyird.utils import irdstream
    
    basedir = pathlib.Path('~/pyird/data/20210317/').expanduser()
    datadir = basedir/'target/'
    anadir = basedir/'reduc/'
    
    id_demo = irdstream.Stream2D("id_demo", datadir, anadir, fitsid=[41510])
    
    id_demo.fitsid_increment()


.. parsed-literal::

    fitsid: [41510]
    fitsid incremented:  [41511]


.. code:: ipython3

    id_demo.fitsid_decrement()


.. parsed-literal::

    fitsid decremented:  [41510]


Automatic Increment for H-band Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    from pyird.utils import irdstream
    
    id_demo = irdstream.Stream2D("id_demo", datadir, anadir, fitsid=[41510], band="h")


.. parsed-literal::

    fitsid: [41510]
    fitsid incremented:  [41511]


.. _read-hotpixmask:

Read Hotpixel Mask
------------------

``PyIRD`` includes hotpixel masks created from dark data obtained in
Oct. 2022. You can use these masks if dark data from your observation
date is unavailable.

.. code:: ipython3

    import importlib
    from pyird.io.read_hotpix import read_hotpix
    
    band = "h"
    
    if band=='h':
        path=importlib.resources.files('pyird').joinpath('data/hotpix_mask_h_202210_180s.fits')
    elif band=='y':
        path=importlib.resources.files('pyird').joinpath('data/hotpix_mask_y_202210_180s.fits')
    hotpix_mask=read_hotpix(path)

.. _output-format:

Change Output Format
--------------------

The output format of ``pandas.DataFrame()`` is determined by
``tocsvargs`` in the class ``Stream2D()``. By default,
``self.tocsvargs = {"header": False, "index": False, "sep": " "}``, but
you can modify these settings as needed!

.. code:: ipython3

    from pyird.utils import irdstream
    
    output_demo = irdstream.Stream2D("id_demo", datadir, anadir, fitsid=[41510])
    print("default: ", output_demo.tocsvargs)
    
    output_demo.tocsvargs = {"header": False, "index": True, "sep": ","}
    print("modified: ", output_demo.tocsvargs)


.. parsed-literal::

    fitsid: [41510]
    default:  {'header': False, 'index': False, 'sep': ' '}
    modified:  {'header': False, 'index': True, 'sep': ','}


.. _reach-reduction:

REACH Data Reduction
--------------------

For REACH data, set ``inst = REACH`` in ``Stream2D()``.

This option changes default aperture width to 5 pixels (ranging from -2
to 3).

.. code:: ipython3

    from pyird.utils import irdstream
    
    reach_demo = irdstream.Stream2D("reach_demo", datadir, anadir, fitsid=[41510], inst="REACH")


.. parsed-literal::

    fitsid: [41510]

