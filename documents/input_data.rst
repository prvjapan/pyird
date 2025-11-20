Input Data
============

.. image:: ./figures/reduc_flowchart.png

This page describes the input data required for the data reduction process (shown in the green boxes in the image above).  
Please also refer to :doc:`./userguide/terms` for the naming conventions of observation data and detailed explanations of calibration datasets.

- Flat_comb  
    - Frame IDs in the sample data: IRDA00041704.fits - IRDA00041803.fits  
    - **Required** for aperture tracing (for **both** fibers, by default)
    - **Required** for flat fielding and normalizing the spectrum *of the comb fiber*
    
- Flat_star  
    - Frame IDs in the sample data: IRDA00041804.fits - IRDA00041903.fits  
    - **Optionally** used for aperture tracing *for the star fiber*
    - **Required** for flat fielding and  normalizing the spectrum *of the star fiber*


.. tip::
    See also :ref:`terms-calibration` for flat data used for aperture tracing.

- Dark  
    - Frame IDs in the sample data: IRDA00041504.fits  
    - **Optionally** used for identifying hot pixels

.. tip::
    See also :ref:`read-hotpixmask` for another option for handling hot pixels.

- ThAr  
    - Frame IDs in the sample data: IRD[B/A]D00014632.fits - IRD[B/A]D00014732.fits  
    - **Required** for wavelength calibration

- Target (Science)  
    - Frame IDs in the sample data: IRDA00041510.fits, IRDA00041511.fits  
    - Main observation frames

.. note:: 
    Supported instruments in PyIRD v1.1.0:

        - ✅ IRD  
        - ✅ REACH  
        - ☑️ IRCS (available on the develop branch)
