Names and Terms of IRD/REACH data
==================================

Terms Related to Instruments
------------------------------

- ``band``: Y/J-band detector or H-band detector
    - The origin of the detector images below is at the lower right.
    - If the images in your data set are rotated in different orientation, use the ``rotate`` and/or ``inverse`` instance in ``Stream2D``. 

.. list-table:: IRD Detectors
  :widths: 15 10
  :header-rows: 1

  * - Detector Image
    - Note
  * - .. image:: ../figures/IRDA00041704_raw.png
    - - Y/J-band
      - 51 orders x 2 fibers = 102 apertures
  * - .. image:: ../figures/IRDA00041705_raw.png
    - - H-band
      - 21 orders x 2 fibers = 42 apertures

- ``mmf``: Multi Mode Fiber used for IRD
    - mmf1 (or comb-fiber): A fiber to inject the light from Laser Frequency Comb (LFC) for simultaneous wavelength reference.
    - mmf2 (or star-fiber): A fiber to inject the light of a star.

- ``smf``: Single Mode Fiber used for REACH

Names of Calibration Dataset
-----------------------------

- ``thar (Th-Ar)``: Thrium-Argon hollow cathode lamp for wavelength reference

- ``FLAT_COMB, FLAT_STAR``: 'FLAT' means the light from a white lamp for flat fielding. '_COMB' and '_STAR' denote the injection to the comb-fiber (mmf1) and star-fiber (mmf2), respectively.
    - In the FLAT_COMB image, flat light injected on the comb fiber is also leaking into the star fiber, and two flat lights, one strong and one weak, are seen at each order.
    - In the pipeline, by utilizing this feature, all the apertures are identified in the FLAT-COMB, and then selected the ones to be used.
    - The advantage of this method is that there is no need to change frame IDs when switching fiber to analyze.

- ``ThAr-ThAr``: Calibration image with ThAr light inject to both fibers.

- ``COMB-COMB``: Calibration image with LFC light inject to both fibers.