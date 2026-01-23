Names and Terms of IRD/REACH data
==================================

For information about each instrument, visit
`IRD official page <https://ird.mtk.nao.ac.jp/IRDpub/index_tmp.html>`_
or `REACH official page <https://secondearths.sakura.ne.jp/reach/>`_.

Terms Related to IRD/REACH
--------------------------------

- ``band``: Y/J-band detector or H-band detector
    - The origin of the detector images below is at the lower left.
    - If the images in your data set are rotated in different orientation, use the ``rotate`` and/or ``inverse`` instance in ``Stream2D``. 

.. list-table:: IRD Detectors
  :widths: 15 10
  :header-rows: 1

  * - Detector Image
    - Note
  * - .. image:: ../figures/IRDA00041704_raw.png
          :width: 30%  
    - - Y/J-band
      - 51 orders x 2 fibers = 102 apertures
  * - .. image:: ../figures/IRDA00041705_raw.png
          :width: 30%
    - - H-band
      - 21 orders x 2 fibers = 42 apertures

- ``mmf``: Multi Mode Fiber used for IRD
    - mmf1 (or comb-fiber): A fiber to inject the light from Laser Frequency Comb (LFC) for simultaneous wavelength reference.
    - mmf2 (or star-fiber): A fiber to inject the light of a star.

- ``smf``: Single Mode Fiber used for REACH

.. _terms-filename:

Naming Conventions for Observed Data with IRD/REACH
--------------------------------------------------------------

Each exposure taken with IRD and REACH produces a pair of FITS files: one for the Y/J-band detector and the other for the H-band detector. 

The prefix of each file name varies depending on the observation mode:

- ``IRDA``: Data obtained via Gen2, the observation control system of the Subaru Telescope.
    - The *even-numbered* FITS file (e.g., ``IRDA00041510.fits``) contains the Y/J-band data.
    - The *odd-numbered* FITS file (e.g., ``IRDA00041511.fits``) contains the H-band data.

- ``IRDBD`` or ``IRDAD``: Data obtained directly from the detector PC, mainly used for calibration purposes.
    - ``IRDBDxxxxxxxx.fits``: Y/J-band data
    - ``IRDADxxxxxxxx.fits``: H-band data 

.. _terms-calibration:

Names of Calibration Dataset
-----------------------------

- ``thar (Th-Ar)``: Thrium-Argon hollow cathode lamp for wavelength reference

- ``FLAT_COMB, FLAT_STAR``: 'FLAT' means the light from a white lamp for flat fielding. '_COMB' and '_STAR' denote the injection to the comb-fiber (mmf1) and star-fiber (mmf2), respectively.
    - In the FLAT_COMB image, flat light injected on the comb fiber is also leaking into the star fiber, and two flat lights, one strong and one weak, are seen at each order.
    - In the pipeline, by utilizing this feature, all the apertures are identified in the FLAT-COMB, and then selected the ones to be used.
    - The advantage of this method is that there is no need to change frame IDs when switching fiber to analyze.

- ``ThAr-ThAr``: Calibration image with ThAr light inject to both fibers.

- ``COMB-COMB``: Calibration image with LFC light inject to both fibers.