PYIRD
===========

Free from messy directory, file management, wavelength calibration, IRAF stuffs, and so on, when analyzing the Subaru/IRD data.


Examples
------------------

Creating 1D spectra from raw data.
- pyird/examples/python/IRD_stream.py -- IRD data; read [the docs](https://secondearths.sakura.ne.jp/pyird/tutorials/IRD_stream.html).
- pyird/examples/python/REACH_stream.py -- REACH data


Install
------------------

```
pip install pyird
```

or

```
python setup.py install
```


Core Pipelines
------------------

- pyird.utils -- setting us free from file name managements
- pyird.image -- reducing systematics in the image level and extracting raw spectra
- pyird.spec -- assigning wavelength reference (wavref) to pixel coordinates in the raw spectra


Classes
------------------

- fitsset.FitsSet --  sets of fits files
- irdstream.Stream2D -- 2D fits stream from raw images
- aperture.TraceAperture -- aperture class for trace
- image_widget.image_1Dand2D -- interactive figures of 1D spectrum and 2D detector image

Aperture
------------------------------

Mask: n=42 (H), 102 (YJ) apertures extracted from FLAT for aperture mask.

Star: `aperture.TraceAperture.mmf2()` sets n=21 (H) apertures of star fiber.

Comb: `aperture.TraceAperture.mmf1()` sets n=51 (YJ) apertures of comb fiber.
