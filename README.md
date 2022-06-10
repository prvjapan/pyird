# PYIRD

Free from messy directory and file management of IRD analysis.
Currently under heavily construction. Use develop branch for dev.


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

Aperture
------------------------------

Mask: n=42 (H) apertures extracted from FLAT for aperture mask.

Star: `aperture.TraceAperture.mmf2()` sets n=21 (H) apertures of star fiber.
