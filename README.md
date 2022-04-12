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

The followings are out dated info.

Aperture
------------------------------

For n=51 or 52 (YJ) and n=21 (H)
Mask 104 (YJ) 


Help tool for ECIDENTIFY
--------------------------

- calref.py

s option: The number of the orders varies depending on environment. s option gives an offset that defines the order numbering. Try it if you couldn't find good match. In particular, for SMF/YJ, try s = 1.

xorder 4
yorder 3
