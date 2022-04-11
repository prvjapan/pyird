import numpy as np
import pkg_resources
from pyird.utils import irdstream
import pathlib
from pyird.image.bias import bias_subtract_image
from pyird.image.hotpix import identify_hotpix
from pyird.image.trace_function import trace_legendre
from pyird.image.aptrace import aptrace
import astropy.io.fits as pyf

# aperture extraction
datadir = pathlib.Path('/home/kawahara/pyird/data/flat/')
anadir = pathlib.Path('/home/kawahara/pyird/data/flat/')
flat_smf66=irdstream.Stream2D("flat_smf66",datadir,anadir)
flat_smf66.fitsid=list(range(47224,47246,2)) #no light in the speckle fiber
flat_smf66.fitsid_increment()
flatmedian=flat_smf66.immedian()
cutrow = 1000  ## 501 ~ 1549
nap = 102 ## 42 for H band, 102 for YJ band
if nap==42:
    flatmedian = flatmedian[::-1,::-1]
y0, xmin, xmax, coeff = aptrace(flatmedian,cutrow,nap)

# hotpixel mask
datadir = pathlib.Path('/home/kawahara/pyird/data/dark/')
anadir = pathlib.Path('/home/kawahara/pyird/data/dark/')
dark = irdstream.Stream2D('targets', datadir, anadir)
dark.fitsid = [41018]
for data in dark.rawpath:
    im = pyf.open(str(data))[0].data
im_subbias = bias_subtract_image(im)
hotpix_mask, obj = identify_hotpix(im_subbias)

# Load data
datadir = pathlib.Path('/home/kawahara/pyird/data/samples/REACH/')
anadir = pathlib.Path('/home/kawahara/pyird/data/samples/REACH/')
target = irdstream.Stream2D(
    'targets', datadir, anadir, fitsid=[47078, 47080, 47082])
target.info = True  # show detailed info

# clean pattern
pathC = (pkg_resources.resource_filename('pyird', 'data/samples/aprefC'))
path_c = (pkg_resources.resource_filename('pyird', 'data/samples/apref_c'))
target.clean_pattern(extin='', extout='_cp', trace_path_list=[
                     pathC, path_c], hotpix_mask=hotpix_mask)

# flatten
path_trace_flatten = (pkg_resources.resource_filename(
    'pyird', 'data/samples/aprefB'))
target.flatten(path_trace_flatten)

# load ThAr raw image
datadir = pathlib.Path('/home/kawahara/pyird/data/samples/REACH/')
anadir = pathlib.Path('/home/kawahara/pyird/data/samples/REACH/')

#wavelength calibration
pathB = (pkg_resources.resource_filename('pyird', 'data/samples/aprefB'))
thar=irdstream.Stream2D("thar",datadir,anadir,rawtag="IRDBD000",fitsid=list(range(15480,15530)))
thar.clean_pattern(extin='', extout='_cp', trace_path_list=[pathC, path_c], hotpix_mask=hotpix_mask)
wavsol, data=thar.calibrate_wavlength(pathB)
