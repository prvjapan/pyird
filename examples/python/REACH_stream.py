import numpy as np
import pkg_resources
from pyird.utils import irdstream
import pathlib
from pyird.image.bias import bias_subtract_image
from pyird.image.hotpix import identify_hotpix
import astropy.io.fits as pyf

# base directory
basedir = pathlib.Path('~/pyird/data/REACH/').expanduser()
print (basedir)

# aperture extraction
datadir = basedir/pathlib.Path('flat')
anadir = basedir/pathlib.Path('flat')
flat_smf66=irdstream.Stream2D("flat_smf66",datadir,anadir)
flat_smf66.fitsid=list(range(47224,47246,2)) #no light in the speckle fiber
flat_smf66.fitsid_increment() # when you use H-band
trace_smf66=flat_smf66.aptrace(cutrow = 1000,nap=42) #TraceAperture instance

import matplotlib.pyplot as plt
plt.imshow(trace_smf66.mask()) #apeture mask plot
plt.show()

# hotpixel mask
datadir = basedir/pathlib.Path('dark')
anadir = basedir/pathlib.Path('dark')
dark = irdstream.Stream2D('targets', datadir, anadir)
dark.fitsid = [41018]
for data in dark.rawpath:
    im = pyf.open(str(data))[0].data
im_subbias = bias_subtract_image(im)
hotpix_mask, obj = identify_hotpix(im_subbias)

# Load data
datadir = basedir/pathlib.Path('samples')#/pathlib.Path('REACH')
anadir = basedir/pathlib.Path('samples')#/pathlib.Path('REACH')
target = irdstream.Stream2D(
    'targets', datadir, anadir, fitsid=[47078, 47080, 47082]) # YJ band?
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
datadir = basedir/pathlib.Path('samples')#/pathlib.Path('REACH')
anadir = basedir/pathlib.Path('samples')#/pathlib.Path('REACH')

#wavelength calibration
pathB = (pkg_resources.resource_filename('pyird', 'data/samples/aprefB'))
thar=irdstream.Stream2D("thar",datadir,anadir,rawtag="IRDBD000",fitsid=list(range(15480,15530)))
thar.clean_pattern(extin='', extout='_cp', trace_path_list=[pathC, path_c], hotpix_mask=hotpix_mask)
thar.calibrate_wavlength(pathB)
