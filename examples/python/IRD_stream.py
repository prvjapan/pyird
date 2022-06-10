import numpy as np
#import pkg_resources
from pyird.utils import irdstream
import pathlib
from pyird.image.bias import bias_subtract_image
from pyird.image.hotpix import identify_hotpix
import astropy.io.fits as pyf

# path
basedir = pathlib.Path('~/pyird/data/20210317/').expanduser()

### FOR CALIBRATION ###
# aperture extraction
datadir = basedir/'flat/'
anadir = basedir/'flat/'
flat_mmf=irdstream.Stream2D("flat_mmf",datadir,anadir)
flat_mmf.fitsid=list(range(41704,41804,2)) #no light in the speckle fiber
flat_mmf.fitsid_increment() # when you use H-band
trace_mmf=flat_mmf.aptrace(cutrow = 1000,nap=42) #TraceAperture instance

import matplotlib.pyplot as plt
plt.imshow(trace_mmf.mask()) #apeture mask plot
plt.show()

# hotpixel mask
datadir = basedir/'dark/'
anadir = basedir/'dark/'
dark = irdstream.Stream2D('dark', datadir, anadir)
dark.fitsid = [41505]
for data in dark.rawpath:
    im = pyf.open(str(data))[0].data
im_subbias = bias_subtract_image(im)
hotpix_mask, obj = identify_hotpix(im_subbias)

# load ThAr raw image
datadir = basedir/'thar'
anadir = basedir/'thar'

#wavelength calibration
thar=irdstream.Stream2D("thar",datadir,anadir,rawtag="IRDAD000",fitsid=list(range(14632,14732)))
thar.trace = trace_mmf
trace_mmf.mmf2() # choose mmf2 (star fiber)
thar.clean_pattern(extin='', extout='_cp', hotpix_mask=hotpix_mask)
thar.calibrate_wavlength()

### TARGET ###
# Load data
datadir = basedir/'target/'
anadir = basedir/'target/'
target = irdstream.Stream2D(
    'targets', datadir, anadir, fitsid=[41511])
target.info = True  # show detailed info
target.trace = trace_mmf
# clean pattern
target.clean_pattern(extin='', extout='_cp', hotpix_mask=hotpix_mask)
# flatten
target.flatten()
# assign reference spectra & resample
target.dispcor()

### FLAT (for blaze function) ###
flat_mmf.trace = trace_mmf
flat_mmf.clean_pattern(extin='', extout='_cp', hotpix_mask=hotpix_mask)
flat_mmf.imcomb = True # median combine
flat_mmf.flatten()
flat_mmf.dispcor()

# combine & normalize
target.normalize1D(flatid=flat_mmf.streamid)

"""### mmfmmf (test) ###
datadir = basedir/'mmfmmf/'
anadir = basedir/'mmfmmf/'
mmfmmf=irdstream.Stream2D("mmfmmf",datadir,anadir,rawtag="IRDAD000",fitsid=list(range(14733,14833)))
mmfmmf.trace = trace_mmf
mmfmmf.clean_pattern(extin='', extout='_cp', hotpix_mask=hotpix_mask)
mmfmmf.imcomb = True
mmfmmf.flatten()
mmfmmf.dispcor()
"""
