import numpy as np
import pkg_resources
from pyird.utils import irdstream
import pathlib
from pyird.image.bias import bias_subtract_image
from pyird.image.hotpix import identify_hotpix
import astropy.io.fits as pyf

# path
basedir = pathlib.Path('~/pyird/data/20211110_REACH/').expanduser()

### For REACH ###
inst = 'REACH'

### FOR CALIBRATION ###
# aperture extraction
datadir = basedir/'flat/'
anadir = basedir/'flat/'
flat=irdstream.Stream2D("flat",datadir,anadir,inst=inst)
flat.fitsid=list(range(53235,53334,2)) #no light in the speckle fiber
################################
### SELECT H band or YJ band ###
################################
flat.band='h' #'h' or 'y'
print(flat.band,' band')
if flat.band=='h':
    #flat.fitsid_increment() # when you use H-band
    trace_smf=flat.aptrace(cutrow = 800,nap=42) #TraceAperture instance
elif flat.band=='y':
    trace_smf=flat.aptrace(cutrow = 1000,nap=102) #TraceAperture instance

import matplotlib.pyplot as plt
plt.imshow(trace_smf.mask()) #apeture mask plot
plt.show()

# hotpixel mask
datadir = basedir/'dark/'
anadir = basedir/'dark/'
dark = irdstream.Stream2D('dark', datadir, anadir,fitsid=[47269],inst=inst)
#if flat.band=='h':
#    dark.fitsid_increment() # when you use H-band
for data in dark.rawpath:
    im = pyf.open(str(data))[0].data
im_subbias = bias_subtract_image(im)
hotpix_mask, obj = identify_hotpix(im_subbias)

###########################
### SELECT mmf2 or mmf1 ###
###########################
#trace_smf.mmf2() #mmf2 (star fiber)
trace_smf.mmf1() #mmf1 (comb fiber)

# load ThAr raw image
datadir = basedir/'thar'
anadir = basedir/'thar'
if flat.band=='h':
    rawtag='IRDAD000'
elif flat.band=='y':
    rawtag='IRDBD000'

#wavelength calibration
thar=irdstream.Stream2D("thar",datadir,anadir,fitsid=list(range(53411,53460,2)),inst=inst) #rawtag=rawtag,
thar.trace = trace_smf
thar.clean_pattern(extin='', extout='_cp', hotpix_mask=hotpix_mask)
thar.calibrate_wavelength()

### TARGET ###
# Load data
datadir = basedir/'target/'
anadir = basedir/'target/'
target = irdstream.Stream2D(
    'targets', datadir, anadir, fitsid=[53205],inst=inst)
#if flat.band=='h':
#    target.fitsid_increment() # when you use H-band
target.info = True  # show detailed info
target.trace = trace_smf
# clean pattern
target.clean_pattern(extin='', extout='_cp', hotpix_mask=hotpix_mask)
# flatten
target.flatten()
# assign reference spectra & resample
target.dispcor(master_path=thar.anadir)

### FLAT (for blaze function) ###
flat.trace = trace_smf
if flat.band == 'h':
    flat.clean_pattern(extin='', extout='_cp', hotpix_mask=hotpix_mask)
flat.imcomb = True # median combine
flat.flatten()
flat.dispcor(master_path=thar.anadir)

# combine & normalize
target.normalize1D(flatid=flat.streamid,master_path=flat.anadir)
