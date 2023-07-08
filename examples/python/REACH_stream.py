import numpy as np
import pkg_resources
from pyird.utils import irdstream
import pathlib
from pyird.image.bias import bias_subtract_image
from pyird.image.hotpix import identify_hotpix_sigclip
import astropy.io.fits as pyf

# path
basedir = pathlib.Path('~/pyird/data/20211110_REACH/').expanduser()

### For REACH ###
inst = 'REACH'

### FOR CALIBRATION ###
# aperture extraction
datadir = basedir/'flat/'
anadir = basedir/'flat/'
flat_star=irdstream.Stream2D("flat",datadir,anadir,inst=inst)
flat_star.fitsid=list(range(53235,53334,2)) #no light in the speckle fiber
################################
### SELECT H band or YJ band ###
################################
flat_star.band='h' #'h' or 'y'
print(flat_star.band,' band')
if flat_star.band=='h' and flat_star.fitsid[0]%2==0:
    flat_star.fitsid_increment() # when you use H-band
    trace_smf=flat_star.aptrace(cutrow = 800,nap=42) #TraceAperture instance
elif flat_star.band=='y':
    trace_smf=flat_star.aptrace(cutrow = 1000,nap=102) #TraceAperture instance
trace_mask = trace_smf.mask()

import matplotlib.pyplot as plt
plt.imshow(trace_smf.mask()) #apeture mask plot
plt.show()

# hotpixel mask: See pyird/io/read_hotpix.py for reading fixed mask (Optional)
datadir = basedir/'dark/'
anadir = basedir/'dark/'
dark = irdstream.Stream2D('dark', datadir, anadir,fitsid=[47269],inst=inst)
if flat_star.band=='h' and dark.fitsid[0]%2==0:
    dark.fitsid_increment() # when you use H-band
median_image = dark.immedian()
for data in dark.rawpath:
    im = pyf.open(str(data))[0].data
im_subbias = bias_subtract_image(im)
hotpix_mask = identify_hotpix_sigclip(im_subbias)

###########################
### SELECT mmf2 or mmf1 ###
###########################
trace_smf.mmf2() #mmf2 (star fiber)
#trace_smf.mmf1() #mmf1 (comb fiber)

# load ThAr raw image
datadir = basedir/'thar'
anadir = basedir/'reduc_fixblaze'
if flat_star.band=='h':
    rawtag='IRDAD000'
elif flat_star.band=='y':
    rawtag='IRDBD000'

#wavelength calibration
thar=irdstream.Stream2D("thar",datadir,anadir,fitsid=list(range(53347,53410,2)),inst=inst)
thar.trace = trace_smf
thar.clean_pattern(trace_mask=trace_mask,extin='', extout='_cp', hotpix_mask=hotpix_mask)
thar.calibrate_wavelength()

### FLAT ###
flat_star.trace = trace_smf
flat_star.clean_pattern(trace_mask=trace_mask,extin='', extout='_cp', hotpix_mask=hotpix_mask)
flat_star.imcomb = True # median combine
flat_star.flatten(hotpix_mask=hotpix_mask)
df_flatn = flat_star.apnormalize()

### TARGET ###
# Load data
datadir = basedir/'target/'
anadir = basedir/'target/'
target = irdstream.Stream2D(
    'targets', datadir, anadir, fitsid=[53205],inst=inst)
if flat_star.band=='h' and target.fitsid[0]%2==0:
    target.fitsid_increment() # when you use H-band
target.info = True  # show detailed info
target.trace = trace_smf
# clean pattern
target.clean_pattern(trace_mask=trace_mask,extin='', extout='_cp', hotpix_mask=hotpix_mask)
# flatten
target.apext_flatfield(df_flatn,hotpix_mask=hotpix_mask)
# assign reference spectra & resample
target.dispcor(master_path=thar.anadir,extin='_flnhp')#_hp

# blaze function
flat_star.apext_flatfield(df_flatn,hotpix_mask=hotpix_mask)
flat_star.dispcor(master_path=thar.anadir)

# combine & normalize
target.normalize1D(master_path=flat_star.anadir,skipLFC=True)

# save noramlized flat for fringe removal
flat_star.dispcor(master_path=thar.anadir,blaze=False)
flat_star.normalize1D(master_path=flat_star.anadir,skipLFC=True)
