import numpy as np
#import pkg_resources
from pyird.utils import irdstream
import pathlib
from pyird.image.bias import bias_subtract_image
from pyird.image.hotpix import identify_hotpix_sigclip
import astropy.io.fits as pyf

# path
basedir = pathlib.Path('~/pyird/data/20210317/').expanduser()

### FOR CALIBRATION ###
# aperture extraction
datadir = basedir/'flat/'
anadir = basedir/'flat/'
flat_comb=irdstream.Stream2D("flat_comb",datadir,anadir)
flat_comb.fitsid=list(range(41704,41804,2)) ##FLAT_COMB
################################
### SELECT H band or YJ band ###
################################
flat_comb.band='h' #'h' or 'y'
print(flat_comb.band,' band')
if flat_comb.band=='h':
    flat_comb.fitsid_increment() # when you use H-band
    trace_mmf=flat_comb.aptrace(cutrow = 800,nap=42) #TraceAperture instance
elif flat_comb.band=='y':
    trace_mmf=flat_comb.aptrace(cutrow = 1000,nap=102) #TraceAperture instance
trace_mask = trace_mmf.mask()

import matplotlib.pyplot as plt
plt.imshow(trace_mmf.mask()) #apeture mask plot
plt.show()

# hotpixel mask: See pyird/io/read_hotpix.py for reading fixed mask (Optional)
datadir = basedir/'dark/'
anadir = basedir/'dark/'
dark = irdstream.Stream2D('dark', datadir, anadir,fitsid=[41504]) # Multiple file is ok
if flat_comb.band=='h':
    dark.fitsid_increment() # when you use H-band
median_image = dark.immedian()
im_subbias = bias_subtract_image(median_image)
hotpix_mask = identify_hotpix_sigclip(im_subbias)

###########################
### SELECT mmf2 or mmf1 ###
###########################
trace_mmf.mmf2() #mmf2 (star fiber)
#trace_mmf.mmf1() #mmf1 (comb fiber)

# load ThAr raw image
datadir = basedir/'thar'
anadir = basedir/'thar'
if flat_comb.band=='h':
    rawtag='IRDAD000'
elif flat_comb.band=='y':
    rawtag='IRDBD000'

#wavelength calibration
thar=irdstream.Stream2D("thar",datadir,anadir,rawtag=rawtag,fitsid=list(range(14632,14732))) #ThAr-ThAr
thar.trace = trace_mmf
thar.clean_pattern(trace_mask=trace_mask,extin='', extout='_cp', hotpix_mask=hotpix_mask)
thar.calibrate_wavelength()


### FLAT ###
datadir = basedir/'flat/'
anadir = basedir/'flat/'
flat_star=irdstream.Stream2D("flat_star",datadir,anadir)
flat_star.fitsid=list(range(41804,41904,2)) ##FLAT_STAR
flat_star.trace = trace_mmf
flat_star.band='h' #'h' or 'y'
if flat_star.band == 'h':
    flat_star.fitsid_increment() # when you use H-band
    flat_star.clean_pattern(trace_mask=trace_mask,extin='', extout='_cp', hotpix_mask=hotpix_mask)
flat_star.imcomb = True # median combine
flat_star.flatten(hotpix_mask=hotpix_mask)
df_flatn = flat_star.apnormalize()

### TARGET ###
# Load data
datadir = basedir/'target/'
anadir = basedir/'target/'
target = irdstream.Stream2D(
    'targets', datadir, anadir, fitsid=[41510])
if flat_comb.band=='h':
    target.fitsid_increment() # when you use H-band
target.info = True  # show detailed info
target.trace = trace_mmf
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
target.normalize1D(master_path=flat_star.anadir,skipLFC=True)#,flatid=flat_comb.streamid)

"""
### FOR RV MEASUREMENTS ###
### mmfmmf (test) ###
datadir = basedir/'mmfmmf/'
anadir = basedir/'reduc_rv/'
mmfmmf=irdstream.Stream2D("mmfmmf",datadir,anadir,fitsid=list(range(14734,14832)),rawtag=rawtag)
mmfmmf.trace = trace_mmf
mmfmmf.clean_pattern(trace_mask=trace_mask,extin='', extout='_cp', hotpix_mask=hotpix_mask)
mmfmmf.imcomb = True
mmfmmf.flatten()
mmfmmf.dispcor(master_path=thar.anadir)
mmfmmf.normalize1D(master_path=flat.anadir)

### hotpix (test) ###
from pyird.image.hotpix import hotpix_fits_to_dat
dark.trace = trace_mmf
dark.imcomb = True
dark.flatten(hotpix_mask=hotpix_mask)
file = dark.anadir/('%s_%s_%s.fits'%(dark.streamid,dark.band,dark.trace.mmf))
save_path = dark.anadir/('%s_%s_%s.dat'%(dark.streamid,dark.band,dark.trace.mmf))
hotpix_fits_to_dat(file,save_path)
"""
