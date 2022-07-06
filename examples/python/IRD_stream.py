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
flat=irdstream.Stream2D("flat",datadir,anadir)
flat.fitsid=list(range(41704,41804,2))
################################
### SELECT H band or YJ band ###
################################
flat.band='h' #'h' or 'y'
print(flat.band,' band')
if flat.band=='h':
    flat.fitsid_increment() # when you use H-band
    trace_mmf=flat.aptrace(cutrow = 1500,nap=42) #TraceAperture instance
elif flat.band=='y':
    trace_mmf=flat.aptrace(cutrow = 1000,nap=102) #TraceAperture instance

import matplotlib.pyplot as plt
plt.imshow(trace_mmf.mask()) #apeture mask plot
plt.show()

# hotpixel mask
datadir = basedir/'dark/'
anadir = basedir/'dark/'
dark = irdstream.Stream2D('dark', datadir, anadir,fitsid=[41504])
if flat.band=='h':
    dark.fitsid_increment() # when you use H-band
for data in dark.rawpath:
    im = pyf.open(str(data))[0].data
im_subbias = bias_subtract_image(im)
hotpix_mask, obj = identify_hotpix(im_subbias)

###########################
### SELECT mmf2 or mmf1 ###
###########################
trace_mmf.mmf2() #mmf2 (star fiber)
#trace_mmf.mmf1() #mmf1 (comb fiber)

# load ThAr raw image
datadir = basedir/'thar'
anadir = basedir/'thar'
if flat.band=='h':
    rawtag='IRDAD000'
elif flat.band=='y':
    rawtag='IRDBD000'

#wavelength calibration
thar=irdstream.Stream2D("thar",datadir,anadir,rawtag=rawtag,fitsid=list(range(14632,14732)))
thar.trace = trace_mmf
thar.clean_pattern(extin='', extout='_cp', hotpix_mask=hotpix_mask)
thar.calibrate_wavelength()

### TARGET ###
# Load data
datadir = basedir/'target/'
anadir = basedir/'target/'
target = irdstream.Stream2D(
    'targets', datadir, anadir, fitsid=[41510])
if flat.band=='h':
    target.fitsid_increment() # when you use H-band
target.info = True  # show detailed info
target.trace = trace_mmf
# clean pattern
target.clean_pattern(extin='', extout='_cp', hotpix_mask=hotpix_mask)
# flatten
target.flatten()
# assign reference spectra & resample
target.dispcor(master_path=thar.anadir)

### FLAT (for blaze function) ###
flat.trace = trace_mmf
if flat.band == 'h':
    flat.clean_pattern(extin='', extout='_cp', hotpix_mask=hotpix_mask)
flat.imcomb = True # median combine
flat.flatten()
flat.dispcor(master_path=thar.anadir)

# combine & normalize
target.normalize1D(flatid=flat.streamid,master_path=flat.anadir)

"""
### FOR RV MEASUREMENTS ###
### mmfmmf (test) ###
datadir = basedir/'mmfmmf/'
anadir = basedir/'mmfmmf/'
mmfmmf=irdstream.Stream2D("mmfmmf",datadir,anadir,rawtag=rawtag,fitsid=list(range(15106,15206)))
mmfmmf.trace = trace_mmf
mmfmmf.clean_pattern(extin='', extout='_cp')#, hotpix_mask=hotpix_mask)
mmfmmf.imcomb = True
mmfmmf.flatten()
mmfmmf.dispcor(master_path=thar.anadir)
mmfmmf.normalize1D(master_path=flat.anadir)

### hotpix (test) ###
dark.trace = trace_mmf
dark.imcomb = True
dark.flatten(mask=hotpix_mask)
#dark.dispcor()
import pandas as pd
input = basedir/'reduc/dark_h_m2.fits' ## check!!!
hdu = pyf.open(input)[0]
spec_m2 = hdu.data
wspec = pd.DataFrame([],columns=['wav','order','flux'])
for i in range(len(spec_m2[0])):
    wav = range(1,2049)#reference[:,i]
    order = np.ones(len(wav))
    order[:] = i+1
    data_order = [wav,order,spec_m2[:,i]]
    df_order = pd.DataFrame(data_order,index=['wav','order','flux']).T
    wspec = pd.concat([wspec,df_order])
wspec = wspec.fillna(0)
save_path = basedir/'reduc/hotpix_h_m2.dat'
wspec.to_csv(save_path,header=False,index=False,sep=' ')

#dark.normalize1D()
"""
