import numpy as np
from pyird.utils import irdstream
import pathlib
from pyird.image.bias import bias_subtract_image
from pyird.image.hotpix import identify_hotpix_sigclip
import astropy.io.fits as pyf

#--------SETTINGS--------#
basedir = pathlib.Path('~/pyird/data/20210317/').expanduser()

band = 'h' #'h' or 'y'
mmf = 'mmf2' #'mmf1' (comb fiber) or 'mmf2' (star fiber)
skipLFC = False #if False, uncertainties are output. mmf1 of y band must be reduced first.

#--------FOR CALIBRATION--------#
## FLAT_COMB
# please change directry names and fits ids
datadir = basedir/'flat/'
anadir = basedir/'reduc/'
flat_comb=irdstream.Stream2D("flat_comb",datadir,anadir)
flat_comb.fitsid=list(range(41704,41804,2)) 
# aperture extraction
flat_comb.band=band
print(flat_comb.band,' band')
if band=='h' and flat_comb.fitsid[0]%2==0:
    flat_comb.fitsid_increment() 
    trace_mmf=flat_comb.aptrace(cutrow = 1200,nap=42) 
elif band=='y':
    trace_mmf=flat_comb.aptrace(cutrow = 1000,nap=102) 
trace_mask = trace_mmf.mask()

# apeture mask plot
import matplotlib.pyplot as plt
plt.imshow(trace_mmf.mask()) 
plt.show()

## HOTPIXEL MASK: 
# See pyird/io/read_hotpix.py for reading fixed mask (Optional)
## DARK
datadir = basedir/'dark/'
anadir = basedir/'reduc/'
dark = irdstream.Stream2D('dark', datadir, anadir,fitsid=[41504]) # Multiple file is ok
if band=='h' and dark.fitsid[0]%2==0:
    dark.fitsid_increment()
median_image = dark.immedian()
im_subbias = bias_subtract_image(median_image)
hotpix_mask = identify_hotpix_sigclip(im_subbias)

# reduce mmf1 or mmf2
if mmf=='mmf2':
    trace_mmf.choose_mmf2_aperture() #mmf2 (star fiber)
elif mmf=='mmf1':
    trace_mmf.choose_mmf1_aperture() #mmf1 (comb fiber)

## THAR (ThAr-ThAr)
datadir = basedir/'thar'
anadir = basedir/'reduc'
if band=='h':
    rawtag='IRDAD000'
elif band=='y':
    rawtag='IRDBD000'
thar=irdstream.Stream2D("thar",datadir,anadir,rawtag=rawtag,fitsid=list(range(14632,14732))) 

#wavelength calibration
thar.trace = trace_mmf
thar.clean_pattern(trace_mask=trace_mask,extin='', extout='_cp', hotpix_mask=hotpix_mask)
thar.calibrate_wavelength()

## FLAT
if mmf=='mmf2':
    ## FLAT_STAR
    datadir = basedir/'flat/'
    anadir = basedir/'reduc/'
    flat_star=irdstream.Stream2D("flat_star",datadir,anadir)
    flat_star.fitsid=list(range(41804,41904,2)) 
    flat_star.trace = trace_mmf
    flat_star.band=band 
    if band == 'h' and flat_star.fitsid[0]%2==0:
        flat_star.fitsid_increment() 
    flat_star.clean_pattern(trace_mask=trace_mask,extin='', extout='_cp', hotpix_mask=hotpix_mask)
    flat_star.imcomb = True # median combine
    flat_star.flatten(hotpix_mask=hotpix_mask)
    df_flatn = flat_star.apnormalize()
elif mmf=='mmf1':
    flat_comb.trace = trace_mmf
    flat_comb.clean_pattern(trace_mask=trace_mask,extin='', extout='_cp', hotpix_mask=hotpix_mask)
    flat_comb.imcomb = True # median combine
    flat_comb.flatten(hotpix_mask=hotpix_mask)
    df_flatn = flat_comb.apnormalize()

#--------FOR TARGET--------#
datadir = basedir/'target/'
anadir = basedir/'reduc/'
target = irdstream.Stream2D(
    'targets', datadir, anadir, fitsid=[41510])
if band=='h' and target.fitsid[0]%2==0:
    target.fitsid_increment() 
target.info = True  # show detailed info
target.trace = trace_mmf
# clean pattern
target.clean_pattern(trace_mask=trace_mask,extin='', extout='_cp', hotpix_mask=hotpix_mask)
# flatten
target.apext_flatfield(df_flatn,hotpix_mask=hotpix_mask)
# assign reference spectra & resample
target.dispcor(master_path=thar.anadir,extin='_flnhp')#_hp

if mmf=='mmf2':
    # blaze function
    flat_star.apext_flatfield(df_flatn,hotpix_mask=hotpix_mask)
    flat_star.dispcor(master_path=thar.anadir)

    # combine & normalize
    target.normalize1D(master_path=flat_star.anadir,skipLFC=skipLFC)#,flatid=flat_comb.streamid)
elif mmf=='mmf1':
    # blaze function
    flat_comb.apext_flatfield(df_flatn,hotpix_mask=hotpix_mask)
    flat_comb.dispcor(master_path=thar.anadir)

    # combine & normalize
    target.normalize1D(master_path=flat_comb.anadir,skipLFC=skipLFC)#,flatid=flat_comb.streamid)

"""
#--------FOR RV MEASUREMENTS--------#
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
