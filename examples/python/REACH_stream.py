import numpy as np
from pyird.utils import irdstream
import pathlib
from pyird.image.bias import bias_subtract_image
from pyird.image.hotpix import identify_hotpix_sigclip
import astropy.io.fits as pyf

#--------SETTINGS--------#
basedir = pathlib.Path('~/pyird/data/202106_REACH/').expanduser()

inst = 'REACH'
band = 'y' #'h' or 'y'
mmf = 'mmf2' #'mmf1' (comb fiber) or 'mmf2' (star fiber)
mask_hotpix = False #True
hotpix_mode = 'nan' # 'nan' or 'interp'
readout_noise_mode = "default" #'real' or 'default'

datadir_flat = basedir/'flat/'
datadir_dark = basedir/'dark/'
datadir_thar = basedir/'thar'
datadir_target = basedir/'target/'
anadir = basedir/'reduc/'

# last 5 digits of FITS file numbers: [start, end file number]
flat_star_id = [48024, 48123]  # flat image for star/target
dark_id = [47074] # dark image
thar_id = [15480, 15529] # ThAr image
target_id = [47998, 48000] # target image

#--------FOR CALIBRATION--------#
## FLAT_STAR
# please change directry names and fits ids
flat_star=irdstream.Stream2D("flat", 
                             datadir_flat, 
                             anadir, 
                             fitsid=list(range(flat_star_id[0], flat_star_id[-1]+1, 2)), 
                             band=band, 
                             inst=inst)
# aperture extraction
print(flat_star.band,' band')
if band=='h':
    trace_smf=flat_star.aptrace(cutrow = 800,nap=42) 
elif band=='y':
    trace_smf=flat_star.aptrace(cutrow = 600,nap=102) 
trace_mask = trace_smf.mask()
trace_smf.choose_aperture(fiber=mmf)

## HOTPIXEL MASK: 
# See pyird/io/read_hotpix.py for reading fixed mask (Optional)
## DARK
dark = irdstream.Stream2D('dark', 
                          datadir_dark, 
                          anadir, 
                          fitsid=list(range(dark_id[0], dark_id[-1]+1, 2)), 
                          band=band, 
                          inst=inst)
median_image = dark.immedian()
im_subbias = bias_subtract_image(median_image)
hotpix_mask = identify_hotpix_sigclip(im_subbias)
if mask_hotpix:
    kwargs_hotpix = {"hotpix_mask": hotpix_mask, "hotpix_mode": hotpix_mode}
    extin_dispcor = '_flnhp'
else:
    kwargs_hotpix = {"hotpix_mask": None, "hotpix_mode": None}
    extin_dispcor = '_fln'

## THAR
if band=='h':
    rawtag='IRDAD000' #for calib data taken by Gen2
    ignore_orders=None
elif band=='y':
    rawtag='IRDBD000' #for calib data taken by Gen2
    ignore_orders = list(np.arange(1,32))+list([50,51])
thar=irdstream.Stream2D("thar", 
                        datadir_thar, 
                        anadir, 
                        fitsid=list(range(thar_id[0], thar_id[-1]+1)), 
                        band=band, 
                        inst=inst, 
                        rawtag=rawtag)

#wavelength calibration
thar.trace = trace_smf
thar.clean_pattern(trace_mask=trace_mask,extin='', extout='_cp', hotpix_mask=hotpix_mask)
thar.calibrate_wavelength()

## FLAT
flat_star.trace = trace_smf
flat_star.clean_pattern(trace_mask=trace_mask, extin='', extout='_cp', hotpix_mask=hotpix_mask)
flat_star.imcomb = True # median combine
flat_star.flatten(**kwargs_hotpix)
df_flatn = flat_star.apnormalize(ignore_orders=ignore_orders)

#--------FOR TARGET--------#
target = irdstream.Stream2D('targets', 
                            datadir_target, 
                            anadir, 
                            fitsid=list(range(target_id[0], target_id[-1]+1, 2)), 
                            band=band, 
                            inst=inst)
target.info = True  # show detailed info
target.trace = trace_smf
target.check_aperture_position()
# clean pattern
target.clean_pattern(trace_mask=trace_mask, extin='', extout='_cp', hotpix_mask=hotpix_mask)
# flatten
target.apext_flatfield(df_flatn, **kwargs_hotpix)
# assign reference spectra & resample
target.dispcor(master_path=thar.anadir, extin=extin_dispcor)

# blaze function
flat_star.apext_flatfield(df_flatn,  **kwargs_hotpix)
flat_star.dispcor(master_path=thar.anadir)

# combine & normalize
target.normalize1D(master_path=flat_star.anadir, 
                   readout_noise_mode=readout_noise_mode, 
                   outputs=["w", "nw", "ncw"])

# save noramlized flat for fringe removal
flat_star.dispcor(master_path=thar.anadir, blaze=False)
flat_star.normalize1D(master_path=flat_star.anadir, 
                      readout_noise_mode=readout_noise_mode, 
                      outputs=["w", "nw", "ncw"])
