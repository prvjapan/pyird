from pyird.utils import irdstream
import pathlib
from pyird.image.bias import bias_subtract_image
from pyird.image.hotpix import identify_hotpix_sigclip

#--------SETTINGS--------#
basedir = pathlib.Path('~/pyird/data/20210317/').expanduser()

band = 'y' #'h' or 'y'
mmf = 'mmf2' #'mmf1' (comb fiber) or 'mmf2' (star fiber)
readout_noise_mode = "default" #'real' or 'default'

datadir_flat = basedir/'flat/'
datadir_dark = basedir/'dark/'
datadir_thar = basedir/'thar'
datadir_target = basedir/'target/'
anadir = basedir/'reduc/'

fitsid_flat_comb = list(range(41704,41804,2)) 
fitsid_flat_star = list(range(41804,41904,2)) 
fitsid_dark = [41504]
fitsid_thar = list(range(14632,14732))
fitsid_target = [41510]
#-------------------------#

#--------FOR CALIBRATION--------#
## FLAT_COMB
flat_comb=irdstream.Stream2D("flat_comb", datadir_flat, anadir, fitsid=fitsid_flat_comb, band=band)
# aperture extraction
print(flat_comb.band,' band')
if band=='h':
    trace_mmf=flat_comb.aptrace(cutrow = 1300,nap=42) 
elif band=='y':
    trace_mmf=flat_comb.aptrace(cutrow = 1000,nap=102) 
trace_mask = trace_mmf.mask()
trace_mmf.choose_aperture(fiber=mmf)

## HOTPIXEL MASK: 
# See pyird/io/read_hotpix.py for reading fixed mask (Optional)
## DARK
dark = irdstream.Stream2D('dark', datadir_dark, anadir, fitsid=fitsid_dark, band=band) # Multiple file is ok
median_image = dark.immedian()
im_subbias = bias_subtract_image(median_image)
hotpix_mask = identify_hotpix_sigclip(im_subbias)

## THAR (ThAr-ThAr)
if band=='h':
    rawtag='IRDAD000'
elif band=='y':
    rawtag='IRDBD000'
thar=irdstream.Stream2D("thar", datadir_thar, anadir, rawtag=rawtag, fitsid=fitsid_thar, band=band) 
thar.info = True

#wavelength calibration
thar.trace = trace_mmf
thar.clean_pattern(trace_mask=trace_mask,extin='', extout='_cp', hotpix_mask=hotpix_mask)
thar.calibrate_wavelength()

## FLAT
if mmf=='mmf2':
    ## FLAT_STAR
    flat_star=irdstream.Stream2D("flat_star", datadir_flat, anadir, fitsid=fitsid_flat_star, band=band)
    flat_star.trace = trace_mmf
    flat_star.clean_pattern(trace_mask=trace_mask, extin='', extout='_cp', hotpix_mask=hotpix_mask)
    flat_star.imcomb = True # median combine
    flat_star.flatten(hotpix_mask=hotpix_mask)
    df_flatn = flat_star.apnormalize()
elif mmf=='mmf1':
    flat_comb.trace = trace_mmf
    flat_comb.clean_pattern(trace_mask=trace_mask, extin='', extout='_cp', hotpix_mask=hotpix_mask)
    flat_comb.imcomb = True # median combine
    flat_comb.flatten(hotpix_mask=hotpix_mask)
    df_flatn = flat_comb.apnormalize()

#--------FOR TARGET--------#
target = irdstream.Stream2D('targets', datadir_target, anadir, fitsid=fitsid_target, band=band)
target.info = True  # show detailed info
target.trace = trace_mmf
# clean pattern
target.clean_pattern(trace_mask=trace_mask, extin='', extout='_cp', hotpix_mask=hotpix_mask)
# flatten
target.apext_flatfield(df_flatn, hotpix_mask=hotpix_mask)
# assign reference spectra & resample
target.dispcor(master_path=thar.anadir, extin='_flnhp')#_hp

if mmf=='mmf2':
    # blaze function
    flat_star.apext_flatfield(df_flatn, hotpix_mask=hotpix_mask)
    flat_star.dispcor(master_path=thar.anadir)

    # combine & normalize
    target.normalize1D(master_path=flat_star.anadir, readout_noise_mode=readout_noise_mode)
elif mmf=='mmf1':
    # blaze function
    flat_comb.apext_flatfield(df_flatn, hotpix_mask=hotpix_mask)
    flat_comb.dispcor(master_path=thar.anadir)

    # combine & normalize
    target.normalize1D(master_path=flat_comb.anadir, readout_noise_mode=readout_noise_mode)

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
