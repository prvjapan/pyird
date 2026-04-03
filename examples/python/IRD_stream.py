from pyird.utils import irdstream
import pathlib
from pyird.image.bias import bias_subtract_image
from pyird.image.hotpix import identify_hotpix_sigclip

#--------SETTINGS--------#
basedir = pathlib.Path('~/pyird/data/20210317/').expanduser()

band = 'h' #'h' or 'y'
mmf = 'mmf2' #'mmf1' (comb fiber) or 'mmf2' (star fiber)
mask_hotpix = True
hotpix_mode = 'nan' # 'nan' or 'interp'
readout_noise_mode = 'default' #'real' or 'default'

datadir_flat = basedir/'flat/'
datadir_dark = basedir/'dark/'
datadir_thar = basedir/'thar'
datadir_target = basedir/'target/'
anadir = basedir/'reduc/'

# last 5 digits of FITS file numbers: [start, end file number]
flat_comb_id = [41704, 41803]  # flat image for comb
flat_star_id = [41804, 41903]  # flat image for star/target
dark_id = [41504, 41505] # dark image
thar_id = [14632, 14731] # ThAr image
target_id = [41510, 41511] # target image
#-------------------------#

#--------FOR CALIBRATION--------#
## FLAT_COMB
flat_comb=irdstream.Stream2D("flat_comb", 
                             datadir_flat, 
                             anadir, 
                             fitsid=list(range(flat_comb_id[0], flat_comb_id[-1], 2)), 
                             band=band)
# aperture extraction
if band=='h':
    trace_mmf=flat_comb.aptrace(search_start_row=1300, num_aperture=42) 
elif band=='y':
    trace_mmf=flat_comb.aptrace(search_start_row=1000, num_aperture=102) 
trace_mask = trace_mmf.mask()
trace_mmf.choose_aperture(fiber=mmf)

## HOTPIXEL MASK: 
# See pyird/io/read_hotpix.py for reading fixed mask (Optional)
## DARK
dark = irdstream.Stream2D('dark', 
                        datadir_dark, 
                        anadir, 
                        fitsid=list(range(dark_id[0], dark_id[-1], 2)), 
                        band=band) # Multiple file is ok
median_image = dark.immedian()
im_subbias = bias_subtract_image(median_image)
hotpix_mask = identify_hotpix_sigclip(im_subbias)
if mask_hotpix:
    kwargs_hotpix = {"hotpix_mask": hotpix_mask, "hotpix_mode": hotpix_mode}
    extin_dispcor = '_flnhp'
else:
    kwargs_hotpix = {"hotpix_mask": None, "hotpix_mode": None}
    extin_dispcor = '_fln'

## THAR (ThAr-ThAr)
if band=='h':
    rawtag='IRDAD000'
elif band=='y':
    rawtag='IRDBD000'
thar=irdstream.Stream2D("thar", 
                        datadir_thar, 
                        anadir, 
                        rawtag=rawtag, 
                        fitsid=list(range(thar_id[0], thar_id[-1]+1)), 
                        band=band) 
thar.info = True

#wavelength calibration
thar.trace = trace_mmf
thar.clean_pattern(trace_mask=trace_mask,extin='', extout='_cp', hotpix_mask=hotpix_mask)
thar.calibrate_wavelength()

## FLAT
if mmf=='mmf2':
    ## FLAT_STAR
    flat_star=irdstream.Stream2D("flat_star", 
                                 datadir_flat, 
                                 anadir, 
                                 fitsid=list(range(flat_star_id[0], flat_star_id[-1], 2)), 
                                 band=band)
    flat_star.trace = trace_mmf
    flat_star.clean_pattern(trace_mask=trace_mask, extin='', extout='_cp', hotpix_mask=hotpix_mask)
    flat_star.imcomb = True # median combine
    flat_star.flatten(**kwargs_hotpix)
    df_flatn = flat_star.apnormalize()
elif mmf=='mmf1':
    flat_comb.trace = trace_mmf
    flat_comb.clean_pattern(trace_mask=trace_mask, extin='', extout='_cp', hotpix_mask=hotpix_mask)
    flat_comb.imcomb = True # median combine
    flat_comb.flatten(**kwargs_hotpix)
    df_flatn = flat_comb.apnormalize()

#--------FOR TARGET--------#
target = irdstream.Stream2D('targets', 
                            datadir_target, 
                            anadir, 
                            fitsid=list(range(target_id[0], target_id[-1], 2)), 
                            band=band)
target.info = True  # show detailed info
target.trace = trace_mmf
# clean pattern
target.clean_pattern(trace_mask=trace_mask, extin='', extout='_cp', hotpix_mask=hotpix_mask)
# flatten
target.apext_flatfield(df_flatn, **kwargs_hotpix)
# assign reference spectra & resample
target.dispcor(master_path=thar.anadir, extin=extin_dispcor)

if mmf=='mmf2':
    # blaze function
    flat_star.apext_flatfield(df_flatn, **kwargs_hotpix)
    flat_star.dispcor(master_path=thar.anadir)

    # combine & normalize
    target.normalize1D(master_path=flat_star.anadir, readout_noise_mode=readout_noise_mode)
elif mmf=='mmf1':
    # blaze function
    flat_comb.apext_flatfield(df_flatn, **kwargs_hotpix)
    flat_comb.dispcor(master_path=thar.anadir)

    # combine & normalize
    target.normalize1D(master_path=flat_comb.anadir, readout_noise_mode=readout_noise_mode)


#--------FOR WAVELENGTH RE-CALIBRATION--------#
datadir_comb = basedir/'comb/'
comb_id = [14733, 14832]

comb_master = irdstream.Stream2D("comb_master",
                            datadir_comb,
                            anadir,
                            fitsid=list(range(comb_id[0], comb_id[-1]+1)),
                            rawtag=rawtag,
                            band=band)
comb_master.trace = trace_mmf
comb_master.clean_pattern(trace_mask=trace_mask,
                     extin='', 
                     extout='_cp', 
                     hotpix_mask=hotpix_mask)
comb_master.imcomb = True
comb_master.apext_flatfield(df_flatn, **kwargs_hotpix, median_filter=False)
comb_master.dispcor(master_path=thar.anadir, blaze=False)

"""### hotpix (test) ###
from pyird.image.hotpix import hotpix_fits_to_dat
dark.trace = trace_mmf
dark.imcomb = True
dark.flatten(hotpix_mask=hotpix_mask)
file = dark.anadir/('%s_%s_%s.fits'%(dark.streamid,dark.band,dark.trace.mmf))
save_path = dark.anadir/('%s_%s_%s.dat'%(dark.streamid,dark.band,dark.trace.mmf))
hotpix_fits_to_dat(file,save_path)
"""
