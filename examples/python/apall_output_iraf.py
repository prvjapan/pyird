from pyird.utils import irdstream
import pathlib

#--------SETTINGS--------#
basedir = pathlib.Path('~/pyird/data/20210317/').expanduser()

band = 'h' #'h' or 'y'
mmf = 'mmf2' #'mmf1' (comb fiber) or 'mmf2' (star fiber)

datadir_flat = basedir/'flat/'
anadir = basedir/'reduc/'

# last 5 digits of FITS file numbers: [start, end file number]
flat_comb_id = [41704, 41803]  # flat image for comb
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

# IRAF-like output
trace_mmf.info=True
rot180 = True if band=='h' else False
trace_mmf.write_iraf(stream=flat_comb, rot180=rot180)
