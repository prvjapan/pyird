import numpy as np
from pyird.utils import irdstream
import pathlib
from pyird.image.bias import bias_subtract_image
from pyird.image.hotpix import identify_hotpix_sigclip

#import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('tkagg')

#---Same as the usual reduction process (e.g. IRD_stream) up to the middle.----#
#--------SETTINGS--------#
basedir = pathlib.Path('~/pyird/data/20210317/').expanduser()

band = 'h' #'h' or 'y'
mmf = 'mmf2' #'mmf1' (comb fiber) or 'mmf2' (star fiber)

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
# aperture extraction
flat_comb=irdstream.Stream2D("flat_comb", 
                             datadir_flat, 
                             anadir, 
                             fitsid=list(range(flat_comb_id[0], flat_comb_id[-1], 2)), 
                             band=band)

if band=='h':
    trace_mmf=flat_comb.aptrace(search_start_row=1300, num_aperture=42) 
elif band=='y':
    trace_mmf=flat_comb.aptrace(search_start_row=1000, num_aperture=102) 
trace_mask = trace_mmf.mask()
trace_mmf.choose_aperture(fiber=mmf)

# hotpixel mask
dark = irdstream.Stream2D('dark', 
                        datadir_dark, 
                        anadir, 
                        fitsid=list(range(dark_id[0], dark_id[-1], 2)), 
                        band=band) # Multiple file is ok
median_image = dark.immedian()
im_subbias = bias_subtract_image(median_image)
hotpix_mask = identify_hotpix_sigclip(im_subbias)

# load ThAr raw image
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


#---The following is the process of displaying the figure----------------------#
from pyird.utils.image_widget import image_1Dand2D
import tkinter as tk

### SET PARAMETERS ###
hotpix_mask = None # comment out this if you want to show hotpixels
target.imcomb = False # set 'True' if you want to median combine images.
wavcal_path = thar.anadir/('thar_%s_%s.fits'%(thar.band,thar.trace.mmf))

## additional parameters for plot
vmin = -10
vmax = 50
scale = 'linear' # 'linear' or 'log'
params = {'vmin':vmin,'vmax':vmax,'scale':scale}

orders=[10,12] # orders to be plotted
#######################

## Values needed for plotting
rsd,wav,mask,pixcoord,rotim,iys_plot,iye_plot = target.flatten(extin='_cp',check=True,master_path=wavcal_path)

## show 1d spectrum and 2d image
for order in orders:
    print(order)
    ## draw window
    root = tk.Tk()
    window = image_1Dand2D(root,order=order,band=flat_comb.band)
    window.show_spec_to_image(rsd,wav,mask,pixcoord,rotim,iys_plot,iye_plot,master_path=wavcal_path,hotpix_mask=hotpix_mask,**params)
root.mainloop()

"""## show positions of emissions on a detector image
for order in orders:
    ## draw window
    root2 = tk.Tk()
    window2 = image_1Dand2D(root2,order=order,band=flat.band)
    window2.show_emission_position(target,rsd,wav,mask,pixcoord,rotim,iys_plot,iye_plot,wavcal_path=wavcal_path,hotpix_mask=hotpix_mask,**params)
root2.mainloop()
"""
