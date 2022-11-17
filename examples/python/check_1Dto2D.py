import numpy as np
from pyird.utils import irdstream
import pathlib
from pyird.image.bias import bias_subtract_image
from pyird.image.hotpix import identify_hotpix_sigclip
import astropy.io.fits as pyf

#import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('tkagg')

#---Same as the usual reduction process (e.g. IRD_stream) up to the middle.----#
# path
basedir = pathlib.Path('~/pyird/data/20210317/').expanduser()

### FOR CALIBRATION ###
# aperture extraction
datadir = basedir/'flat/'
anadir = basedir/'reduc/'
flat=irdstream.Stream2D("flat",datadir,anadir)
flat.fitsid=list(range(41704,41804,2)) ##FLAT_COMB

################################
### SELECT H band or YJ band ###
################################
flat.band='h' #'h' or 'y'
print(flat.band,' band')
if flat.band=='h':
    flat.fitsid_increment() # when you use H-band
    trace_mmf=flat.aptrace(cutrow = 800,nap=42) #TraceAperture instance
elif flat.band=='y':
    trace_mmf=flat.aptrace(cutrow = 1000,nap=102) #TraceAperture instance

import matplotlib.pyplot as plt
plt.imshow(trace_mmf.mask()) #apeture mask plot
plt.show()

# hotpixel mask
datadir = basedir/'dark/'
anadir = basedir/'reduc/'
dark = irdstream.Stream2D('dark', datadir, anadir,fitsid=[43814])
if flat.band=='h':
    dark.fitsid_increment() # when you use H-band
for data in dark.rawpath:
    im = pyf.open(str(data))[0].data
im_subbias = bias_subtract_image(im)
hotpix_mask = identify_hotpix_sigclip(im_subbias)

###########################
### SELECT mmf2 or mmf1 ###
###########################
trace_mmf.mmf2() #mmf2 (star fiber)
#trace_mmf.mmf1() #mmf1 (comb fiber)

# load ThAr raw image
datadir = basedir/'thar'
anadir = basedir/'reduc'
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
anadir = basedir/'reduc/'
target = irdstream.Stream2D(
    'targets', datadir, anadir, fitsid=[41510])
if flat.band=='h':
    target.fitsid_increment() # when you use H-band
target.info = True  # show detailed info
target.trace = trace_mmf
# clean pattern
target.clean_pattern(extin='', extout='_cp', hotpix_mask=hotpix_mask)


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
rsd,wav,mask,pixcoord,rotim,iys_plot,iye_plot = target.flatten_check(extin='_cp',wavcal_path=wavcal_path)

## show 1d spectrum and 2d image
for order in orders:
    print(order)
    ## draw window
    root = tk.Tk()
    root.title("Order %d"%(order))
    window = image_1Dand2D(root,order=order)
    window.show_spec_to_image(rsd,wav,mask,pixcoord,rotim,iys_plot,iye_plot,wavcal_path=wavcal_path,hotpix_mask=hotpix_mask,**params)
root.mainloop()

"""## show positions of emissions on a detector image
for order in orders:
    ## draw window
    root2 = tk.Tk()
    root2.title("Order %d"%(order))
    window2 = image_1Dand2D(root2,order=order)
    window2.show_emission_position(target,rsd,wav,mask,pixcoord,rotim,iys_plot,iye_plot,wavcal_path=wavcal_path,hotpix_mask=hotpix_mask,**params)
root2.mainloop()
"""
