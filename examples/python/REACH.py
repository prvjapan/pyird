import numpy as np
import pkg_resources
from pyird.utils import irdstream
import pathlib
import matplotlib.pyplot as plt
from pyird.image.channel import image_to_channel_cube, channel_cube_to_image, eopixel_split, eopixel_combine
from pyird.image.bias import bias_subtract_image
from pyird.image.hotpix import identify_hotpix
import astropy.io.fits as pyf

#hotpixel mask
datadir=pathlib.Path("/home/kawahara/pyird/data/dark/")
anadir=pathlib.Path("/home/kawahara/pyird/data/dark/") 
dark=irdstream.Stream2D("targets",datadir,anadir)
dark.fitsid=[41018]    
for data in dark.rawpath:
    im = pyf.open(str(data))[0].data
im_subbias=bias_subtract_image(im)
hotpix_mask,obj=identify_hotpix(im_subbias)


#Load data
datadir=pathlib.Path("/home/kawahara/pyird/data/samples/REACH/")
anadir=pathlib.Path("/home/kawahara/pyird/data/samples/REACH/")
target=irdstream.Stream2D("targets",datadir,anadir)
#target.fitsid=[47103]
target.fitsid=[47077]

# Load an image
for datapath in target.rawpath:
    im = pyf.open(str(datapath))[0].data

#image for calibration 
calim=np.copy(im)

#read pattern removal
from pyird.io.iraf_trace import read_trace_file
from pyird.image.mask import trace
from pyird.image.trace_function import trace_legendre
from pyird.image.pattern_model import median_XY_profile

pathC=(pkg_resources.resource_filename('pyird', "data/samples/aprefC"))
path_c=(pkg_resources.resource_filename('pyird', "data/samples/apref_c"))
y0, interp_function, xmin, xmax, coeff=read_trace_file([pathC,path_c])
mask=trace(im, trace_legendre, y0, xmin, xmax, coeff)
calim[mask]=np.nan
calim[hotpix_mask]=np.nan
model_im=median_XY_profile(calim)
corrected_im=im-model_im

from pyird.plot.detector import corrected_detector
corrected_detector(im, model_im, corrected_im)

###########################################
#One Dimensionalization
###########################################

#trace
pathA=(pkg_resources.resource_filename('pyird', "data/samples/aprefB"))    
y0, interp_function, xmin, xmax, coeff=read_trace_file(pathA)

#flatten
from pyird.image.oned_extract import flatten
spec=flatten(corrected_im, trace_legendre, y0, xmin, xmax, coeff)


fig=plt.figure()
for esp in spec[5:6]:
    plt.plot(esp)
plt.show()