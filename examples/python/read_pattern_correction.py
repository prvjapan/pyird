import numpy as np
import pkg_resources
from pyird.utils import irdstream
import pathlib
import matplotlib.pyplot as plt
from pyird.image.channel import image_to_channel_cube, channel_cube_to_image, eopixel_split, eopixel_combine
from pyird.image.bias import bias_subtract

mode="faint"
if mode=="dark":
    datadir=pathlib.Path("/home/kawahara/pyird/data/dark/")
    anadir=pathlib.Path("/home/kawahara/pyird/data/dark/")
    target=irdstream.Stream2D("targets",datadir,anadir)
    target.fitsid=[41018]
elif mode=="faint":
    datadir=pathlib.Path("/home/kawahara/pyird/data/samples/REACH/")
    anadir=pathlib.Path("/home/kawahara/pyird/data/samples/REACH/")
    target=irdstream.Stream2D("targets",datadir,anadir)
    #target.fitsid=[47103]
    target.fitsid=[47077]
else:
    print("no mode")
    import sys
    sys.exit()

# Load an image
import astropy.io.fits as pyf
for datapath in target.rawpath:
    im = pyf.open(str(datapath))[0].data

#image for calibration 
calim=np.copy(im)

#aperture
if mode=="faint":
    from pyird.io.iraf_trace import read_trace_file
    pathC=(pkg_resources.resource_filename('pyird', "data/samples/aprefC"))
    path_c=(pkg_resources.resource_filename('pyird', "data/samples/apref_c"))
    y0, interp_function, xmin, xmax, coeff=read_trace_file([pathC,path_c])

    from pyird.image.mask import trace
    from pyird.image.trace_function import trace_legendre
    mask=trace(im, trace_legendre, y0, xmin, xmax, coeff)
    calim[mask]=np.nan
            
    
#############################################################
# REMOVAL CODE

from pyird.image.pattern_model import median_XY_profile
model_im=median_XY_profile(calim)
corrected_im=im-model_im
print(np.sum(corrected_im))

###########################################
fig=plt.figure(figsize=(8,4))
ax1=fig.add_subplot(131)
if mode=="dark":
    vmax=-8.0;vmin=-15.0
elif mode=="faint":
    vmax=10.0;vmin=-15.0

cc=ax1.imshow(im,vmin=vmin,vmax=vmax)
plt.colorbar(cc,shrink=0.55)
ax1.set_title("Raw image")

ax2=fig.add_subplot(132)
cc=ax2.imshow(model_im,vmin=vmin,vmax=vmax)
plt.colorbar(cc,shrink=0.55)
ax2.set_title("Model image")

ax3=fig.add_subplot(133)
cc=ax3.imshow(corrected_im,vmin=-3.0,vmax=4.0)
plt.colorbar(cc,shrink=0.55)
ax3.set_title("pattern corrected")
plt.savefig("pattern_correction.png")
plt.show()
