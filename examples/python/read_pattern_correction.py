import numpy as np
from pyird.utils import irdstream
import pathlib
import matplotlib.pyplot as plt
from pyird.image.channel import image_to_channel_cube, channel_cube_to_image, eopixel_split, eopixel_combine
from pyird.image.bias import bias_subtract

datadir=pathlib.Path("/home/kawahara/pyird/data/dark/")
anadir=pathlib.Path("/home/kawahara/pyird/data/dark/")
dark=irdstream.Stream2D("targets",datadir,anadir)
dark.fitsid=[41018]

import astropy.io.fits as pyf
for datapath in dark.rawpath:
    im = pyf.open(str(datapath))[0].data

channel_cube=image_to_channel_cube(im,revert=True)
bs_channel_cube, bias=bias_subtract(channel_cube)

#eotensor=eopixel_split(bs_channel_cube)
eotensor=eopixel_split(channel_cube)

#### CODE
_,nchan,xsize,ysize=np.shape(eotensor)
channel_median_offset=np.median(eotensor,axis=(2,3))
xprofile_offset_subtracted=np.median(eotensor,axis=3)-channel_median_offset[:,:,np.newaxis]
yprofile_offset_subtracted=np.median(eotensor,axis=2)-channel_median_offset[:,:,np.newaxis]

#median as model
xprofile_offset_subtracted_model=np.median(xprofile_offset_subtracted,axis=1)
yprofile_offset_subtracted_model=np.median(yprofile_offset_subtracted,axis=1)

image_pattern_model=\
    xprofile_offset_subtracted_model[:,np.newaxis,:,np.newaxis]\
    +yprofile_offset_subtracted_model[:,np.newaxis,np.newaxis,:]\
    +channel_median_offset[:,:,np.newaxis,np.newaxis]

#print(np.shape(image_pattern_model))
#import sys
#sys.exit()

fig=plt.figure(figsize=(10,9))
ax=fig.add_subplot(311)
for ichan in range(0,nchan):
    ax.plot(xprofile_offset_subtracted[0,ichan,:],alpha=0.1, color="C0")
ax.plot(xprofile_offset_subtracted_model[0,:],alpha=1.0, color="C0",label="even pixels")
for ichan in range(0,nchan):
    ax.plot(xprofile_offset_subtracted[1,ichan,:],alpha=0.1, color="C1")
ax.plot(xprofile_offset_subtracted_model[1,:],alpha=1.0, color="C1",label="odd pixels")
ax.set_ylabel("X-profile")
plt.legend()

ax2=fig.add_subplot(312)
for ichan in range(0,nchan):
    ax2.plot(yprofile_offset_subtracted[0,ichan,:],alpha=0.1, color="C0")
ax2.plot(yprofile_offset_subtracted_model[0,:],alpha=1.0, color="C0",label="even pixels")
for ichan in range(0,nchan):
    ax2.plot(yprofile_offset_subtracted[1,ichan,:],alpha=0.1, color="C1")
ax2.plot(yprofile_offset_subtracted_model[1,:],alpha=1.0, color="C1",label="odd pixels")
ax2.set_ylabel("Y-profile")
ax3=fig.add_subplot(313)
for ichan in range(0,nchan):
    ax3.plot(yprofile_offset_subtracted[0,ichan,:],alpha=0.1, color="C0")
ax3.plot(yprofile_offset_subtracted_model[0,:],alpha=1.0, color="C0",label="even pixels")
for ichan in range(0,nchan):
    ax3.plot(yprofile_offset_subtracted[1,ichan,:],alpha=0.1, color="C1")
ax3.plot(yprofile_offset_subtracted_model[1,:],alpha=1.0, color="C1",label="odd pixels")
ax3.set_ylabel("Y-profile")
plt.xlim(100,140)
plt.ylim(-3,3)
plt.xlabel("pixel")
plt.show()

eotensor_corrected = eotensor - image_pattern_model

####
channel_cube_corrected=eopixel_combine(eotensor_corrected)
corrected_im=channel_cube_to_image(channel_cube_corrected)

fig=plt.figure(figsize=(8,4))
ax1=fig.add_subplot(121)
cc=ax1.imshow(im,vmin=-15.0,vmax=-8.0)
plt.colorbar(cc,shrink=0.55)
ax1.set_title("Raw dark image")
ax2=fig.add_subplot(122)
cc=ax2.imshow(corrected_im,vmin=-3.0,vmax=4.0)
plt.colorbar(cc,shrink=0.55)
ax2.set_title("pattern corrected")
plt.savefig("pattern_correction.png")
plt.show()
