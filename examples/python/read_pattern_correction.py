import numpy as np
import pkg_resources
from pyird.utils import irdstream
import pathlib
import matplotlib.pyplot as plt
from pyird.image.channel import image_to_channel_cube, channel_cube_to_image, eopixel_split, eopixel_combine
from pyird.image.bias import bias_subtract
from pyird.image.trace_function import trace_legendre
from pyird.io.iraf_trace import read_trace_file

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

#aperture
if mode=="faint":
    pathC=(pkg_resources.resource_filename('pyird', "data/samples/aprefC"))
    path_c=(pkg_resources.resource_filename('pyird', "data/samples/apref_c"))
    y0, interp_function, xmin, xmax, coeff=read_trace_file([pathC,path_c])

    #trace
    x=[]
    for i in range(len(y0)):
        x.append(list(range(xmin[i],xmax[i]+1)))
    tl=trace_legendre(x, y0, xmin, xmax, coeff)

    import tqdm
    mask=np.zeros_like(im,dtype=bool)    
    width=2
    nx,ny=np.shape(im)
    for i in tqdm.tqdm(range(len(y0))):
        tl_tmp=np.array(tl[i],dtype=int)
        for j,ix in enumerate(x[i]):
            iys=np.max([0,tl_tmp[j]-width])
            iye=np.min([ny,tl_tmp[j]+width+2])
            mask[ix,iys:iye]=True
            
    calim=np.copy(im[::-1,::-1])
    calim[mask]=np.nan
    calim=calim[::-1,::-1]
    
#############################################################
# REMOVAL CODE


cal_channel_cube=image_to_channel_cube(calim,revert=True)
cal_eotensor=eopixel_split(cal_channel_cube)
_,nchan,xsize,ysize=np.shape(cal_eotensor)
channel_median_offset=np.nanmedian(cal_eotensor,axis=(2,3))
xprofile_offset_subtracted=np.nanmedian(cal_eotensor,axis=3)-channel_median_offset[:,:,np.newaxis]
yprofile_offset_subtracted=np.nanmedian(cal_eotensor,axis=2)-channel_median_offset[:,:,np.newaxis]

#median as model--------------------
xprofile_offset_subtracted_model=np.nanmedian(xprofile_offset_subtracted,axis=1)
yprofile_offset_subtracted_model=np.nanmedian(yprofile_offset_subtracted,axis=1)

image_pattern_model=\
    xprofile_offset_subtracted_model[:,np.newaxis,:,np.newaxis]\
    +yprofile_offset_subtracted_model[:,np.newaxis,np.newaxis,:]\
    +channel_median_offset[:,:,np.newaxis,np.newaxis]
#-----------------------------------

#original image
channel_cube=image_to_channel_cube(im,revert=True)
eotensor=eopixel_split(channel_cube)


eotensor_corrected = eotensor - image_pattern_model
channel_cube_corrected=eopixel_combine(eotensor_corrected)
corrected_im=channel_cube_to_image(channel_cube_corrected)
#############################################################



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


fig=plt.figure(figsize=(8,4))
ax1=fig.add_subplot(121)
if mode=="dark":
    vmax=-8.0;vmin=-15.0
elif mode=="faint":
    vmax=10.0;vmin=-15.0

cc=ax1.imshow(im,vmin=vmin,vmax=vmax)
plt.colorbar(cc,shrink=0.55)
ax1.set_title("Raw dark image")
ax2=fig.add_subplot(122)
cc=ax2.imshow(corrected_im,vmin=-3.0,vmax=4.0)
plt.colorbar(cc,shrink=0.55)
ax2.set_title("pattern corrected")
plt.savefig("pattern_correction.png")
plt.show()
