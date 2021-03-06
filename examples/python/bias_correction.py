from pyird.image.bias import bias_subtract
from pyird.image.channel import image_to_channel_cube, channel_cube_to_image
import matplotlib.pyplot as plt
import astropy.io.fits as pyf
from pyird.utils import irdstream
import pathlib
datadir = pathlib.Path('/home/kawahara/pyird/data/dark/')
anadir = pathlib.Path('/home/kawahara/pyird/data/dark/')
dark = irdstream.Stream2D('targets', datadir, anadir)
dark.fitsid = [41018]


for datapath in dark.rawpath:
    im = pyf.open(str(datapath))[0].data

plt.imshow(im, vmin=-15, vmax=-8)
plt.title('Raw dark image')
plt.show()

channel_cube = image_to_channel_cube(im)
plt.imshow(channel_cube[0, :, :], vmin=-15, vmax=-8)
plt.title('Channel 0')
plt.show()

bias_subtracted_channel_cube, bias = bias_subtract(channel_cube)

bias_subtracted_im = channel_cube_to_image(bias_subtracted_channel_cube)

fig = plt.figure(figsize=(8, 4))
ax1 = fig.add_subplot(121)
cc = ax1.imshow(im, vmin=-15.0, vmax=-8.0)
plt.colorbar(cc, shrink=0.55)
ax1.set_title('Raw dark image')
ax2 = fig.add_subplot(122)
cc = ax2.imshow(bias_subtracted_im, vmin=-3.0, vmax=4.0)
plt.colorbar(cc, shrink=0.55)
ax2.set_title('Bias corrected')
plt.savefig('bias_correction.png')
plt.show()
