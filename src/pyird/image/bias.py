"""Bias subtraction.

- originally developed by M. Kuzuhara
"""
import numpy as np
import scipy.stats as stats
from pyird.image.channel import bias_region
from pyird.image.channel import image_to_channel_cube, channel_cube_to_image
import matplotlib.pyplot as plt


def bias_subtract(channel_cube):
    """Bias subtraction for channel cube.

    Args:
       channel_cube: channel cube

    Returns:
       unbiased channel cube
       bias for channels

    Example:
        >>> #this is equivalent to image_rmbias=bias_subtract_image(im)
        >>> channel_cube=image_to_channel_cube(im)
        >>> bias_subtracted_channel_cube, channel_bias=bias_subtract(channel_cube)
        >>> image_rmbias=channel_cube_to_image(bias_subtracted_channel_cube)
    """
    snclip = 3.0  # clipping value
    ref = bias_region(channel_cube, margin=4)
    meancp = []
    for ch_num in range(np.shape(channel_cube)[0]):
        cp = stats.sigmaclip(ref[ch_num, :, :], snclip, snclip)[0]
        meancp.append(np.nanmedian(cp))
    meancp = np.array(meancp)
    unbias_channel_cube = channel_cube-meancp[:, np.newaxis, np.newaxis]
    return unbias_channel_cube, meancp


def bias_subtract_image(im):
    """Bias subtraction for image.

    Args:
       im: image

    Returns:
       bias-subtracted image

    Examples:
        >>> image_rmbias=bias_subtract_image(im)
    """
    bias_subtracted_channel_cube, channel_bias = bias_subtract(
        image_to_channel_cube(im))
    return channel_cube_to_image(bias_subtracted_channel_cube)


if __name__ == '__main__':
    import numpy as np
    from pyird.image.channel import image_to_channel_cube, channel_cube_to_image
    from pyird.utils import irdstream
    import astropy.io.fits as pyf
    import tqdm
    import pathlib

    mode = 'YJ'
    datadir = pathlib.Path('/home/kawahara/pyird/data/dark/')
    anadir = pathlib.Path('/home/kawahara/pyird/data/dark/')
    dark = irdstream.Stream2D('targets', datadir, anadir)
    dark.fitsid = [41018]
    if mode == 'H':
        dark.fitsid_increment()

    print(dark.rawpath)
    for data in tqdm.tqdm(dark.rawpath):
        im = pyf.open(str(data))[0].data
        hd = pyf.open(str(data))[0].header

    import matplotlib.pyplot as plt
    from scipy.stats import median_absolute_deviation as mad
    # channel_cube=image_to_channel_cube(im)
    #bias_subtracted_channel_cube, channel_bias=bias_subtract(channel_cube)
    # image_rmbias=channel_cube_to_image(bias_subtracted_channel_cube)
    image_rmbias = bias_subtract_image(im)

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    cc = ax1.imshow(im, vmin=-15.0, vmax=-8.0)
    plt.colorbar(cc)
    ax2 = fig.add_subplot(122)
    cc = ax2.imshow(image_rmbias, vmin=-3.0, vmax=4.0)
    plt.colorbar(cc)
    print(mad(im.flatten()), '->', mad(image_rmbias.flatten()))
    plt.show()
