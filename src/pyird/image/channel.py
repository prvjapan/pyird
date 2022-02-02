import numpy as np


def image_to_channel_cube(image, revert=True, Nch=32):
    """ conversion of an image to a channel cube
    Args:
       image: 2D image
       revert: if True reverting odd channels
       Nch: Number of the channels of a detector. For IRD Nch=32

    Returns:
       channel cube (Nch, ch_pix_num, 2048)

    """
    channel_cube = np.array(np.split(image, Nch))  # (Nch, ch_pix_num, 2048)
    if revert:
        channel_cube = revert_channel_cube(channel_cube)
    return channel_cube


def channel_cube_to_image(channel_cube, revert=True):
    """conversion of a channel cube to an image.

    Args:
       channel_cube: channel_cube
       revert: if True reverting odd channels

    Returns:
       image: 2D image
    """
    Nch, ch_pix_num, ysize = np.shape(channel_cube)
    if revert:
        image = (revert_channel_cube(channel_cube).reshape(
            Nch*ch_pix_num, ysize))  # recover
    else:
        image = (channel_cube.reshape(Nch*ch_pix_num, ysize))  # recover

    return image


def revert_channel_cube(channel_cube):
    """Revert y direction of odd channels.

    Args:
       channel cube

    Returns:
       channel cube (odd channel revereted)
    """
    channel_cube[1::2, 0::, :] = channel_cube[1::2, -1::-1, :]
    return channel_cube


def eopixel_split(channel_cube):
    """Revert y direction of odd channels.

    Args:
       channel cube

    Returns:
       eo tensor (2 [e/o] x Nchannel x xsize x ysize)
    """
    return np.array([channel_cube[:, 0::2, :], channel_cube[:, 1::2, :]])


def eopixel_combine(eop):
    """inverse operation of eopixel_split.

    Args:
       eo tensor (2 [e/o] x Nchannel x xsize x ysize)

    Returns:
       channel cube
    """
    _, l, m, n = np.shape(eop)
    channel_cube = np.zeros((l, 2*m, n))
    channel_cube[:, 0::2, :] = eop[0, :, :, :]
    channel_cube[:, 1::2, :] = eop[1, :, :, :]
    return channel_cube


def bias_region(channel_cube, margin=4):
    _, _, ysize = np.shape(channel_cube)
    ref = np.concatenate([channel_cube[:, :, 0:margin],
                         channel_cube[:, :, ysize-margin:ysize]], axis=2)
    # both margin image (Nch, ch_pix_num, 2*margin)
    return ref


if __name__ == '__main__':
    import numpy as np

    a = np.random.normal(0.0, 1.0, 2048*2048).reshape(2048, 2048)
    channel_cube = image_to_channel_cube(a, revert=True)
    eop = eopixel_split(channel_cube)
    cc = eopixel_combine(eop)
    arec = channel_cube_to_image(cc, revert=True)
    print(np.sum((a-arec)**2))
    import sys
    sys.exit()
    from pyird.image.channel import image_to_channel_cube, channel_cube_to_image
    from pyird.image.bias import bias_subtract
    from pyird.utils import irdstream
    import astropy.io.fits as pyf
    import pathlib
    import matplotlib.pyplot as plt
    datadir = pathlib.Path('/home/kawahara/pyird/data/dark/')
    anadir = pathlib.Path('/home/kawahara/pyird/data/dark/')
    dark = irdstream.Stream2D('targets', datadir, anadir)
    dark.fitsid = [41018]
    for datapath in dark.rawpath:
        im = pyf.open(str(datapath))[0].data
    channel_cube = image_to_channel_cube(im)

    bias_subtracted_channel_cube, bias = bias_subtract(channel_cube)

    if False:
        profile_o = np.median(bias_subtracted_channel_cube[7, :, :], axis=1)
        profile_e = np.median(bias_subtracted_channel_cube[8, :, :], axis=1)
        plt.plot(profile_e)
        plt.plot(profile_o)
        plt.show()

    eop = eopixel_split(bias_subtracted_channel_cube)
    pf = np.median(eop, axis=3)

    med = np.median(pf, axis=2)
    pf = pf-med[:, :, np.newaxis]
    plt.plot(np.median(pf, axis=1)[0, :], lw=2,
             alpha=1, color='C0', label='Even')
    plt.plot(np.median(pf, axis=1)[1, :], lw=2,
             alpha=1, color='C1', label='Odd')

    for i in range(0, int(np.shape(eop)[1])):
        plt.plot(pf[0, i, :], alpha=0.2, color='C0')
    for i in range(0, int(np.shape(eop)[1])):
        plt.plot(pf[1, i, :], alpha=0.2, color='C1')

    plt.legend()
    plt.xlabel('X-direction')
    plt.show()
