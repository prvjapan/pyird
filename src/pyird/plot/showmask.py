import matplotlib.pyplot as plt
import numpy as np


def show_hotpix(obj, im):
    """show hotpixels.

    Args:
       obj: sep object
       im: image or mask
    """
    from matplotlib.patches import Ellipse

    # plot background-subtracted image
    fig, ax = plt.subplots()
    m, s = np.mean(im), np.std(im)
    ax.imshow(im, interpolation='nearest', cmap='gray',
              vmin=m-5*s, vmax=m+5*s, origin='lower')

    # plot an ellipse for each object
    for i in range(len(obj)):
        e = Ellipse(xy=(obj['x'][i], obj['y'][i]),
                    width=6*obj['a'][i],
                    height=6*obj['b'][i],
                    angle=obj['theta'][i] * 180. / np.pi)
        e.set_facecolor('none')
        e.set_edgecolor('red')
        ax.add_artist(e)
    plt.show()

def show_maskedpix(im, mask,vmin=-10,vmax=50):
    """show masked pixels.

    Args:
       im: image
       mask: mask
       vmin: min value of imshow
       vmax: max value of imshow
    """
    fig=plt.figure()
    ax=fig.add_subplot()

    ax.imshow(im,vmin=vmin,vmax=vmax,origin='lower',cmap='OrRd')
    ## frames of masked pixels are colored
    fac=10
    condition = np.kron(mask > 0, np.ones((fac, fac)))
    extent = (-0.5, mask.shape[1]-0.5, -0.5, mask.shape[0]-0.5)
    ax.contour(condition, levels=[0.5], extent=extent,cmap='cool',linewidths=[1],origin='lower')

    plt.show()

def show_mask(maskfits, obj, vmin=None, vmax=None):
    """
    Args:
        maskfits: fitsset of mask
        obj: fitsset of target/flat
    """
    from scipy.stats import median_absolute_deviation as mad
    import numpy as np
    darr = maskfits.data()
    cimg = obj.data()
    dimg = np.copy(cimg)
    mmask = darr == 0
    cimg[mmask] = None
    mmask = darr > 0
    dimg[mmask] = None

    mad((dimg).flatten())
    np.median((dimg).flatten())

    fig = plt.figure(figsize=(20, 10))
    fig.add_subplot(121)
    c = plt.imshow(cimg, vmin=vmin, vmax=vmax)
    plt.colorbar(c, shrink=0.5)
    fig.add_subplot(122)
    d = plt.imshow(dimg, vmin=vmin, vmax=vmax)
    plt.colorbar(d, shrink=0.5)
    plt.show()