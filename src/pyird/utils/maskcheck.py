import matplotlib.pyplot as plt


def plotmask(maskfits, obj, vmin=None, vmax=None):
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
