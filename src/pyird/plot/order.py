import matplotlib.pyplot as plt
import numpy as np

def plot_refthar(wavsol, data, norder, npix=2048):
    """show the reference ThAr pixels and fitted model and their residuals.

    Args:
        wavsol: fitting model
        data: the reference ThAr data
        norder: number of the orders
        npix: number of detector pixels in y direction
    """
    wavsol_2d = wavsol.reshape(npix, norder)
    fig = plt.figure(figsize=(15, 5))
    ax1 = fig.add_subplot(121)
    for i in range(len(wavsol_2d[0])):
        ax1.plot(wavsol_2d[:, i], color='tab:blue')
        data_plot = np.copy(data[:, i])
        data_plot[data_plot == 0] = np.nan
        ax1.plot(data_plot, '.', color='r')
        #ax.text(2050, result_plot[i][-1]-5, i, color='green',fontsize=10)
    ax1.set(xlabel='pixel', ylabel='wavelength [nm]')

    ax2 = fig.add_subplot(122)
    for i, data_tmp in enumerate((data.ravel())):
        if data_tmp != 0:
            pix = i % npix + 1
            res = data_tmp-wavsol[i]
            ax2.plot(pix, res, '.', color='tab:blue')
    ax2.set(xlabel='pixel', ylabel='residual [nm]')
    # plt.savefig('wavcal_final.png')
    plt.show()

def plot_crosssection(pixels,flux,peakind):
    """show an image of flux at any one row extracted on a 2D detector.

    Args:
        pixels: detector pixels (column direction on detector)
        flux: flux counts
        peakind: indexes of detected peaks
    """
    fig=plt.figure(figsize=(15,10))
    ax=fig.add_subplot(111)

    ax.plot(pixels,flux,alpha=0.5)
    ax.plot(pixels[peakind],flux[peakind],'x')

    plt.show()

def plot_tracelines(x,y):
    """show traced apertures.

    Args:
        x: position of traced pixels in row direction
        y: position of traced pixels in column direction
    """
    fig=plt.figure()
    ax=fig.add_subplot(111)

    for i in range(len(x)):
        ax.plot(x[i],y[i])

    plt.show()
