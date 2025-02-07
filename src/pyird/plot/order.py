import matplotlib.pyplot as plt
import numpy as np

def plot_fitresult_thar(wavsol, data, norder, npix=2048):
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
    ax1.legend(['Fitting model','ThAr emission'])
    ax1.set(title='Result of ThAr Fitting',xlabel='pixels', ylabel='wavelength [nm]')

    ax2 = fig.add_subplot(122)
    for i, data_tmp in enumerate((data.ravel())):
        if data_tmp != 0:
            pix = i % npix + 1
            res = data_tmp-wavsol[i]
            ax2.plot(pix, res, '.', color='tab:blue')
    ax2.set(title='Residuals (ThAr Wavelength - Fitted Wavelength)',xlabel='pixels', ylabel='residuals [nm]')
    # plt.savefig('wavcal_final.png')
    plt.show()

def plot_compare_used_thar(data1,wavsol1_2d,data2):
    """ compare between the wav-channel list created in the before and after iteration

    Args:
        data1: ThAr data created before iteration 
        wavsol1_2d: fitted model (i.e., wavelength solution) created by using data1
        data2: ThAr data created after iteration
    """
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)
    for i in range(len(data1[0])):
        ax.plot(wavsol1_2d[:, i], color='gray', alpha=0.5, label='Fitted 2D polynomial')
        data1_plot = np.copy(data1[:, i])
        data1_plot[data1_plot == 0] = np.nan
        ax.plot(data1_plot, '.', color='r',label='ThAr lines used for fitting')
        data2_plot = np.copy(data2[:, i])
        data2_plot[data2_plot == 0] = np.nan
        ax.plot(data2_plot, 'x', color='blue',label='Re-searched ThAr lines')
    ax.set(title='Comparison of used ThAr emission lines',xlabel='pixels', ylabel='wavelength [nm]')
    plt.show()

def plot_crosssection(pixels,flux,peakind,dat,cutrow):
    """show an image of flux at any one row extracted on a 2D detector.

    Args:
        pixels: detector pixels (column direction on detector)
        flux: flux counts
        peakind: indexes of detected peaks
        dat: 2D detector image
        cutrow: row number used to set aperture
    """
    fig, axs = plt.subplots(1, 2, figsize=(15, 5))

    axs[0].plot(pixels,flux,alpha=0.5,label="counts at 'cutrow'\n(red line in right panel)")
    axs[0].plot(pixels[peakind],flux[peakind],'x',label="peaks")

    axs[0].legend()
    axs[0].set(title="Apertures are detected!",xlabel="pixels",ylabel="counts")

    axs[1].imshow(dat,origin="lower",vmin=0,vmax=np.median(dat)*3)
    axs[1].axhline(y=cutrow, color='r', linestyle='--', label="cutrow")
    axs[1].plot(pixels[peakind], np.ones(len(peakind))*cutrow, 'x', color='tab:orange')

    axs[1].legend()
    axs[1].set(title="Detector image",xlabel="pixels",ylabel="pixels")

    plt.show()

def plot_tracelines(x, y, npix=2048):
    """show traced apertures.

    Args:
        x: position of traced pixels in row direction
        y: position of traced pixels in column direction
    """
    fig=plt.figure(figsize=(8,7))
    ax=fig.add_subplot(111)

    for i in range(len(x)):
        ax.plot(y[i],x[i],color="tab:blue")

    ax.set(xlim=(0,npix-1),ylim=(0,npix-1))
    ax.set(title="Traced Apertures on a Detector",xlabel="pixels",ylabel="pixels")
    plt.show()
