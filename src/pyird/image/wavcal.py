from pyird.io.read_linelist import read_linelist
from astropy.io import fits
import numpy as np
from scipy.signal import medfilt
from scipy.optimize import leastsq
import pandas as pd

import matplotlib.pyplot as plt

def pdat_to_wavmat(pdat, npix=2048):
    """conversion channel-wavelength data to wavelength matrix.

    Args:
        pdat: channel-wavelength data
        npix: number of detector pixels in y direction

    Return:
        channel-wavelength matrix (nipx x norder)
    """
    mat = np.zeros((npix,l-j))
    m0 = j+1
    for i in range(len(pdat)):
        order, pixel, wav = pdat.iloc[i].values
        order = order - m0
        mat[pixel,order] = wav
    return mat

def fitfunc(XY,Ni,Nx,*params):
    """calculate 2d Legendre series.

    Args:
        XY: meshgrid of (pixels, orders)
        Ni: order of the fitting function arong each echelle order
        Nx: order of the fitting function with respect to the aperture number
        params: fitting coefficients

    Return:
        wavelength of each pixels (flattened from npix x norder matrix)
    """
    pix = XY[0][0]
    m = XY[1][:,0]
    ms = XY[1]
    p= np.array(params).reshape((Ni,Nx))
    f = np.polynomial.legendre.leggrid2d(m,pix,p)/ms
    f = f.T
    return f.ravel()

def errfunc(p,XY,data,Ni,Nx):
    """calculate error function.

    Args:
        p: fitting coefficients
        XY: meshgrid of (pixels, orders)
        data: fitted data
        Ni: order of the fitting function arong each echelle order
        Nx: order of the fitting function with respect to the aperture number

    Return:
        residuals of data and fitting model
    """
    ind_data = data.ravel()!=0
    data_on = data.ravel()[ind_data]
    model = fitfunc(XY,Ni,Nx,p)
    model_on = model[ind_data]
    return data_on - model_on

def fit_wav_solution(XY,data,Ni,Nx):
    """optimize the fitting by using least-square method.
    Args:
        XY: meshgrid of (pixels, orders)
        data: fitted data
        Ni: order of the fitting function arong each echelle order
        Nx: order of the fitting function with respect to the aperture number

    Return:
        best fit parameters (coefficients of 2d legendre series)
    """
    p0 = np.ones(Ni*Nx)
    p1, cov = leastsq(errfunc,p0,args=((X,Y),data,Ni,Nx))
    p1 = p1.reshape(Ni,Nx)
    return p1

def sigmaclip(data,wavsol,N=3):
    """clipping outliers.

    Args:
        data: the reference ThAr data
        wavsol: best-fit model
        N: the number of stds to use for both the lower and upper clipping limit

    Returns:
        residuals, drop_ind
    """
    zeroind = data.ravel()==0
    residuals = (data.ravel()-wavsol.ravel())[~zeroind]
    Nsigma = N*np.std(residuals)
    n = 0
    drop_ind = []
    for i,data_tmp in enumerate((data.ravel())):
        if data_tmp != 0:
            res = data_tmp-wavsol.ravel()[i]
            if np.abs(res)>Nsigma:
                #print(n, data_tmp)
                drop_ind.append(n)
            n += 1
    return residuals,drop_ind

def plot_result(wavsol,data,npix=2048):
    """show the reference ThAr pixels and fitted model and their residuals

    Args:
        wavsol: fitting model
        data: the reference ThAr data
        npix: number of detector pixels in y direction
    """
    wavsol_2d = wavsol.reshape(npix,l-j)
    fig=plt.figure(figsize=(15,5))
    ax1=fig.add_subplot(121)
    for i in range(len(wavsol_2d[0])):
        ax1.plot(wavsol_2d[:,i],color="tab:blue")
        data_plot = np.copy(data[:,i])
        data_plot[data_plot==0] = np.nan
        ax1.plot(data_plot, ".",color="r")
        #ax.text(2050, result_plot[i][-1]-5, i, color='green',fontsize=10)
    ax1.set(xlabel="pixel",ylabel="wavelength [nm]")

    ax2=fig.add_subplot(122)
    for i,data_tmp in enumerate((data.ravel())):
        if data_tmp != 0:
            pix = i%npix + 1
            res = data_tmp-wavsol[i]
            ax2.plot(pix,res,".",color="tab:blue")
    ax2.set(xlabel="pixel",ylabel="residual [nm]")
    #plt.savefig('wavcal_final.png')
    plt.show()

if __name__ == '__main__':
    # Load a Th-Ar spectrum
    tharfile = '../../../data/thar_mmf2_a_H_20210317.fits'
    hd = fits.open(tharfile)
    dat = (hd[0].data)

    # Load a channel list
    if np.shape(dat)[0] == 21:
        chanfile = '../../../data/channel_H.list'
        norder = 21
        print('H band')
    elif np.shape(dat)[0] > 45:
        chanfile = '../../../data/channel_YJ.list'
        norder = np.shape(dat)[0]
        print('YJ band')

    pdat0 = pd.read_csv(chanfile, delimiter=',')
    j = 0
    l = norder
    offset = 0
    npix = 2048

    # allocate the line positions in pixcoord
    pdat1 = pd.DataFrame([],columns=['ORDER','CHANNEL','WAVELENGTH'])
    for i in range(j, l):
        porder = pdat0['ORDER']
        pchan = pdat0['CHANNEL']
        wav = pdat0['WAVELENGTH']

        med = medfilt(dat[i, :], kernel_size=3)

        for k, po in enumerate(porder):
            if po == i+1-offset:
                paround = 5  ## search area
                plocal_med = max(med[pchan[k]-paround:pchan[k]+paround])
                pchanlocal_med = np.where(med==plocal_med)[0]
                pchanlocal_dat = np.where(dat[i,:]==max(dat[i,:][pchanlocal_med]))[0][0]
                channel_tmp = pchanlocal_dat+1
                data = [i+1-offset,channel_tmp,wav[k]]
                df_tmp = pd.DataFrame([data],columns=pdat1.columns)
                pdat1 = pd.concat([pdat1,df_tmp],ignore_index=True)

    data1 = pdat_to_wavmat(pdat1)
    # set matrix
    pixels = np.arange(1,npix+1,1)
    orders = np.arange(1,l+1-j,1)
    X,Y = np.meshgrid(pixels,orders)
    Ni, Nx = 5,4

    # wavlength solution (rough estimation)
    p1 = fit_wav_solution((X,Y),data1,Ni,Nx)
    wavsol1 = fitfunc((X,Y),Ni,Nx,p1)
    wavsol1_2d = wavsol1.reshape(npix,l-j)

    zeroind = data1.ravel()==0
    residuals = (data1.ravel() - wavsol1)[~zeroind]
    print('standard deviation of residuals of 1st iteration = %.5f'%np.std(residuals))

    plot_result(wavsol1,data1)

    # read thar_ird2.dat
    wavref = read_linelist('../data/thar_ird2.dat')

    # add lines
    pdat2 = pd.DataFrame([],columns=['ORDER','CHANNEL','WAVELENGTH'])
    for i in range(j,l):
        pdat1_order = pdat1[pdat1['ORDER']==i+1-offset]
        wavsol1_order = wavsol1_2d[:,i]
        wavref_order = wavref[(min(wavsol1_order)<=wavref) & (wavref<=max(wavsol1_order))]

        med = medfilt(dat[i, :], kernel_size=3)
        waround = 2 + int(80000/(np.mean(wavsol1_order)/(3*np.std(residuals))))
        for ref in wavref_order:
            ind = np.searchsorted(wavsol1_order,ref)
            ind_low, ind_upp = max(ind - waround,0), min(ind+waround+1,npix)
            if np.nanmax(med[ind_low:ind_upp])>15:
                pix_med = np.where(med==np.nanmax(med[ind_low:ind_upp]))[0]
                pix_dat = np.where(dat[i,:]==max(dat[i,pix_med]))[0][0]
                pix_tmp = pix_dat + 1
                wav = ref
                data = [i+1-offset,pix_tmp,wav]
                df_tmp = pd.DataFrame([data],columns=pdat2.columns)
                pdat2 = pd.concat([pdat2,df_tmp],ignore_index=True)

    data2 = pdat_to_wavmat(pdat2)
    """ # compare previous wav-channel list & new wav-channel list
    fig=plt.figure(figsize=(7,5))
    ax=fig.add_subplot(111)
    for i in range(len(data1[0])):
        ax.plot(wavsol1_2d[:,i],color="gray",alpha=0.5)
        data1_plot = np.copy(data1[:,i])
        data1_plot[data1_plot==0] = np.nan
        ax.plot(data1_plot, ".",color="r")
        data2_plot = np.copy(data2[:,i])
        data2_plot[data2_plot==0] = np.nan
        ax.plot(data2_plot, "x",color="blue")
        #ax.text(2050, result_plot[i][-1]-5, i, color='green',fontsize=10)
    ax.set(xlabel="pixel",ylabel="wavelength [nm]")
    plt.show()
    """

    # iterate fitting polynomial function
    Ni, Nx = 5,4
    std,iter = 1,1
    while (std > 0.005) and (iter < 10):
        print(iter)
        pdat2 = pdat2[~pdat2.duplicated(keep=False,subset='CHANNEL')] # reject duplicated channel
        pdat2 = pdat2.reset_index(drop=True)
        data2 = pdat_to_wavmat(pdat2)
        p2 = fit_wav_solution((X,Y),data2,Ni,Nx)
        wavsol2 = fitfunc((X,Y),Ni,Nx,p2)
        wavsol2_2d = wavsol2.reshape(npix,l-j)
        residuals,drop_ind = sigmaclip(data2.T,wavsol2_2d.T)
        std = np.std(residuals)
        print(std)
        iter += 1
        #print(pdat2.iloc[pdat2.index[drop_ind]])
        if len(drop_ind)!=0:
            pdat2 = pdat2.drop(pdat2.index[drop_ind])

    plot_result(wavsol2,data2)
