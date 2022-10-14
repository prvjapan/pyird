import numpy as np
import pandas as pd
from scipy.signal import medfilt, find_peaks
from scipy.optimize import leastsq
from tqdm import tqdm

from pyird.plot.order import plot_crosssection, plot_tracelines

def cross_section(dat,nrow,nap):
    """extract cross section in row direction and search peaks.

    Args:
        dat: flat data
        nrow: row number to be extracted
        nap: number of total apertures to be traced

    Returns:
        masked data and peaks at the cross section
    """
    onerow = dat[nrow]
    # search peaks
    heights=np.arange(50,0,-5)*np.median(onerow)
    i = 0
    diffstd = 100
    peakind = []
    while (len(peakind)<nap or (diffstd>10)) and (i<len(heights)):
        height = heights[i]
        peakind,_ = find_peaks(onerow,height=height,distance=5)
        diffstd = np.std(np.diff(peakind)[::2])
        i += 1
    # mask bad pixels
    maskind = np.ones(len(onerow),dtype=bool)
    maskind[peakind] = False
    peakind_cut=[]
    for i in range(len(onerow)):
        if i in peakind:
            if (onerow[i-1]>height) or (onerow[i+1]>height):
                maskind[i]=True
                peakind_cut.append(i)
    onerow_masked = np.copy(onerow)
    onerow_masked[~maskind] = np.nan
    return onerow_masked, peakind_cut

def set_aperture(dat,cutrow,nap,plot=True):
    """search aperture

    Args:
        dat: flat data
        cutrow: row number used to set aperture
        nap: number of total apertures to be traced
        plot: show figure of selected apertures

    Returns:
        peak position in the cross section at cutrow
    """
    # Search peaks in the cross section at cutrow
    # cutrow is selected so that the number of peaks and nap match
    cutrow_min=500
    cutrow_max=1550
    peakind_cut = []
    cutrow_lim = True
    prange=False
    while ((len(peakind_cut)!=nap) or prange) and cutrow_lim:
        onerow_masked,peakind_cut = cross_section(dat,cutrow,nap)
        #print('cross-section: row ',cutrow)
        if nap==42: # h band
            prange = (peakind_cut[-1]>1500) or (peakind_cut[1]>40)
        elif nap>=102: # yj band
            prange = (peakind_cut[0]<250) or (peakind_cut[-2]<2000)
        cutrow_lim = (cutrow>cutrow_min) and (cutrow<cutrow_max)
        cutrow += 1
    if not cutrow_lim:
        print('Error: %d apertures could not be found.'%(nap))
        print('Please change "cutrow" or "nap" or confirm the orientation of the image.')
        return
    print('cross-section: row ',cutrow)
    if plot==True:
        # Plot to confirm the selected aperture
        pixels = np.arange(0,len(onerow_masked))
        plot_crosssection(pixels,onerow_masked,peakind_cut)
    return peakind_cut, cutrow

def trace_pix(dat,cutrow,peakind,npix=2048):
    """trace apertures

    Args:
        dat: flat data of an order
        cutrow: row number used to set aperture
        peakind: aperture (peak position) in the cross section at cutrow
        npix: number of pixels

    Returns:
        traced pixel data
    """
    def set_newind(dat,nrow,ind):
        onerow_masked = dat[nrow]
        around = 2
        aroundpeak = onerow_masked[ind-around:ind+around+1]
        if (len(aroundpeak)==0) or np.isnan(aroundpeak[0]):
            return -1
        ind_tmp = np.argmax(aroundpeak)
        if max(aroundpeak)<3:
            return -2
        ind = ind-around+ind_tmp
        return ind

    traceind2d = pd.DataFrame([])
    ind = peakind
    for nrow in range(cutrow,npix):
        ind_new = set_newind(dat,nrow,ind)
        if ind_new==-1:
            continue
        if ind_new==-2:
            break
        ind = ind_new
        data = [nrow,ind]
        df_tmp = pd.DataFrame([data],columns=['row','column'])
        traceind2d = pd.concat([traceind2d,df_tmp],ignore_index=True)
    ind = peakind
    for nrow in range(cutrow-1,0,-1):
        ind_new = set_newind(dat,nrow,ind)
        if ind_new==-1:
            continue
        if ind_new==-2:
            break
        ind = ind_new
        data = [nrow,ind]
        df_tmp = pd.DataFrame([data],columns=['row','column'])
        traceind2d = pd.concat([traceind2d,df_tmp],ignore_index=True)

    traceind2d = traceind2d.sort_values('row',ignore_index=True)
    row = traceind2d.row.values
    column = traceind2d.column.values
    diff = np.diff(column,prepend=column[0]-2)
    if peakind < 40:
        useind = ((0<=diff) & (diff<=1)) & (row<2000) & ~((column<25) & (row<cutrow)) ## check!!
    elif peakind > 2000:
        useind = ((0<=diff) & (diff<=1)) & (row<2000) & ~((column>2020) & (row>cutrow))
    else:
        useind = ((0<=diff) & (diff<=1)) & (row<2000)
    x_ord = row[useind]
    y_ord = column[useind]

    xmid = (x_ord.min()+x_ord.max())//2
    midind = np.searchsorted(x_ord,xmid)
    y0_ord=y_ord[midind]
    return x_ord, y_ord, y0_ord

def fitfunc(x, y0, coeff):
    """Legendre function to trace.

    Args:
        x: x-array
        y0: y-offset
        coeff: Legendre polynomial coefficients

    Returns:
        Legendre function
    """
    from numpy.polynomial.legendre import legval

    x_ = np.array(x)
    xmax_ = x_.max()
    xmin_ = x_.min()
    nx = (2.*x_ - (xmax_+xmin_))/(xmax_-xmin_)
    f = legval(nx, coeff)+y0-1
    trace_lines = f
    return trace_lines

def errfunc(coeff, x, y0, data):
    """calculate error function.

    Args:
        coeff: coefficients
        x: x-array
        y0: y-offset
        data: traced pixel data to be fitted

    Returns:
        residuals of data and fitting model
    """
    model = fitfunc(x, y0, coeff)
    res = np.array(data) - np.array(model)
    return res


def fit_ord(x, y0, data):
    """optimize the fitting by using least-square method.

    Args:
        x: x-array
        y0: y-offset
        data: traced pixel data to be fitted

    Returns:
        best fit parameters
    """
    coeff = np.ones(3)
    p1, cov = leastsq(errfunc, coeff, args=(x, y0, data))
    return p1

def aptrace(dat,cutrow,nap,plot=True):
    """trace apertures by a polynomial function.

    Args:
        dat: flat data
        cutrow: row number used to set aperture
        nap: number of total apertures to be traced
        plot: show figure of traced apertures

    Returns:
        parameters of a polynomial to trace apertures
    """
    if nap!=42 and nap!=102 and nap!=104:
        print('Warning: Please set nap = 42 (for H band) or nap = 102 or 104 (for YJ band).')

    peakind_cut,row = set_aperture(dat,cutrow,nap,plot=plot)

    # Trace each peak
    x, y, y0 = [], [], []
    for peakind in tqdm(peakind_cut):
        x_ord, y_ord, y0_ord = trace_pix(dat, row, peakind)
        x.append(list(x_ord))
        y.append(list(y_ord))
        y0.append(y0_ord)
    if plot==True:
        plot_tracelines(x,y)

    # Fit each aperture
    coeff=[]
    xmin, xmax = [], []
    for i in range(len(x)):
        coeff_tmp = fit_ord(x[i],y0[i],y[i])
        coeff.append(coeff_tmp)
        x_sort = sorted(x[i])
        xmin.append(np.min(x[i]))
        xmax.append(np.max(x[i]))
    return y0, xmin, xmax, coeff

if __name__ == '__main__':
    from pathlib import Path
    from pyird.image.trace_function import trace_legendre
    import matplotlib.pyplot as plt

    # Load median combined flat
    fflat = '../../../data/flat_combine_202105.dat'
    dat = pd.read_csv(fflat,delim_whitespace=True,header=None).values
    dat = dat[::-1,::-1]

    # Trace aperture
    cutrow = 1205
    nap = 42
    y0, xmin, xmax, coeff = aptrace(dat,cutrow,nap)

    # Plot to confirm traced aperture
    x_fit = []
    for i in range(len(xmin)):
        x_tmp = np.arange(xmin[i],xmax[i]+1,1)
        x_fit.append(x_tmp)
    fit = trace_legendre(x_fit, y0, xmin, xmax, coeff)

    fig=plt.figure()
    ax=fig.add_subplot()
    for i in range(len(fit)):
        ax.plot(x_fit[i],fit[i])

    plt.show()
