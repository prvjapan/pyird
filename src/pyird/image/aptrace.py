import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.optimize import leastsq
from tqdm import tqdm
import warnings
from pyird.plot.order import plot_crosssection, plot_tracelines


def cross_section(dat, nrow, nap):
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
    heights = np.arange(50, 0, -5) * np.abs(np.median(onerow))
    i = 0
    diffstd = 100
    peakind = []
    while (len(peakind) < nap or (diffstd > 10)) and (i < len(heights)):
        height = heights[i]
        peakind, _ = find_peaks(onerow, height=height, distance=5)
        diffstd = np.std(np.diff(peakind)[::2])
        i += 1
    # mask bad pixels
    maskind = np.ones(len(onerow), dtype=bool)
    maskind[peakind] = False
    peakind_cut = []
    for i in range(len(onerow)):
        if i in peakind:
            if (onerow[i - 1] > height) or (onerow[i + 1] > height):
                maskind[i] = True
                peakind_cut.append(i)
    onerow_masked = np.copy(onerow)
    onerow_masked[~maskind] = np.nan
    return onerow_masked, peakind_cut


def set_aperture(dat, cutrow, nap, ign_ord=[], plot=True):
    """search aperture

    Args:
        dat: flat data
        cutrow: row number used to set aperture
        nap: number of total apertures to be traced
        ign_ord: orders to be ignored
        plot: show figure of selected apertures

    Returns:
        peak position in the cross section at cutrow
    """
    if len(ign_ord) >= 30:
       raise ValueError("length of ign_ord should be < 30.")

    # Search peaks in the cross section at cutrow
    # cutrow is selected so that the number of peaks and nap match
    npix = dat.shape[0]
    cutrow_min = int(500 * npix/2048)
    cutrow_max = int(1550 * npix/2048)
    if (cutrow < cutrow_min) or (cutrow_max < cutrow):
        raise ValueError('Error: please set the value of "cutrow" between %d and %d.' % (cutrow_min, cutrow_max))
    
    peakind_cut = []
    cutrow_lim = True
    prange = False
    while ((len(peakind_cut) != nap) or prange) and cutrow_lim:
        onerow_masked, peakind_cut = cross_section(dat, cutrow, nap)
        # mask peaks for selected apertures
        if len(ign_ord)>0:
            mask = np.ones_like(peakind_cut, dtype=bool)
            mask[np.array(ign_ord)-1] = False
            peakind_cut = np.array(peakind_cut)[mask]
        if nap in [21, 42]:  # h band
            prange = (peakind_cut[-1] > 1500) or (peakind_cut[0] > 40)
        elif nap in [51, 102]:  # yj band
            prange = (peakind_cut[0] < 250) or (peakind_cut[-1] < 2000)
        cutrow_lim = (cutrow > cutrow_min) and (cutrow < cutrow_max)
        cutrow += 1
    if not cutrow_lim:
        warnings.warn(
            'Please change "cutrow" or "nap" or confirm the orientation of the image.',
            UserWarning,
        )
        raise ValueError("Error: %d apertures could not be found." % (nap))

    print("cross-section: row ", cutrow)
    if plot == True:
        # Plot to confirm the selected aperture
        pixels = np.arange(0, len(onerow_masked))
        plot_crosssection(pixels, onerow_masked, peakind_cut, dat, cutrow)
    return peakind_cut, cutrow


def trace_pix(dat, cutrow, peakind, npix=2048, trace_lim=[0, 1]):
    """trace apertures

    Args:
        dat: flat data of an order
        cutrow: row number used to set aperture
        peakind: aperture (peak position) in the cross section at cutrow
        npix: number of pixels
        trace_lim: limit of the difference between pixels in the traced aperture

    Returns:
        traced pixel data
    """

    def set_newind(dat, nrow, ind):
        onerow_masked = dat[nrow]
        around = 2
        aroundpeak = onerow_masked[ind - around : ind + around + 1]
        if (len(aroundpeak) == 0) or np.isnan(aroundpeak[0]):
            return -1
        ind_tmp = np.argmax(aroundpeak)
        if max(aroundpeak) < 3:
            return -2
        ind = ind - around + ind_tmp
        return ind

    traceind2d = pd.DataFrame([])
    ind = peakind
    for nrow in range(cutrow, npix):
        ind_new = set_newind(dat, nrow, ind)
        if ind_new == -1:
            continue
        if ind_new == -2:
            break
        ind = ind_new
        data = [nrow, ind]
        df_tmp = pd.DataFrame([data], columns=["row", "column"])
        traceind2d = pd.concat([traceind2d, df_tmp], ignore_index=True)
    ind = peakind
    for nrow in range(cutrow - 1, -1, -1):
        ind_new = set_newind(dat, nrow, ind)
        if ind_new == -1:
            continue
        if ind_new == -2:
            break
        ind = ind_new
        data = [nrow, ind]
        df_tmp = pd.DataFrame([data], columns=["row", "column"])
        traceind2d = pd.concat([traceind2d, df_tmp], ignore_index=True)

    traceind2d = traceind2d.sort_values("row", ignore_index=True)
    row = traceind2d.row.values
    column = traceind2d.column.values
    diff = np.diff(column, prepend=column[0] - 1)
    if peakind < 40:
        useind = (
            ((0 <= diff) & (diff <= 1))
            & (row < 2000)
            & ~((column < 25) & (row < cutrow))
        )  ## check!!
    elif peakind > 2000:
        useind = (
            ((0 <= diff) & (diff <= 1))
            & (row < 2000)
            & ~((column > 2020) & (row > cutrow))
        )
    else:
        useind = ((trace_lim[0] <= diff) & (diff <= trace_lim[1])) & (row < 2000)
    x_ord = row[useind]
    y_ord = column[useind]

    xmid = (x_ord.min() + x_ord.max()) // 2
    midind = np.searchsorted(x_ord, xmid)
    y0_ord = y_ord[midind]
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
    nx = (2.0 * x_ - (xmax_ + xmin_)) / (xmax_ - xmin_)
    f = legval(nx, coeff) + y0 - 1
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


def aptrace(dat, cutrow, nap, ign_ord=[], plot=True):
    """trace apertures by a polynomial function.

    Args:
        dat: flat data
        cutrow: row number used to set aperture
        nap: number of total apertures to be traced
        ign_ord: orders to be ignored, using when nap is set to a non-default value
        plot: show figure of traced apertures

    Returns:
        parameters of a polynomial to trace apertures

    Note:
        the apertures begin its search at ``cutrow``, which is the row number of the detector (wavelength direction),
        and continue in the direction of increasing numbers until it matches the appropriate number of apertures.
        You may as well change the value of ``cutrow`` if the aperture trace is failed.
        
    """
    if nap in [42, 102]:
        print("default nap value")
    elif nap in [21, 51]:
        warnings.warn("Looks a single fiber aperture on the detector.", UserWarning)
    else:
        warnings.warn(
            '"nap" is not default value. default: nap = 42 for H / 102 for YJ if you analyse IRD or REACH data.',
            UserWarning,
        )

    peakind_cut, row = set_aperture(dat, cutrow, nap, ign_ord=ign_ord, plot=plot)

    # Trace each peak
    npix = dat.shape[0]
    x, y, y0 = [], [], []
    for peakind in tqdm(peakind_cut):
        x_ord, y_ord, y0_ord = trace_pix(dat, row, peakind, npix=npix)
        x.append(list(x_ord))
        y.append(list(y_ord))
        y0.append(y0_ord)
    if plot == True:
        plot_tracelines(x, y, npix=npix)
    # Fit each aperture
    coeff = []
    xmin, xmax = [], []
    for i in range(len(x)):
        coeff_tmp = fit_ord(x[i], y0[i], y[i])
        coeff.append(coeff_tmp)
        x_sort = sorted(x[i])
        xmin.append(np.min(x[i]))
        xmax.append(np.max(x[i]))
    return y0, xmin, xmax, coeff


if __name__ == "__main__":
    from pathlib import Path
    from pyird.image.trace_function import trace_legendre
    import matplotlib.pyplot as plt

    # Load median combined flat
    fflat = "../../../data/flat_combine_202105.dat"
    dat = pd.read_csv(fflat, delim_whitespace=True, header=None).values
    dat = dat[::-1, ::-1]

    # Trace aperture
    cutrow = 1205
    nap = 42
    y0, xmin, xmax, coeff = aptrace(dat, cutrow, nap)

    # Plot to confirm traced aperture
    x_fit = []
    for i in range(len(xmin)):
        x_tmp = np.arange(xmin[i], xmax[i] + 1, 1)
        x_fit.append(x_tmp)
    fit = trace_legendre(x_fit, y0, xmin, xmax, coeff)

    fig = plt.figure()
    ax = fig.add_subplot()
    for i in range(len(fit)):
        ax.plot(x_fit[i], fit[i])

    plt.show()
