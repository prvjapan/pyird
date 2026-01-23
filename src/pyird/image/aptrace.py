import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.optimize import leastsq
from tqdm import tqdm
import warnings
from pyird.plot.order import plot_crosssection, plot_tracelines


def cross_section(dat, nrow, num_aperture):
    """extract cross section in row direction and search peaks.

    Args:
        dat: flat data
        nrow: row number to be extracted
        num_aperture: number of total apertures to be traced

    Returns:
        masked data and peaks at the cross section
    """
    onerow = dat[nrow]
    # search peaks
    heights = np.arange(50, 0, -5) * np.median(onerow)
    i = 0
    diffstd = 100
    peakind = []
    while (len(peakind) < num_aperture or (diffstd > 10)) and (i < len(heights)):
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


def set_aperture(dat, search_start_row, num_aperture, ign_ord=[], plot=True):
    """search aperture

    Args:
        dat: flat data
        search_start_row: starting row number to search apertures
        num_aperture: number of total apertures to be traced
        ign_ord: orders to be ignored
        plot: show figure of selected apertures

    Returns:
        peak position in the cross section at search_row
    """
    if len(ign_ord) >= 30:
       raise ValueError("length of ign_ord should be < 30.")

    # Search peaks in the cross section at search_row
    # search_row is selected so that the number of peaks and num_aperture match
    search_row_min = 500
    search_row_max = 1550
    peakind_cut = []
    search_row_lim = True
    prange = False
    search_row = search_start_row
    while ((len(peakind_cut) != num_aperture) or prange) and search_row_lim:
        onerow_masked, peakind_cut = cross_section(dat, search_row, num_aperture)
        # mask peaks for selected apertures
        if len(ign_ord)>0:
            mask = np.ones_like(peakind_cut, dtype=bool)
            mask[np.array(ign_ord)-1] = False
            peakind_cut = np.array(peakind_cut)[mask]
        if num_aperture in [21, 42]:  # h band
            prange = (peakind_cut[-1] > 1500) or (peakind_cut[0] > 40)
        elif num_aperture in [51, 102]:  # yj band
            prange = (peakind_cut[0] < 250) or (peakind_cut[-1] < 2000)
        search_row_lim = (search_row > search_row_min) and (search_row < search_row_max)
        search_row += 1
    if not search_row_lim:
        warnings.warn(
            'Please change "search_row" or "num_aperture" or confirm the orientation of the image.',
            UserWarning,
        )
        raise ValueError("Error: %d apertures could not be found." % (num_aperture))

        return
    print(f"Successfully detected the required number of apertures on detector row {search_row}.")
    if plot == True:
        # Plot to confirm the selected aperture
        pixels = np.arange(0, len(onerow_masked))
        plot_crosssection(pixels, onerow_masked, peakind_cut, dat, search_row)
    return peakind_cut, search_row


def trace_pix(dat, search_row, peakind, npix=2048):
    """trace apertures

    Args:
        dat: flat data of an order
        search_row: row number used to set aperture
        peakind: aperture (peak position) in the cross section at search_row
        npix: number of pixels

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
    for nrow in range(search_row, npix):
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
    for nrow in range(search_row - 1, 0, -1):
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
    diff = np.diff(column, prepend=column[0] - 2)
    if peakind < 40:
        useind = (
            ((0 <= diff) & (diff <= 1))
            & (row < 2000)
            & ~((column < 25) & (row < search_row))
        )  ## check!!
    elif peakind > 2000:
        useind = (
            ((0 <= diff) & (diff <= 1))
            & (row < 2000)
            & ~((column > 2020) & (row > search_row))
        )
    else:
        useind = ((0 <= diff) & (diff <= 1)) & (row < 2000)
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


def aptrace(dat, search_start_row, num_aperture, ign_ord=[], plot=True):
    """trace apertures by a polynomial function.

    Args:
        dat: flat data
        search_start_row: starting row number to search apertures
        num_aperture: number of total apertures to be traced
        ign_ord: orders to be ignored, using when num_aperture is set to a non-default value
        plot: show figure of traced apertures

    Returns:
        parameters of a polynomial to trace apertures

    Note:
        the apertures begin its search at ``search_row``, which is the row number of the detector (wavelength direction),
        and continue in the direction of increasing numbers until it matches the appropriate number of apertures.
        You may as well change the value of ``search_row`` if the aperture trace is failed.
        
    """
    if num_aperture in [42, 102]:
        print("Searching for apertures using the default num_aperture value.")
    elif num_aperture in [21, 51]:
        print(f"Searching for apertures using num_aperture = {num_aperture}.")
        print("This num_aperture value assumes a single-fiber aperture configuration on the detector.")
    else:
        warnings.warn(
            "num_aperture is not set to a default value for IRD/REACH. default values: num_aperture = 42 for H / 102 for YJ.",
            UserWarning,
        )

    peakind_cut, row = set_aperture(dat, search_start_row, num_aperture, ign_ord=ign_ord, plot=plot)

    # Trace each peak
    x, y, y0 = [], [], []
    for peakind in tqdm(peakind_cut):
        x_ord, y_ord, y0_ord = trace_pix(dat, row, peakind)
        x.append(list(x_ord))
        y.append(list(y_ord))
        y0.append(y0_ord)
    if plot == True:
        plot_tracelines(x, y)

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
    search_start_row = 1205
    num_aperture = 42
    y0, xmin, xmax, coeff = aptrace(dat, search_start_row, num_aperture)

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
