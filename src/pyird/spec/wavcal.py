import numpy as np
from scipy.optimize import leastsq
import matplotlib.pyplot as plt

import pandas as pd
from pyird.plot.order import plot_fitresult_thar
from pyird.io.read_linelist import read_linelist
from astropy.io import fits
from scipy.signal import medfilt

import importlib

def wavcal_thar(dat, W, Ni=5, Nx=4, maxiter=10, std_threshold=0.005, channelfile_path=None, ign_ord=[]):
    """wavelegth calibration for ThAr spectrum.

    Args:
        dat: ThAr spectrum (norder x npix matrix)
        W: matrix of weights
        Ni: order of the fitting function arong each echelle order
        Nx: order of the fitting function with respect to the aperture number
        maxiter: maximum number of iterations
        std_threshold: When the std of fitting residuals reaches this value, the iteration is terminated.
        channelfile_path: path to the channel file
        ign_ord: orders to be ignored
        
    Returns:
        final results of the wavlength solution
        data of ThAr signals used for fitting

    Examples:
        >>> wavlength_solution, data = wavcal_thar(thar)
    """

    norder, npix = np.shape(dat)
    orders, channelfile, order_ref = identify_channel_mode(norder, channelfile_path, ign_ord)

    if W.shape != dat.T.shape:
        raise ValueError('The shape of weight matrix does not match the shape of the data.')

    #----- 1st identification by using the prepared channel file -----#
    # map pixels to wavelengths
    df_pixwavmap_obs = first_identification(dat, channelfile, order_ref)
    data = pixel_df_to_wav_mat(df_pixwavmap_obs, order_ref) # convert wavelength matrix

    # initialize pixel matrix
    pixels = np.arange(1, npix+1, 1)
    X, Y = np.meshgrid(pixels, orders)

    # wavlength solution from the 1st guess
    coeffs = fit_wav_solution((X, Y), data, W, Ni, Nx)
    wavlength_solution = fit_polynomial((X, Y), Ni, Nx, coeffs)
    wavlength_solution_matrix = wavlength_solution.reshape(npix, norder)

    # calculate std of residuals
    residuals = calculate_residuals(data, wavlength_solution)
    print('standard deviation of residuals (1st identification) = %.5f' %
          np.std(residuals))

    #----- 2nd identification by using the ThAr line list -----#
    # map pixels to wavelengths
    df_pixwavmap_obs = second_identification(dat, wavlength_solution_matrix, residuals, npix, norder, order_ref)
    data = pixel_df_to_wav_mat(df_pixwavmap_obs,order_ref) # convert to wavelength matrix

    #----- Iterations -----#
    print("Start iterations of ThAr fitting:")
    wavlength_solution, data = iterate_fitting(
                    X, Y, df_pixwavmap_obs, W, Ni, Nx, maxiter, std_threshold, npix, norder, order_ref)
            
    # plot result from the final iteration
    plot_fitresult_thar(wavlength_solution, data, norder)
    return wavlength_solution, data

def pixel_df_to_wav_mat(df_pixwavmap, order_ref, npix=2048):
    """conversion channel-wavelength data to wavelength matrix.

    Args:
        df_pixwavmap: channel-wavelength data
        order_ref: conventional orders

    Returns:
        channel-wavelength matrix (nipx x norder)
    """
    l_ord = len(order_ref)
    mat = np.zeros((npix, l_ord))
    for i in range(len(df_pixwavmap)):
        order, pixel, wav = df_pixwavmap.iloc[i].values
        idx_order = np.where(order_ref==order)[0][0]
        mat[pixel, idx_order] = wav
    return mat

def fit_polynomial(XY, Ni, Nx, params, poly='chebyshev'):
    """calculate 2d polynomial series.

    Args:
        XY: meshgrid of (pixels, orders)
        Ni: order of the fitting function arong each echelle order
        Nx: order of the fitting function with respect to the aperture number
        params: fitting coefficients
        poly: 'chebyshev' or 'legendre' for fitting polynomial series

    Returns:
        wavelength of each pixels (flattened from npix x norder matrix)
    """
    pix = XY[0][0]
    m = XY[1][:, 0]
    ms = XY[1]
    p = np.array(params).reshape((Ni, Nx))
    if poly=='chebyshev':
        f = np.polynomial.chebyshev.chebgrid2d(m, pix, p)/ms
    elif poly=='legendre':
        f = np.polynomial.legendre.leggrid2d(m, pix, p)/ms
    f = f.T
    return f.ravel()

def errfunc(p, XY, data, W, Ni, Nx):
    """calculate error function.

    Args:
        p: fitting coefficients
        XY: meshgrid of (pixels, orders)
        data: fitted data
        W: matrix of weights
        Ni: order of the fitting function arong each echelle order
        Nx: order of the fitting function with respect to the aperture number

    Returns:
        residuals of data and fitting model
    """
    ind_data = data.ravel() != 0
    data_on = data.ravel()[ind_data]
    model = fit_polynomial(XY, Ni, Nx, p)
    model_on = model[ind_data]
    W_on = W.ravel()[ind_data]
    return (data_on - model_on)*W_on

def fit_wav_solution(XY, data, W, Ni, Nx):
    """optimize the fitting by using least-square method.

    Args:
        XY: meshgrid of (pixels, orders)
        data: fitted data
        W: matrix of weights
        Ni: order of the fitting function arong each echelle order
        Nx: order of the fitting function with respect to the aperture number

    Returns:
        best fit parameters (coefficients of 2d legendre series)
    """
    p0 = np.ones(Ni*Nx)
    p1, cov = leastsq(errfunc, p0, args=(XY, data, W, Ni, Nx))
    p1 = p1.reshape(Ni, Nx)
    return p1

def calculate_residuals(data, wavlength_solution):
    """calculate residuals for none-zero values

    Args:
        data: fitted data
        wavelength_solution: best-fit model

    Returns:
        residuals for none-zero values
    """
    zeroind = data.ravel() == 0
    residuals = (data.ravel() - wavlength_solution.ravel())[~zeroind]
    return residuals

def sigmaclip(data, wavlength_solution, N=3):
    """clipping outliers.

    Args:
        data: the reference ThAr data
        wavlength_solution: best-fit model
        N: the number of stds to use for both the lower and upper clipping limit

    Returns:
        residuals, drop_ind
    """
    residuals = calculate_residuals(data, wavlength_solution)
    Nsigma = N*np.std(residuals)
    n = 0
    drop_ind = []
    for i, data_tmp in enumerate((data.ravel())):
        if data_tmp != 0:
            res = data_tmp-wavlength_solution.ravel()[i]
            if np.abs(res) > Nsigma:
                #print(n, data_tmp)
                drop_ind.append(n)
            n += 1
    return residuals, drop_ind

def identify_channel_mode(norder, channelfile_path=None, ign_ord=[]):
    """identify the channel model based on data shape

    Args:
        dat: ThAr spectrum (norder x npix matrix)
        channelfile_path: path to the channel file
        ign_ord: orders to be ignored

    Returns:
        diffraction orders, the reference channel-wavelength map file, and conventional orders
    """
    # Load a channel list
    if norder == 21:
        channelfile = (importlib.resources.files('pyird').joinpath('data/channel_H.list')) 
        orders = np.arange(104, 83, -1)
        order_ref = np.arange(1, 22)
        print("H band")
    elif norder == 51:
        channelfile = (importlib.resources.files('pyird').joinpath('data/channel_YJ.list'))
        orders = np.arange(158, 107, -1)
        order_ref = np.arange(1, 52)
        print("YJ band")
    else:
        if channelfile_path is None:
            raise ValueError("Cannot identify H or YJ mode. Please define the channel file for norder=", norder)
        else:
            orders, channelfile, order_ref = check_channelfile(channelfile_path, ign_ord)

    return orders, channelfile, order_ref

def check_channelfile(channelfile_path, ign_ord):
    """check the channel file format and set orders
    
    Args:
        channelfile_path: path to the channel file (user-defined)
        ign_ord: orders to be ignored
    """
    # read the channel file assuming it has the same format as data/chennel_?.lst
    try:
        df_pixwavmap_ref = pd.read_csv(channelfile_path)
    except:
        raise ValueError("The channel file format does not match that of data/channel_?.list.")
    
    # identify the band based on the wavelength range
    wav = df_pixwavmap_ref["WAVELENGTH"]
    if wav[0] > 1400:
        band = "h"
        order_default = np.arange(1, 22)
        orders = np.arange(104, 83, -1)
        print("H band")
    elif wav[0] < 1400:
        band = "y"
        order_default = np.arange(1, 52)
        orders = np.arange(158, 107, -1)
        print("YJ band")

    # set orders
    mask = np.ones_like(order_default, dtype=bool)
    mask[np.array(ign_ord)-1] = False
    order_ref = order_default[mask]
    orders = order_default[mask]
    return orders, channelfile_path, order_ref

def first_identification(dat, channelfile, order_ref, pixel_search_area=5, kernel_size=3):
    """map pixels to wavelengths by using the previously identified data with ecidentify(IRAF)

    Args:
        dat: ThAr spectrum (norder x npix matrix)
        channelfile: reference channel-wavelength map
        order_ref: conventional orders
        pixel_search_area: pixel area to search peaks around a ThAr emission in a reference spectrum
        kernel_size: kernel size for median filter

    Returns:
        channel(pixel)-wavelength map
    """
    # pixel-wavelength mapping list of ThAr emissions prepared as a reference
    df_pixwavmap_ref = pd.read_csv(channelfile, delimiter=',')
    # pixel-wavelength mapping list for the observed ThAr spectrum
    df_pixwavmap_obs = pd.DataFrame([], columns=['ORDER', 'CHANNEL', 'WAVELENGTH'])

    for k, order_tmp in enumerate(df_pixwavmap_ref['ORDER']):
        if order_tmp < 0: # order=-7 in channel_YJ.list (?)
            continue
        idx_order_tmp = np.where(order_ref == order_tmp)[0][0]
        dat_order = dat[idx_order_tmp, :]
        filtered_dat = medfilt(dat_order, kernel_size=kernel_size)
        # search peaks in the current spectrum
        ind_low = df_pixwavmap_ref['CHANNEL'][k] - pixel_search_area
        ind_upp = df_pixwavmap_ref['CHANNEL'][k] + pixel_search_area
        filtered_local_peak = np.nanmax(filtered_dat[ind_low:ind_upp])
        ind_filtered_local_peak = np.where(filtered_dat == filtered_local_peak)[0]
        if len(ind_filtered_local_peak)==0:
            continue
        local_peak = np.nanmax(dat_order[ind_filtered_local_peak])
        ind_local_peak = np.where(dat_order == local_peak)[0][0] # pixel position of a peak

        data = [order_tmp, ind_local_peak, df_pixwavmap_ref['WAVELENGTH'][k]]
        df_tmp = pd.DataFrame([data], columns=df_pixwavmap_obs.columns)
        df_pixwavmap_obs = pd.concat([df_pixwavmap_obs, df_tmp], ignore_index=True)
    return df_pixwavmap_obs

def second_identification(dat, wavlength_solution_matrix, residuals, npix, norder, order_ref, pixel_search_area=None, kernel_size=3,detect_level=80):
    """detect additional ThAr lines in the observed data with referencing the line list

    Args:
        dat: ThAr spectrum (norder x npix matrix)
        wavelength_solution_matrix: best-fit model
        residuals: residuals between data and wavelength solution
        npix: number of pixels
        norder: number of orders
        order_ref: conventional orders
        pixel_search_area: pixel area to search peaks around a ThAr emission in a line list
        kernel_size: kernel size for median filter
        detect_level: determine the lower limit of what percentage of the top data should be detected as peaks

    Returns:
        channel(pixel)-wavelength map
    """
    # read ThAr emission list (thar_ird2.dat)
    thar_linelist = (importlib.resources.files('pyird').joinpath('data/thar_ird2.dat'))
    wavref = read_linelist(thar_linelist)

    search_default = pixel_search_area is None

    # update pixel-wavelength mapping list
    df_pixwavmap = pd.DataFrame([], columns=['ORDER', 'CHANNEL', 'WAVELENGTH'])
    
    for order_tmp in order_ref:
        idx_order_tmp = np.where(order_ref == order_tmp)[0][0]
        dat_order = dat[idx_order_tmp, :]
        # filter nan values at order edges
        nanind = np.where(np.isnan(dat_order))[0]
        nanind_small_edge = nanind[nanind < 1900].max()
        nanind_large_edge = nanind[nanind > 1900].min()
        # median filter
        filtered_dat = medfilt(dat_order, kernel_size=kernel_size)
        wavlength_solution_order = wavlength_solution_matrix[:, idx_order_tmp]
        wavref_order = wavref[(min(wavlength_solution_order) <= wavref) & (wavref <= max(wavlength_solution_order))]

        if search_default:
            scale_quality = 3
            first_ident_quality = scale_quality * np.std(residuals) 
            wav_dependence = 100 * (1700 - 900)/np.mean(wavlength_solution_order) # arbitrarily tuned? (wavmin~900, wavmax~1700)
            pixel_search_area = 2 + int(wav_dependence * first_ident_quality)
        for wavref_tmp in wavref_order:
            ind_wavref_neighbor = np.searchsorted(wavlength_solution_order, wavref_tmp)
            ind_low = max(ind_wavref_neighbor - pixel_search_area, nanind_small_edge + 1)
            ind_upp = min(ind_wavref_neighbor + pixel_search_area + 1, nanind_large_edge)
            if ind_low >= ind_upp:
                continue
            # add wavelengths of a distinct peak to the list
            filtered_local_peak = np.nanmax(filtered_dat[ind_low:ind_upp])
            if filtered_local_peak > np.nanpercentile(filtered_dat[filtered_dat>0], detect_level): 
                ind_filtered_local_peak = np.where(filtered_dat == filtered_local_peak)[0]
                ind_local_peak = np.where(dat_order == np.nanmax(dat_order[ind_filtered_local_peak]))[0][0]
                data = [order_tmp, ind_local_peak, wavref_tmp]
                df_tmp = pd.DataFrame([data], columns=df_pixwavmap.columns)
                df_pixwavmap = pd.concat([df_pixwavmap, df_tmp], ignore_index=True)
    return df_pixwavmap

def iterate_fitting(X, Y, df_pixwavmap, W, Ni, Nx, maxiter, std_threshold, npix, norder, order_ref, Nsigma=1.5):
    """iterate the fitting until the std of residuals become lower than std_threshold or the number of iteration become maxiter 

    Args:
        X: grid of pixels
        Y: grid of orders
        df_pixwavmap: channel-wavelength data
        W: matrix of weights
        Ni: order of the fitting function arong each echelle order
        Nx: order of the fitting function with respect to the aperture number
        maxiter: maximum number of iterations
        std_threshold: When the std of fitting residuals reaches this value, the iteration is terminated.
        npix: number of pixels
        norder: number of orders
        order_ref: conventional orders
        Nsigma: the number of stds to use for both the lower and upper clipping limit

    Returns:
        wavelength solution and data of ThAr signals used for fitting
    """
    #Ni, Nx = 5, 4
    std, iter = 1, 1
    while (std > std_threshold) and (iter < maxiter):
        # reject duplicated channel
        df_pixwavmap = df_pixwavmap[~df_pixwavmap.duplicated(keep=False, subset='CHANNEL')]
        df_pixwavmap = df_pixwavmap.reset_index(drop=True)
        data = pixel_df_to_wav_mat(df_pixwavmap,order_ref)
        # fit again
        coeffs = fit_wav_solution((X, Y), data, W, Ni, Nx)
        wavlength_solution = fit_polynomial((X, Y), Ni, Nx, coeffs)
        wavlength_solution_matrix = wavlength_solution.reshape(npix, norder)
        residuals, drop_ind = sigmaclip(data.T, wavlength_solution_matrix.T,N=Nsigma)
        std = np.std(residuals)
        print("#",iter,"standard dev=",std)
        iter += 1
        if len(drop_ind) != 0:
            df_pixwavmap = df_pixwavmap.drop(df_pixwavmap.index[drop_ind])
    return wavlength_solution, data

def make_weight(): 
    """ weight matrix

    Note:
        REVIEW: there may be other appropreate weights...
    """
    path = (importlib.resources.files('pyird').joinpath('data/IP_fwhms_h.dat'))
    fwhms = np.loadtxt(path)
    calc_wtmp = lambda fwhm: 1/(fwhm)
    w = []
    for order in range(len(fwhms)):
        w_ord = []
        for part in range(len(fwhms[0])):
            if part != len(fwhms[0])-1:
                w_tmp = [calc_wtmp(fwhms[order][part])]*108
            else:
                w_tmp = [calc_wtmp(fwhms[order][part])]*104
            w_ord.extend(w_tmp)
        w.append(w_ord)
    w = np.array(w).T
    return w

if __name__ == '__main__':
    from pyird.plot.order import plot_fitresult_thar
    # Load a Th-Ar spectrum
    path = (importlib.resources.files('pyird').joinpath('data/thar_mmf2_a_H_20210317.fits'))
    hd = fits.open(path)
    dat = (hd[0].data)

    wavlength_solution2, data2 = wavcal_thar(dat)

    plot_fitresult_thar(wavlength_solution2, data2, np.shape(dat)[0])
