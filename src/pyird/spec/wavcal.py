import numpy as np
from scipy.optimize import leastsq
import matplotlib.pyplot as plt

import pandas as pd
from pyird.plot.order import plot_fitresult_thar
from pyird.io.read_linelist import read_linelist
from astropy.io import fits
from scipy.signal import medfilt

import pkg_resources

def wavcal_thar(dat, W, Ni=5, Nx=4, maxiter=10, std_threshold=0.005):
    """wavelegth calibration for ThAr spectrum.

    Args:
        dat: ThAr spectrum (norder x npix matrix)
        W: matrix of weights
        Ni: order of the fitting function arong each echelle order
        Nx: order of the fitting function with respect to the aperture number
        maxiter: maximum number of iterations
        std_threshold: When the std of fitting residuals reaches this value, the iteration is terminated.
        
    Returns:
        final results of the wavlength solution
        data of ThAr signals used for fitting

    Examples:
        >>> wavlength_solution, data = wavcal_thar(thar)
    """

    npix, norder, orders, channelfile = identify_channel_mode(dat)

    if W.shape != dat.T.shape:
        print('Error: weights does not match data.')
        return

    #----- 1st identification by using the prepared channel file -----#
    # map pixels to wavelengths
    df_pixel_wav_map1 = first_identification(dat, channelfile)
    data1 = pixel_df_to_wav_mat(df_pixel_wav_map1,0,norder) # convert wavelength matrix

    # initialize pixel matrix
    pixels = np.arange(1, npix+1, 1)
    X, Y = np.meshgrid(pixels, orders)

    # wavlength solution from the 1st guess
    coeffs1 = fit_wav_solution((X, Y), data1, W, Ni, Nx)
    wavlength_solution1 = fit_polynomial((X, Y), Ni, Nx, coeffs1)
    wavlength_solution1_matrix = wavlength_solution1.reshape(npix, norder)

    # calculate std of residuals
    residuals = calculate_residuals(data1, wavlength_solution1)
    print('standard deviation of residuals (1st guess) = %.5f' %
          np.std(residuals))

    #----- 2nd identification by using the ThAr line list -----#
    # map pixels to wavelengths
    df_pixel_wav_map2 = second_identification(dat, wavlength_solution1_matrix, residuals, npix, norder)
    data2 = pixel_df_to_wav_mat(df_pixel_wav_map2,0,norder) # convert to wavelength matrix

    #----- Iterations -----#
    print("Start iterations of ThAr fitting:")
    wavlength_solution2, data2 = iterate_fitting(
                    X, Y, df_pixel_wav_map2, W, Ni, Nx, maxiter, std_threshold, npix, norder)
            
    # plot result from the final iteration
    plot_fitresult_thar(wavlength_solution2, data2, norder)
    return wavlength_solution2, data2

def pixel_df_to_wav_mat(df_pixel_wav_map, j, l, npix=2048):
    """conversion channel-wavelength data to wavelength matrix.

    Args:
        df_pixel_wav_map: channel-wavelength data
        npix: number of detector pixels in y direction

    Returns:
        channel-wavelength matrix (nipx x norder)
    """
    mat = np.zeros((npix, l-j))
    m0 = j+1
    for i in range(len(df_pixel_wav_map)):
        order, pixel, wav = df_pixel_wav_map.iloc[i].values
        order = order - m0
        mat[pixel, order] = wav
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

def identify_channel_mode(dat):
    """identify the channel model based on data shape

    Args:
        dat: ThAr spectrum (norder x npix matrix)

    Returns:
        diffraction orders and the reference channel-wavelength map file
    """
    # Load a channel list
    npix = np.shape(dat)[1]
    if np.shape(dat)[0] == 21:
        channelfile = (pkg_resources.resource_filename(
            'pyird', 'data/channel_H.list'))
        norder = 21
        orders = np.arange(104, 83, -1)
        print('H band')
    elif np.shape(dat)[0] > 45:
        channelfile = (pkg_resources.resource_filename(
            'pyird', 'data/channel_YJ.list'))
        norder = np.shape(dat)[0]
        orders = np.arange(158, 107, -1)
        print('YJ band')
    else:        
        raise ValueError("Cannot identify H or YJ mode.")
    return npix, norder, orders, channelfile

def first_identification(dat, channelfile):
    """map pixels to wavelengths by using the previously identified data with ecidentify(IRAF)

    Args:
        dat: ThAr spectrum (norder x npix matrix)
        channelfile: reference channel-wavelength map

    Returns:
        channel(pixel)-wavelength map
    """
    # pixel-wavelength mapping list of ThAr emissions prepared as a reference
    df_pixel_wav_map0 = pd.read_csv(channelfile, delimiter=',')
    # pixel-wavelength mapping list for the observed ThAr spectrum
    df_pixel_wav_map1 = pd.DataFrame([], columns=['ORDER', 'CHANNEL', 'WAVELENGTH'])

    for k, order_tmp in enumerate(df_pixel_wav_map0['ORDER']):
        dat_order = dat[order_tmp - 1, :]
        filtered_dat = medfilt(dat_order, kernel_size=3)
        # search peaks in the current spectrum
        pixel_search_area = 5
        ind_low = df_pixel_wav_map0['CHANNEL'][k] - pixel_search_area
        ind_upp = df_pixel_wav_map0['CHANNEL'][k] + pixel_search_area
        filtered_local_peak = np.nanmax(filtered_dat[ind_low:ind_upp])
        ind_filtered_local_peak = np.where(filtered_dat == filtered_local_peak)[0]
        if len(ind_filtered_local_peak)==0:
            continue
        local_peak = np.nanmax(dat_order[ind_filtered_local_peak])
        ind_local_peak = np.where(dat_order == local_peak)[0][0] # pixel position of a peak

        data = [order_tmp, ind_local_peak, df_pixel_wav_map0['WAVELENGTH'][k]]
        df_tmp = pd.DataFrame([data], columns=df_pixel_wav_map1.columns)
        df_pixel_wav_map1 = pd.concat([df_pixel_wav_map1, df_tmp], ignore_index=True)
    return df_pixel_wav_map1

def second_identification(dat, wavlength_solution_matrix, residuals, npix, norder):
    """detect additional ThAr lines in the observed data with referencing the line list

    Args:
        dat: ThAr spectrum (norder x npix matrix)
        wavelength_solution_matrix: best-fit model
        residuals: residuals between data and wavelength solution
        npix: number of pixels
        norder: number of orders

    Returns:
        channel(pixel)-wavelength map
    """
    # read ThAr emission list (thar_ird2.dat)
    thar_linelist = (pkg_resources.resource_filename(
        'pyird', 'data/thar_ird2.dat'))
    wavref = read_linelist(thar_linelist)

    # update pixel-wavelength mapping list
    df_pixel_wav_map = pd.DataFrame([], columns=['ORDER', 'CHANNEL', 'WAVELENGTH'])
    
    for order_tmp in range(1,norder+1):
        dat_order = dat[order_tmp - 1, :]
        filtered_dat = medfilt(dat_order, kernel_size=3)
        wavlength_solution1_order = wavlength_solution_matrix[:, order_tmp - 1]
        wavref_order = wavref[(min(wavlength_solution1_order) <= wavref) & (wavref <= max(wavlength_solution1_order))]

        pixel_search_area = 2 + int(80000/(np.mean(wavlength_solution1_order)/(3*np.std(residuals)))) # arbitrarily tuned
        for wavref_tmp in wavref_order:
            ind_wavref_neighbor = np.searchsorted(wavlength_solution1_order, wavref_tmp)
            ind_low = max(ind_wavref_neighbor - pixel_search_area, 0)
            ind_upp = min(ind_wavref_neighbor + pixel_search_area + 1, npix)
            # add wavelengths of a distinct peak to the list
            if np.nanmax(filtered_dat[ind_low:ind_upp]) > np.nanpercentile(filtered_dat[filtered_dat>0], 80): #CHECK!!
                ind_filtered_local_peak = np.where(filtered_dat == np.nanmax(filtered_dat[ind_low:ind_upp]))[0]
                ind_local_peak = np.where(dat_order == np.nanmax(dat_order[ind_filtered_local_peak]))[0][0]
                data = [order_tmp, ind_local_peak, wavref_tmp]
                df_tmp = pd.DataFrame([data], columns=df_pixel_wav_map.columns)
                df_pixel_wav_map = pd.concat([df_pixel_wav_map, df_tmp], ignore_index=True)
    return df_pixel_wav_map

def iterate_fitting(X, Y, df_pixel_wav_map, W, Ni, Nx, maxiter, std_threshold, npix, norder):
    """iterate the fitting until the std of residuals become lower than std_threshold or the number of iteration become maxiter 

    Args:
        X: grid of pixels
        Y: grid of orders
        df_pixel_wav_map: channel-wavelength data
        W: matrix of weights
        Ni: order of the fitting function arong each echelle order
        Nx: order of the fitting function with respect to the aperture number
        maxiter: maximum number of iterations
        std_threshold: When the std of fitting residuals reaches this value, the iteration is terminated.
        npix: number of pixels
        norder: number of orders

    Returns:
        wavelength solution and data of ThAr signals used for fitting
    """
    #Ni, Nx = 5, 4
    std, iter = 1, 1
    while (std > std_threshold) and (iter < maxiter):
        # reject duplicated channel
        df_pixel_wav_map = df_pixel_wav_map[~df_pixel_wav_map.duplicated(keep=False, subset='CHANNEL')]
        df_pixel_wav_map = df_pixel_wav_map.reset_index(drop=True)
        data = pixel_df_to_wav_mat(df_pixel_wav_map,0,norder)
        # fit again
        coeffs = fit_wav_solution((X, Y), data, W, Ni, Nx)
        wavlength_solution = fit_polynomial((X, Y), Ni, Nx, coeffs)
        wavlength_solution_matrix = wavlength_solution.reshape(npix, norder)
        residuals, drop_ind = sigmaclip(data.T, wavlength_solution_matrix.T,N=1.5)
        std = np.std(residuals)
        print("#",iter,"standard dev=",std)
        iter += 1
        if len(drop_ind) != 0:
            df_pixel_wav_map = df_pixel_wav_map.drop(df_pixel_wav_map.index[drop_ind])
    return wavlength_solution, data

def make_weight(): 
    """ weight matrix

    Note:
        REVIEW: there may be other appropreate weights...
    """
    path = (pkg_resources.resource_filename('pyird', 'data/IP_fwhms_h.dat'))
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
    # Load a Th-Ar spectrum
    path = (pkg_resources.resource_filename(
        'pyird', 'data/thar_mmf2_a_H_20210317.fits'))
    hd = fits.open(path)
    dat = (hd[0].data)

    wavlength_solution2, data2 = wavcal_thar(dat)

    plot_refthar(wavlength_solution2, data2, 21)
