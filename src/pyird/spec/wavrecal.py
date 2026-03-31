import numpy as np
from numpy.linalg import lstsq
import pandas as pd
import astropy.constants as const
from scipy.optimize import curve_fit
import tqdm

def generate_theoretical_wavelengths_of_lfc_lines(lam_low=950.003, dnu_comb=12.5, nline=12000):
    """Generate theoretical wavelengths of LFC lines.

    Args
    ----------
    lam_low : float
        Wavelength of the lowest comb line in nm.
    dnu_comb : float
        Frequency spacing of the comb lines in GHz.
    nline : int
        Number of comb lines to generate.

    Returns
    -------
    list of float
        Theoretical wavelengths of the comb lines in nm.
    """
    c_value = const.c.value #m/s
    
    nu0 = c_value / lam_low
    n = np.arange(nline + 1)
    nu = nu0 - n * dnu_comb
    
    lambda_th = c_value / nu
    return lambda_th.tolist()

def _gaussian_with_offset(x, y0, x0, amp, sigma):
    return y0 + amp * np.exp(-(x - x0)**2 / (2 * sigma**2))

def _fit_gaussian(wavelength, flux, params_init, flux_err=None, lower=None, upper=None, maxfev=1000):
    if lower is None:
        lower = [-np.inf, -np.inf, -np.inf, 1e-6]
    if upper is None:
        upper = [np.inf, np.inf, np.inf, 1.]

    try:
        popt, _ = curve_fit(
            _gaussian_with_offset,
            wavelength,
            flux,
            p0=params_init,
            sigma=flux_err,
            absolute_sigma=(flux_err is not None),
            bounds=(lower, upper),
            maxfev=maxfev
        )
        return popt
    except Exception:
        return params_init
    
def _fit_comb_oneord(grouped, lambda_th, ord):
    df_ord = grouped.get(ord)

    wl = df_ord['wav'].to_numpy()
    flux = df_ord['flux'].to_numpy()

    jmin = np.searchsorted(lambda_th, wl[1]) + 1
    jmax = np.searchsorted(lambda_th, wl[-1]) - 6

    jmin = max(jmin, 1)
    jmax = min(jmax, len(lambda_th) - 2)

    results_ord = []
    for j in range(jmin, jmax+1):
        lambda_th_low = (lambda_th[j-1]+lambda_th[j])/2
        lambda_th_high = (lambda_th[j]+lambda_th[j+1])/2
        mask = (wl > lambda_th_low) & (wl < lambda_th_high)
        if not np.any(mask):
            continue

        wl_j = wl[mask]
        flux_j = flux[mask]

        fit_flg = True
        params_init = [0., 
                       lambda_th[j], 
                       np.max(flux_j), 
                       (wl_j[-1] - wl_j[0]) / 5.]
        y0, x0, amp, sigma = _fit_gaussian(wl_j, flux_j, params_init)
        if (x0 < lambda_th_low) or (x0 > lambda_th_high) or (np.abs(sigma) > 0.1):
            y0, x0, amp, sigma = params_init
            fit_flg = False

        results_ord.append([ord, j, x0, amp, sigma, fit_flg])
    return results_ord

def fit_comb(df_obs, lambda_th):
    """Fit Gaussian profiles to the comb lines in each order and return a DataFrame with the fit results.
    Args
    ----------
    df_obs : pandas.DataFrame
        DataFrame containing the observed spectrum with columns 'wav', 'order', and 'flux'.
    lambda_th : list of float
        Theoretical wavelengths of the comb lines in nm.
    Returns
    -------
    pandas.DataFrame
        DataFrame with columns 'order', 'j' (the index of the comb line), 'x0' (the fitted wavelength), 'amp' (fitted amplitude), 'sigma' (fitted width), 'fit_flg' (whether the fit was successful), and 'lambda_th' (the theoretical wavelength of the comb line), 
    """
    orders = df_obs['order'].unique()

    # fit Gaussian to each comb line in each order
    grouped = {ord_: g.sort_values('wav') for ord_, g in df_obs.groupby('order')}
    results = []
    for ord in orders:
        results_ord = _fit_comb_oneord(grouped, lambda_th, ord)
        results.extend(results_ord)

    df_comblines = pd.DataFrame(results, columns=['order', 'j', 'x0', 'amp', 'sigma', 'fit_flg'])

    # Add theoretical wavelengths to the DataFrame
    lambda_fit = [lambda_th[j] for j in df_comblines['j']]
    df_comblines['lambda_fit'] = lambda_fit
    return df_comblines

def mk_weight(df_comblines):
    """Create a new DataFrame with the same structure as df_comblines, but with negative amplitudes set to zero and failed fits filled with neighboring valid amplitudes.

    Args
    ----------
    df_comblines : pandas.DataFrame
        DataFrame containing the fit results for the comb lines, with columns 'order', 'j', 'x0', 'amp', 'sigma', 'fit_flg', and 'lambda_th'.
    Returns
    -------
    pandas.DataFrame
        A new DataFrame with the same columns as df_comblines, but with negative amplitudes set to zero and failed fits filled with neighboring valid amplitudes.
    Notes
    -----
    - Negative amplitudes are set to zero because they are not physically meaningful in this context.
    - Failed fits (where 'fit_flg' is False) are filled with the average of the nearest valid amplitudes from the same order. If there are no valid neighbors, the amplitude is set to zero.
    """
    df_linelist = df_comblines.copy()

    amp = df_linelist['amp']
    fit = df_linelist['fit_flg']

    # Set negative amplitudes to zero
    df_linelist.loc[fit & (amp < 0), 'amp'] = 0

    # Fill failed fits with neighboring valid amplitudes
    amp_prev = df_comblines['amp'].shift(1)
    amp_next = df_comblines['amp'].shift(-1)
    fit_prev = fit.shift(1, fill_value=False)
    fit_next = fit.shift(-1, fill_value=False)

    fail = ~fit
    prev_valid = fit_prev & (amp_prev > 0)
    next_valid = fit_next & (amp_next > 0)

    filled_amp = pd.Series(0.0, index=df_comblines.index)
    filled_amp = filled_amp.mask(next_valid, amp_next)
    filled_amp = filled_amp.mask(prev_valid, amp_prev)
    filled_amp = filled_amp.mask(prev_valid & next_valid, (amp_prev + amp_next) / 2)

    df_linelist.loc[fail, 'amp'] = filled_amp[fail]

    return df_linelist

def _poly_lfit(x, y, sig, ma):
    """Fit a polynomial of degree ma to the data (x, y) with weights sig. 
    The model is: y = a0 + a1*x + a2*x^2 + ... + a_{ma-1}*x^{ma-1}

    Args
    ----------
    x, y, sig : array-like
        Fitting data and their uncertainties. Must have the same length.
    ma : int
        Degree of the polynomial + 1 (i.e., number of coefficients). Must be >= 1.

    Returns
    -------
    coeffs : ndarray
        Fitted polynomial coefficients, where coeffs[0] is the constant term, coeffs[1] is the linear term, etc.
    cov : ndarray
        Covariance matrix of the fitted coefficients.
    chisq : float
        Chi-squared of the fit.
    A : ndarray
        Design matrix used for the fit (Vandermonde matrix of x).
    x_mean : float
        Mean of x used for centering.
    x_scale : float
        Scale of x used for standardization.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    sig = np.asarray(sig, dtype=float)

    if not (len(x) == len(y) == len(sig)):
        raise ValueError("x, y, sig must have the same length")
    if ma < 1:
        raise ValueError("ma must be >= 1")

    # Remove data points where x, y, or sig are not finite, or where sig is non-positive
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(sig) & (sig > 0)

    x = x[mask]
    y = y[mask]
    sig = sig[mask]
    #print(f"Removed {np.sum(~mask)} invalid data points; {len(x)} remain for fitting")

    if len(x) == 0:
        raise ValueError("No valid data points remain after removing nan/inf and invalid sig")
    if len(x) < ma:
        raise ValueError(f"Not enough valid data points for fitting: need at least {ma}, got {len(x)}")

    # standardize x to improve numerical stability
    x_mean = np.mean(x)
    x_scale = np.std(x)

    if not np.isfinite(x_scale) or x_scale == 0.0:
        raise ValueError("x has zero variance after filtering; cannot standardize")

    z = (x - x_mean) / x_scale

    # construct the design matrix A for a polynomial fit, where A[i, j] = z[i]**j
    A = np.vander(z, N=ma, increasing=True)

    Aw = A / sig[:, None]
    yw = y / sig

    coeffs, _, rank, _ = lstsq(Aw, yw, rcond=None)

    if rank < ma:
        raise np.linalg.LinAlgError(
            f"Design matrix is rank deficient: rank={rank}, expected={ma}"
        )

    ATA = Aw.T @ Aw
    cov = np.linalg.inv(ATA)
    chisq = np.sum(((y - A @ coeffs) / sig) ** 2)

    return coeffs, cov, chisq, A, x_mean, x_scale

def _poly_eval_scaled(x, coeffs, x_mean, x_scale):
    x = np.asarray(x, dtype=float)
    z = (x - x_mean) / x_scale
    A = np.vander(z, N=len(coeffs), increasing=True)
    return A @ coeffs

def _reject_outliers_by_group(
    pfit,
    lambda_fit,
    weight,
    x_mean,
    x_scale,
    par,
    n_poly,
    rej_lev1,
    group_size=5,
):
    """
    Reject outliers in the fit by dividing the data into groups of size group_size, calculating the residuals of the fit within each group, and rejecting points that deviate from the median residual by more than rej_lev1 times the MAD of the residuals in that group. The last group may have fewer than group_size points, and will be rejected using the median and MAD from the last complete group.

    Args
    ----------
    pfit : array-like
    lambda_fit : array-like
    weight : array-like
    par : array-like
        par[0] + par[1] * x + ... + par[n_poly] * x**n_poly
    n_poly : int
    rej_lev1 : float
    group_size : int, default=5

    Returns
    -------
    pused : ndarray
    lambda_used : ndarray
    weight_used : ndarray
    lambda_model : ndarray
    used_mask : ndarray of bool
    """
    lambda_fit = np.asarray(lambda_fit, dtype=float)
    weight = np.asarray(weight, dtype=float)
    par = np.asarray(par, dtype=float)

    if group_size < 1:
        raise ValueError("group_size must be >= 1")
    if n_poly < 0:
        raise ValueError("n_poly must be >= 0")
    if len(par) < n_poly + 1:
        raise ValueError("len(par) must be at least n_poly + 1")

    data_fit = len(lambda_fit)
    ngroup = data_fit // group_size
    remainder = data_fit % group_size

    #print(f"data_fit - ngroup*{group_size} = {remainder}")

    # Evaluate the polynomial model at the pfit values
    lambda_model = _poly_eval_scaled(pfit, par[: n_poly + 1], x_mean, x_scale)

    used_mask = np.zeros(data_fit, dtype=bool)

    last_median = None
    last_mad = None

    # Reject outliers in complete groups
    for i in range(ngroup):
        sl = slice(i * group_size, (i + 1) * group_size)

        residual = lambda_model[sl] - lambda_fit[sl]
        median = np.median(residual)
        mad = np.median(np.abs(residual - median))

        last_median = median
        last_mad = mad

        if mad == 0:
            mask = np.isclose(residual, median)
        else:
            mask = np.abs(residual - median) < rej_lev1 * mad
        used_mask[sl] = mask

    # Reject outliers in the last incomplete group using the median and MAD from the last complete group
    if remainder > 0:
        if last_median is None or last_mad is None:
            raise ValueError(
                "remainder exists but there is no complete group to provide median/MAD"
            )

        start = ngroup * group_size
        for idx in range(start, data_fit):
            residual_i = lambda_model[idx] - lambda_fit[idx]
            if last_mad == 0:
                used_mask[idx] = np.isclose(residual_i, last_median)
            else:
                used_mask[idx] = abs(residual_i - last_median) < rej_lev1 * last_mad

    pused = pfit[used_mask]
    lambda_used = lambda_fit[used_mask]
    weight_used = weight[used_mask]

    return pused, lambda_used, weight_used, lambda_model, used_mask

def _convert_lambda_to_pixel(lambda_obs, lambda_peak):
    idx = np.searchsorted(lambda_obs, lambda_peak) - 1
    valid = (idx >= 0) & (idx < len(lambda_obs) - 1)
    idx = idx[valid]
    lambda_peak_valid = lambda_peak[valid]

    pix_obs = 1 + idx + (
        (lambda_peak_valid - lambda_obs[idx]) /
        (lambda_obs[idx + 1] - lambda_obs[idx])
    )
    return pix_obs

def _prepare_df_fit(df_obs, df_linelist):
    df_fit = df_linelist.copy()

    # set weight as 1/sqrt(amp)
    df_fit['weight'] = 1.0 / np.sqrt(np.abs(df_fit['amp']))

    # convert lambda_obs to pixel using the theoretical lambda of the comb lines as reference
    orders = df_obs['order'].unique()
    for ord in orders:
        df_obs_ord = df_obs[df_obs['order'] == ord]
        df_linelist_ord = df_linelist[df_linelist['order'] == ord]
        lambda_obs = df_obs_ord['wav'].to_numpy()
        lambda_peak = df_linelist_ord['x0'].to_numpy()

        pix_obs = _convert_lambda_to_pixel(lambda_obs, lambda_peak)
        df_fit.loc[df_fit['order'] == ord, 'pix'] = pix_obs

    # mask out failed fits and negative amplitudes by setting pix, x0, weight to nan
    mask_fit = df_fit['fit_flg'] & (df_fit['amp']>10.0)
    nan_columns = ['pix', 'x0', 'weight']
    for col in nan_columns:
        df_fit.loc[~mask_fit, col] = np.nan
    return df_fit

def update_wavelength(df_obs, df_linelist, n_poly=6):
    """Update the wavelength solution for each order.
    Args
    ----------
    df_obs : pandas.DataFrame
        DataFrame containing the observed spectrum with columns 'wav', 'order', and 'flux'.
    df_linelist : pandas.DataFrame
        DataFrame containing the fit results for the comb lines, with columns 'order', 'j', 'x0', 'amp', 'sigma', 'fit_flg', and 'lambda_th'.
    n_poly : int
        Degree of the polynomial + 1 (i.e., number of coefficients) to fit for the wavelength solution. Default is 6 (i.e., a 5th degree polynomial). 

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns 'lambda_model' (updated wavelength solution), 'order', and 'flux',
    """
    orders = df_obs['order'].unique()

    df_fit = _prepare_df_fit(df_obs, df_linelist)

    df_recalib = pd.DataFrame(columns=['lambda_model', 'order', 'flux'])
    for ord in tqdm.tqdm(orders):
        df_obs_ord = df_obs[df_obs['order'] == ord]
        lambda_obs = df_obs_ord['wav'].to_numpy()
        flux_obs = df_obs_ord['flux'].to_numpy()
        pix_obs = np.arange(1, len(flux_obs) + 1, dtype=float)

        df_fit_ord = df_fit[df_fit['order'] == ord]
        pix_fit = df_fit_ord['pix'].to_numpy()
        lambda_fit = df_fit_ord['lambda_fit'].to_numpy()
        weight = df_fit_ord['weight'].to_numpy()

        if (19 <= ord <= 50) or (53 <= ord <= 71):
            niter = 0
            while niter<2:
                coeffs, cov, chisq, A, x_mean, x_scale = _poly_lfit(pix_fit, lambda_fit, weight, ma=n_poly+1)
                pix_fit, lambda_fit, weight, lambda_model, used_mask = _reject_outliers_by_group(pix_fit, lambda_fit, weight, x_mean, x_scale, coeffs, n_poly, rej_lev1=4.5)
                niter += 1
            lambda_model = _poly_eval_scaled(pix_obs, coeffs, x_mean, x_scale)
            df_recalib_ord = pd.DataFrame({'lambda_model': lambda_model, 'order': ord, 'flux': flux_obs})

            df_recalib = pd.concat([df_recalib, df_recalib_ord], ignore_index=True)
        else:
            df_recalib_ord = pd.DataFrame({'lambda_model': lambda_obs, 'order': ord, 'flux': flux_obs})
            mask = flux_obs > 1e5
            df_recalib_ord.loc[mask, 'flux'] = 0.0
            df_recalib = pd.concat([df_recalib, df_recalib_ord], ignore_index=True)
    
    return df_recalib