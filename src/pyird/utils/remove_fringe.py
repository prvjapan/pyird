import pandas as pd
import numpy as np
from astropy.timeseries import LombScargle
from astropy.stats import sigma_clip

import argparse


def ls_periodogram(wav,flux,search=None):
    """Lomb-Scargle periodogram

    Args:
        wav: wavelength [nm]
        flux: flux
        search: if True, the output includes the frequencies
                corresponding to the three strongest peaks

    Returns:
        results of periodogram
    """
    ls = LombScargle(wav,flux,normalization='psd')
    frequency,power = ls.autopower(method='fast')
    flux_ones = np.ones(len(wav))
    frequency_window,power_window = LombScargle(wav,flux_ones,fit_mean=False,center_data=False,normalization='psd').autopower(method='fast')

    if not search==None:
        search_min, search_max = search[0],search[1]
        sep_hz = 0.3 ##separate >0.3Hz
        ind = (1/search_max<frequency) & (frequency<1/search_min)
        freq_sort = frequency[ind][np.argsort(power[ind])[::-1]]
        if len(freq_sort)!=0:
            f_max = freq_sort[0]
            f_second = freq_sort[np.abs(freq_sort-freq_sort[0])>sep_hz][0]
            ind_fsec = np.where(freq_sort==f_second)[0][0]
            f_third = freq_sort[ind_fsec:][(np.abs(freq_sort[ind_fsec:]-f_max)>sep_hz) & (np.abs(freq_sort[ind_fsec:]-f_second)>sep_hz)][0]
        else:
            print('skip order ')
            f_max, f_second, f_third = 0, 0, 0
        return ls,frequency,power,frequency_window,power_window,[f_max,f_second,f_third]
    else:
        return ls,frequency,power,frequency_window,power_window

def mk_model(ls,freqs,x_fit):
    """make a sinusoidal model with given frequencies

    Args:
        ls: astropy LombScargle class
        freqs: friquencies
        x_fit: x array for fitting

    Returns:
        model and offset
    """
    offset = ls.offset()
    y_fit = offset
    for freq in freqs:
        theta = ls.model_parameters(freq)
        design_matrix = ls.design_matrix(freq,x_fit)
        y_fit += design_matrix.dot(theta)
    return y_fit,offset

def remove_fringe_order(df_flat,df_target,order,mask=True):
    """remove periodic noise for one order

    Args:
        df_flat: DataFrame of flat (normalized flat)
        df_target: DataFrame of target (normalized flux)
        order: use order
        mask: if True, mask outliers in a flat spectrum

    Returns:
        fringe removed spectrum
    """

    ## search periodic signals in the FLAT spectrum
    useind_flat = df_flat['order']==order
    wav_flat=df_flat['wav'][useind_flat][216:-320]
    flux_flat=df_flat['flux'][useind_flat][216:-320]

    if mask:
        maskind = np.abs(1-flux_flat)<0.05
        wav_flat=wav_flat[maskind]
        flux_flat=flux_flat[maskind]

    search_min, search_max = 0.05, 0.5
    ls,frequency,power,frequency_window,power_window,freqs = ls_periodogram(wav_flat,flux_flat,search=[search_min,search_max])
    if np.sum(freqs)==0:
        print(order)

    ## remove periodic signal(s) from the TARGET spectrum
    useind_target = df_target['order']==order
    wav_target = df_target['wav'][useind_target]
    flux_target = df_target['flux_median'][useind_target]
    flux_err_target = df_target['flux_err'][useind_target]
    flux_target_clipped = sigma_clip(flux_target[216:-320])
    #p = np.polyfit(wav_target,flux_target,1)
    #poly = np.poly1d(p)(wav_target)
    norm = np.median(flux_target[216:-320])

    ls_ori,frequency_ori,power_ori,_,_,freqs_ori = ls_periodogram(wav_target[216:-320],flux_target_clipped/norm,search=[search_min,search_max])
    freqs_ori_use = []
    for freq in freqs_ori:
        diff = np.abs(freqs - freq)
        useind = diff< 0.1 #[nm]#100 [um]
        if True in useind:
            freqs_ori_use.append(freq)
    if len(freqs_ori_use)>0:
        flux_target_fit,offset_target = mk_model(ls_ori,freqs_ori_use,wav_target)
        flux_rmfringe = flux_target/flux_target_fit
        flux_err_target_rmfringe = flux_err_target / flux_target_fit

        #ls_rmfringe,frequency_rmfringe,power_rmfringe,_,_,freqs_rmfringe = ls_periodogram(wav_target,flux_rmfringe)
        return flux_rmfringe, flux_err_target_rmfringe
    else:
        print(f"No prominent frequencies in Order {order}. Periodicity does not removed.")
        return flux_target, flux_err_target
