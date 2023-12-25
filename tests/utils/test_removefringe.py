import pytest
import numpy as np
import pandas as pd
from pyird.utils.remove_fringe import remove_fringe_order,ls_periodogram

def test_removefringe():
    def multi_sin(x,freqs,amps):
        y=1
        for i,freq in enumerate(freqs):
            y += amps[i]*np.sin(2*np.pi*(x*freq))
        return y

    def mk_mockdata(freqs,amps,wavlim,npix,noise=0.01):
        wav = np.linspace(wavlim[0],wavlim[1],npix)
        flux_clean = multi_sin(wav,freqs,amps)
        flux = np.random.normal(flux_clean,noise)
        flux_err = np.ones(npix)*noise
        return wav, flux, flux_err

    wavlim = [1590,1610]
    npix = 2000
    order = 13

    freqs_flat = [13,7,5]
    amps_flat = [0.05,0.07,0.03]
    wav_flat, flux_flat, flux_err_flat = mk_mockdata(freqs_flat,amps_flat,wavlim,npix)
    df_flat = pd.DataFrame(np.array([wav_flat,order*np.ones(npix),flux_flat,flux_err_flat]).T,columns=['wav','order','flux','flux_err'])

    freqs_target = [7000]
    amps_target = [0.04]
    wav_target, flux_target, flux_err_target = mk_mockdata(freqs_target,amps_target,wavlim,npix)
    df_target = pd.DataFrame(np.array([wav_target,order*np.ones(npix),flux_target,flux_err_target]).T,columns=['wav','order','flux_median','flux_err'])

    flux_rmfringe,flux_err_target_rmfringe = remove_fringe_order(df_flat,df_target,order)
    _,_,_,_,_,freqs_rmfringe = ls_periodogram(wav_target,flux_rmfringe,search=[0.05, 0.5])
    diff_abs = np.abs(np.array(freqs_flat) - np.array(freqs_rmfringe))
    assert np.sum(diff_abs)>3

if __name__ == '__main__':
    test_removefringe()
