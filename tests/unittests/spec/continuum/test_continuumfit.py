import pytest
import numpy as np
import pandas as pd
from pyird.spec.continuum_test import ContinuumFit

def test_fit_continuum():
    x = np.arange(10)
    y = - x**2

    continuumfit = ContinuumFit()
    mu, m = continuumfit.fit_continuum(x, y)

    np.testing.assert_array_almost_equal(y, mu, decimal=1)

def test_continuum_rsd():
    pixels = np.arange(2048)
    rsd = np.array([- (pixels-1028)**4/1e12 + 2]*2).T

    continuumfit = ContinuumFit()
    df_continuum = continuumfit.continuum_rsd(rsd)

    # assert by excluding the both ends of order
    np.testing.assert_array_almost_equal(rsd[500:1500,:].T, df_continuum.loc[:,501:1500].values)

def test_continuum_oneord():
    pixels = np.arange(2048)
    wdata = np.array([- (pixels-1028)**4/1e12 + 1]*2).T
    flat = np.array([- (pixels-1028)**4/1e12 + 2]*2).T
    order = np.array([1*np.ones(len(pixels)), 2*np.ones(len(pixels))]).T
    data_target = {'wav': np.array([pixels]*2).T.ravel(), 'flux': wdata.ravel(), 'order': order.ravel()}
    df_wdata = pd.DataFrame.from_dict(data_target)
    data_flat = {'wav': np.array([pixels]*2).T.ravel(), 'flux': flat.ravel(), 'order': order.ravel()}
    df_flat = pd.DataFrame.from_dict(data_flat)
    order_use = 1

    continuumfit = ContinuumFit()
    wav_tmp, flux_tmp, flatflux_tmp, useind, continuum = continuumfit.continuum_oneord(df_wdata, df_flat, order_use)

    # assert by excluding the both ends of order
    np.testing.assert_array_almost_equal(flat[500:1500,order_use].T, continuum[500:1500])


if __name__ == '__main__':
    test_fit_continuum()
    test_continuum_rsd()
    test_continuum_oneord()