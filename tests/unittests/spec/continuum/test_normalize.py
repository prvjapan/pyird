import pytest
import numpy as np
import pandas as pd
from pyird.spec.normalize import SpectrumNormalizer

def test_determine_scale_continuum():
    expected_scale = 1/5
    pixels = np.arange(2048)
    wdata = np.array([- (pixels-1028)**4/1e12 + 1]*2).T
    flat = np.array([- (pixels-1028)**4/1e12 + 1/expected_scale]*2).T
    order = np.array([1*np.ones(len(pixels)), 2*np.ones(len(pixels))]).T
    data_target = {'wav': np.array([pixels]*2).T.ravel(), 'flux': wdata.ravel(), 'order': order.ravel()}
    df_wdata = pd.DataFrame.from_dict(data_target)
    data_flat = {'wav': np.array([pixels]*2).T.ravel(), 'flux': flat.ravel(), 'order': order.ravel()}
    df_flat = pd.DataFrame.from_dict(data_flat)
    standard_order = 1

    spectrum_normalizer = SpectrumNormalizer()
    scale = spectrum_normalizer.determine_scale_continuum(df_wdata, df_flat, standard_order)
    np.testing.assert_array_almost_equal(expected_scale, scale, decimal=1)

def test_trim_nonzero_flux():
    len_zero = 3
    y = [0]*len_zero + list(np.arange(1000)) + [0]*len_zero
    df = pd.DataFrame(y,columns=["flux"])

    spectrum_normalizer = SpectrumNormalizer()
    interp_nonzero_beginning, interp_nonzero_ending = spectrum_normalizer.interp_nonzero_ind
    df_trim = spectrum_normalizer.trim_nonzero_flux(df)

    assert df_trim.index.min() == len_zero + interp_nonzero_beginning + 1
    assert df_trim.index.max() == len(y) - (len_zero+1)  - interp_nonzero_ending - 1

if __name__ == '__main__':
    test_determine_scale_continuum()
    test_trim_nonzero_flux()