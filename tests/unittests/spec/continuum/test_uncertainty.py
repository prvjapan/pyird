import pytest
import numpy as np
import pandas as pd
from io import StringIO
from pyird.spec.uncertainty import FluxUncertainty

def test_determine_readout_noise():
    rng = np.random.default_rng(0)
    wavelength = np.arange(900,1050)
    order = np.ones(len(wavelength))
    flux = rng.standard_normal(len(wavelength))
    data = np.array([wavelength, order, flux]).T
    data_str = '\n'.join([' '.join(map(str, row)) for row in data])
    mock_csv = StringIO(data_str)

    # method = 'mad'
    flux_uncertainty = FluxUncertainty()
    flux_uncertainty.combfile = mock_csv
    flux_uncertainty.method = 'mad'
    flux_uncertainty.determine_readout_noise()
    np.testing.assert_array_almost_equal(flux_uncertainty.readout_noise, 1, decimal=1)

    # method = 'std'
    mock_csv = StringIO(data_str)
    flux_uncertainty = FluxUncertainty()
    flux_uncertainty.combfile = mock_csv
    flux_uncertainty.method = 'std'
    flux_uncertainty.determine_readout_noise()
    np.testing.assert_array_almost_equal(flux_uncertainty.readout_noise, 1, decimal=1)

def test_calc_uncertainty():
    wavelength = np.arange(1500, 1600)
    flux = np.ones(len(wavelength))
    continuum = np.ones(len(wavelength))
    tmp_uncertainty = 0.01 * np.ones(len(wavelength))
    data = {'wav': wavelength, 'flux': flux, 'continuum': continuum, 'tmp_uncertainty': tmp_uncertainty}
    df = pd.DataFrame.from_dict(data)

    flux_uncertainty = FluxUncertainty()
    df_continuum = flux_uncertainty.calc_uncertainty(df)

    gain = flux_uncertainty.gain
    readout_noise = flux_uncertainty.readout_noise
    expected_snr = 1*gain/np.sqrt(1*gain + (gain*readout_noise)**2)
    expected_uncertainty = np.sqrt(1/gain + readout_noise**2)/1

    assert df_continuum['sn_ratio'][0] == expected_snr
    assert df_continuum['uncertainty'][0] == expected_uncertainty

def test_calc_uncertainty_overlap_region():
    def make_mock_df(len_df, snr=10):
        wavelength = np.linspace(1500, 1600, len_df)
        flux = np.ones(len(wavelength))
        continuum = np.ones(len(wavelength))
        tmp_uncertainty =  np.ones(len(wavelength))/snr
        sn_ratio = snr * np.ones(len(wavelength))
        data = {'wav': wavelength, 'flux': flux, 'continuum': continuum, 'sn_ratio': sn_ratio, 'uncertainty': tmp_uncertainty}
        df = pd.DataFrame.from_dict(data)
        return df
    len_head = 101
    len_tail = 200
    df_head = make_mock_df(len_head)
    df_tail = make_mock_df(len_tail, snr=20)

    flux_uncertainty = FluxUncertainty()
    sn_ratio, tmp_uncertainty = flux_uncertainty.calc_uncertainty_overlap_region(df_head, df_tail)

    expected_uncertainty = np.sqrt((1 / flux_uncertainty.gain + flux_uncertainty.readout_noise**2)/2)
    np.testing.assert_array_almost_equal(expected_uncertainty, tmp_uncertainty[len_head//2], decimal=1)


if __name__ == '__main__':
    #test_determine_readout_noise()
    #test_calc_uncertainty()
    test_calc_uncertainty_overlap_region()