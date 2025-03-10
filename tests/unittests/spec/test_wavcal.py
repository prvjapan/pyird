import pytest
import numpy as np
import pandas as pd
from scipy.optimize import leastsq
from scipy.signal import medfilt
import importlib

from pyird.spec.wavcal import (
    pixel_df_to_wav_mat, fit_polynomial, errfunc, fit_wav_solution, calculate_residuals,
    sigmaclip, identify_channel_mode, first_identification, second_identification,
    iterate_fitting
)

@pytest.fixture
def sample_pixwavmap():
    return pd.DataFrame({
        'ORDER': [1, 1, 2, 2],
        'CHANNEL': [10, 15, 5, 25],
        'WAVELENGTH': [5000, 5005, 5010, 5020]
    })


def sample_data(norder=5):
    return np.random.rand(norder, 2048)  # 5 orders and 2048 pixels

@pytest.fixture
def sample_weights():
    return np.ones((5, 2048))

@pytest.fixture
def sample_xy():
    X, Y = np.meshgrid(np.arange(2048), np.arange(5))
    return X, Y

# Test pixel_df_to_wav_mat
def test_pixel_df_to_wav_mat(sample_pixwavmap):
    mat = pixel_df_to_wav_mat(sample_pixwavmap, order_ref=[1, 2], npix=30)
    assert mat.shape == (30, 2)
    assert mat[10, 0] == 5000
    assert mat[15, 0] == 5005
    assert mat[5, 1] == 5010
    assert mat[25, 1] == 5020

# Test fit_polynomial
def test_fit_polynomial(sample_xy):
    X, Y = sample_xy
    params = [1, 0.5, -0.2, 0.1, 1, -0.3]
    result = fit_polynomial((X, Y), Ni=2, Nx=3, params=params)
    assert result.shape == (2048 * 5,)

# Test errfunc
def test_errfunc(sample_xy, sample_weights):
    sample_spec = sample_data()
    X, Y = sample_xy
    params = [1, 0.5, -0.2, 0.1, 1, -0.3]
    residuals = errfunc(params, (X, Y), sample_spec, sample_weights, Ni=2, Nx=3)
    assert len(residuals) == 2048 * 5

# Test fit_wav_solution
def test_fit_wav_solution(sample_xy, sample_weights):
    sample_spec = sample_data()
    X, Y = sample_xy
    coeffs = fit_wav_solution((X, Y), sample_spec, sample_weights, Ni=2, Nx=3)
    assert coeffs.shape == (2, 3)

# Test calculate_residuals
def test_calculate_residuals():
    sample_spec = sample_data()
    wavlength_solution = np.random.rand(5, 2048)
    residuals = calculate_residuals(sample_spec, wavlength_solution)
    assert len(residuals) > 0

# Test sigmaclip
def test_sigmaclip():
    sample_spec = sample_data()
    wavlength_solution = np.random.rand(5, 2048)
    residuals, drop_ind = sigmaclip(sample_spec, wavlength_solution, N=3)
    assert len(residuals) > 0
    assert isinstance(drop_ind, list)

# Test identify_channel_mode
def test_identify_channel_mode_h():
    norder = 21
    orders, channelfile, order_ref = identify_channel_mode(norder)
    assert len(orders) == norder
    assert len(order_ref) == norder

def test_identify_channel_mode_y():
    norder = 51
    orders, channelfile, order_ref = identify_channel_mode(norder)
    assert len(orders) == norder
    assert len(order_ref) == norder

def test_identify_channel_mode_free():
    norder = 20
    channelfile_path = (importlib.resources.files('pyird').joinpath('data/channel_H.list')) 
    orders, channelfile, order_ref = identify_channel_mode(norder, channelfile_path=channelfile_path, ign_ord=[21])
    assert len(orders) == norder
    assert len(order_ref) == norder

# Test first_identification
def test_first_identification():
    sample_spec = sample_data(norder=21)
    channelfile = (importlib.resources.files('pyird').joinpath('data/channel_H.list')) 
    df_pixwavmap_obs = first_identification(sample_spec, channelfile, order_ref=np.arange(1,22))
    assert isinstance(df_pixwavmap_obs, pd.DataFrame)