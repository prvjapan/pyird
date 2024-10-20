import numpy as np
import pytest
import pandas as pd
from pyird.image.aptrace import (
    cross_section,
    set_aperture,
    trace_pix,
    fitfunc,
    errfunc,
    fit_ord,
)

@pytest.fixture
def synthetic_data():
    """
    Create a synthetic 2D array to simulate flat data for testing.
    """
    npix = 2048
    step = 100
    nap = npix//step
    data = np.zeros((npix,npix))

    # Add vertical lines at regular intervals
    for x in range(0,npix,step):
        data[:,x-2:x+2] = 1000
    return data

def test_cross_section(synthetic_data):
    """
    Test the cross_section function to ensure it detects peaks correctly.
    """
    nrow = 1024
    nap = 20
    masked_data, peak_indices = cross_section(synthetic_data, nrow, nap)
    
    # Assert the shape of the output
    assert len(masked_data) == synthetic_data.shape[1], "Masked data length mismatch."
    assert len(peak_indices) == nap, f"Expected {nap} peaks, found {len(peak_indices)}."
    
    # Check if peaks are roughly in the expected positions
    assert all(90 <= idx <= 2010 for idx in peak_indices), "Peak index out of expected range."

def test_set_aperture(synthetic_data):
    """
    Test the set_aperture function to verify correct peak detection at a given row.
    """
    cutrow = 1024
    nap = 20
    peak_indices, final_cutrow = set_aperture(synthetic_data, cutrow, nap, plot=False)
    
    # Assertions
    assert len(peak_indices) == nap, f"Expected {nap} peaks, found {len(peak_indices)}."

def test_trace_pix(synthetic_data):
    """
    Test the trace_pix function to check aperture tracing along rows.
    """
    # Set up initial conditions
    cutrow = 1024
    nap = 20
    peak_indices, _ = set_aperture(synthetic_data, cutrow, nap, plot=False)
    peak_ind = peak_indices[0]
    
    # Trace the pixels
    x_ord, y_ord, y0_ord = trace_pix(synthetic_data, cutrow, peak_ind)
    
    # Assertions
    assert len(x_ord) > 0, "No rows traced."
    assert len(y_ord) > 0, "No columns traced."
    assert y0_ord is not None, "y0_ord is None."

def test_fitfunc():
    """
    Test the fitfunc to verify Legendre polynomial fitting.
    """
    x = np.linspace(0, 2048, 50)
    y0 = 1000
    coeff = [1, 0.1, 0.01]
    
    # Generate fitted values
    fitted_values = fitfunc(x, y0, coeff)
    
    # Assertions
    assert len(fitted_values) == len(x), "Fitted values length mismatch."
    assert np.all(np.isfinite(fitted_values)), "Fitted values contain non-finite numbers."

def test_errfunc():
    """
    Test the errfunc to check residual calculation.
    """
    x = np.linspace(0, 2048, 50)
    y0 = 1000
    data = np.sin(x) + y0
    coeff = [1, 0.1, 0.01]
    
    # Calculate residuals
    residuals = errfunc(coeff, x, y0, data)
    
    # Assertions
    assert len(residuals) == len(x), "Residuals length mismatch."
    assert np.all(np.isfinite(residuals)), "Residuals contain non-finite numbers."

def test_fit_ord():
    """
    Test the fit_ord function to verify least-square fitting.
    """
    x = np.linspace(0, 2048, 50)
    y0 = 1000
    data = fitfunc(x, y0, [1, 0.1, 0.01]) + np.random.normal(0, 0.1, len(x))
    
    # Perform fitting
    best_fit_params = fit_ord(x, y0, data)
    
    # Assertions
    assert len(best_fit_params) == 3, "Expected 3 fit parameters."
    assert np.all(np.isfinite(best_fit_params)), "Fit parameters contain non-finite numbers."

if __name__ == '__main__':
    pytest.main([__file__])
