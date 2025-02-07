import numpy as np
import pandas as pd
import pytest
from pyird.image.oned_extract import flatten, sum_weighted_apertures

def mock_trace_func(x, y0, xmin, xmax, coeff):
    """
    A mock trace function that simulates the output of tracing for testing.
    Returns a simple linear trace based on the provided coefficients.
    """
    return [np.array([y + c[0] * (xi - xmin[i]) for xi in x[i]]) for i, (y, c) in enumerate(zip(y0, coeff))]

@pytest.fixture
def synthetic_image():
    """
    Create a synthetic 2D image for testing.
    """
    np.random.seed(0)
    return np.random.normal(loc=1000, scale=50, size=(2048, 2048))

@pytest.fixture
def trace_parameters():
    """
    Fixture to provide sample trace parameters for testing.
    """
    y0 = [100, 500, 1000]
    xmin = [0, 0, 0]
    xmax = [50, 50, 50]
    coeff = [[0.1], [0.2], [0.3]]  # simple linear coefficients for mock trace function
    return y0, xmin, xmax, coeff

def test_flatten_default(synthetic_image, trace_parameters):
    """
    Test the flatten function with default parameters.
    """
    y0, xmin, xmax, coeff = trace_parameters

    # Generate spectra using the flatten function
    spec, pixcoord, rotim, tl, iys_all, iye_all = flatten(
        synthetic_image, mock_trace_func, y0, xmin, xmax, coeff, onepix=False
    )

    # Assertions
    assert len(spec) == len(y0), "The number of spectra does not match the number of orders."
    assert len(pixcoord) == len(y0), "The number of pixel coordinate arrays does not match the number of orders."
    assert rotim.shape == synthetic_image.shape, "Rotated image shape mismatch."
    assert np.all(np.isfinite(tl)), "Trace line data contains non-finite values."

def test_flatten_onepix(synthetic_image, trace_parameters):
    """
    Test the flatten function with onepix=True to ensure pixel-by-pixel extraction works.
    """
    y0, xmin, xmax, coeff = trace_parameters
    width = [2, 4]

    # Generate pixel-by-pixel spectra using the flatten function
    df_onepix = flatten(
        synthetic_image, mock_trace_func, y0, xmin, xmax, coeff, onepix=True, width=width
    )

    # Assertions
    assert isinstance(df_onepix, dict), "Expected output to be a dictionary of DataFrames."
    assert len(df_onepix) == (width[1] + width[0]), "The number of aperture widths does not match."

def synthetic_df_flatn(width):
    """
    Create a synthetic dictionary for df_flatn with keys of f'ec{width[i]}'.
    """
    orders = range(1, 4)  # Assuming there are 3 orders for the test
    pixels = range(1, 2049)
    df_flatn = {}

    # Generate DataFrames for each width entry
    for i in range(-width[0], width[1]):
        df_flatn[f'ec{i}'] = pd.DataFrame(
            np.random.rand(len(orders), len(pixels)),
            index=orders,
            columns=pixels
        )
    
    return df_flatn

def test_sum_weighted_apertures(synthetic_image, trace_parameters):
    """
    Test the sum_weighted_apertures function to ensure it computes the weighted sum correctly.
    """
    y0, xmin, xmax, coeff = trace_parameters
    width = [2, 4]
    inst = "IRD"

    df_flatn = synthetic_df_flatn(width)
    
    # Generate the summed weighted apertures
    df_sum_wap = sum_weighted_apertures(
        synthetic_image, df_flatn, y0, xmin, xmax, coeff, width, inst
    )

    # Assertions
    assert isinstance(df_sum_wap, pd.DataFrame), "Expected output to be a DataFrame."
    assert df_sum_wap.shape == next(iter(df_flatn.values())).shape, \
        "The output shape does not match the input DataFrame shape."

if __name__ == '__main__':
    pytest.main([__file__])
