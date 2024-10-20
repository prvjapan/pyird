import numpy as np
import pytest
from pyird.image.mask import trace

def mock_trace_func(x, y0, xmin, xmax, coeff):
    """
    A mock trace function that simulates tracing output for testing.
    """
    # Create simple linear traces for testing purposes
    return [[y + c[0] * (xi - xmin[i]) for xi in x[i]] for i, (y, c) in enumerate(zip(y0, coeff))]

@pytest.fixture
def trace_parameters():
    """
    Fixture to provide sample trace parameters for testing.
    """
    y0 = [100, 500, 1000]
    xmin = [0, 0, 0]
    xmax = [100, 100, 100]
    coeff = [[0.1], [0.2], [0.3]]  # simple linear coefficients for mock trace function
    return y0, xmin, xmax, coeff

def test_trace_default(trace_parameters):
    """
    Test the trace function with default parameters (IRD instrument).
    """
    y0, xmin, xmax, coeff = trace_parameters
    mask_shape = (2048, 2048)

    # Generate the mask using the mock trace function
    mask = trace(mock_trace_func, y0, xmin, xmax, coeff, mask_shape=mask_shape)

    # Assertions
    assert mask.shape == mask_shape, "Mask shape does not match the specified shape."
    assert mask.dtype == bool, "Mask data type should be boolean."

    # Check that the mask is not entirely False
    assert np.any(mask), "Mask should contain some True values."

def test_trace_reach_instrument(trace_parameters):
    """
    Test the trace function with the REACH instrument and default widths.
    """
    y0, xmin, xmax, coeff = trace_parameters
    mask_shape = (2048, 2048)

    # Generate the mask for REACH instrument
    mask = trace(mock_trace_func, y0, xmin, xmax, coeff, mask_shape=mask_shape, inst='REACH')

    # Assertions
    assert mask.shape == mask_shape, "Mask shape does not match the specified shape."
    assert mask.dtype == bool, "Mask data type should be boolean."

    # Check that the mask is not entirely False
    assert np.any(mask), "Mask should contain some True values."

def test_trace_custom_width(trace_parameters):
    """
    Test the trace function with custom aperture widths.
    """
    y0, xmin, xmax, coeff = trace_parameters
    mask_shape = (2048, 2048)
    custom_width = [3, 5]

    # Generate the mask with custom width settings
    mask = trace(mock_trace_func, y0, xmin, xmax, coeff, mask_shape=mask_shape, width=custom_width)

    # Assertions
    assert mask.shape == mask_shape, "Mask shape does not match the specified shape."
    assert mask.dtype == bool, "Mask data type should be boolean."

    # Check that the mask is not entirely False
    assert np.any(mask), "Mask should contain some True values."

def test_trace_small_mask(trace_parameters):
    """
    Test the trace function with a small mask shape.
    """
    y0, xmin, xmax, coeff = trace_parameters
    mask_shape = (500, 500) ## must not exceed xmax values 

    # Generate the mask with a smaller shape
    mask = trace(mock_trace_func, y0, xmin, xmax, coeff, mask_shape=mask_shape)

    # Assertions
    assert mask.shape == mask_shape, "Mask shape does not match the specified shape."
    assert mask.dtype == bool, "Mask data type should be boolean."

    # Check that the mask is not entirely False
    assert np.any(mask), "Mask should contain some True values."

def test_trace_no_true_values():
    """
    Test the trace function with y0 values outside the mask range,
    expecting an empty (all False) mask.
    """
    y0 = [3000, 3500]  # Values outside typical mask bounds
    xmin = [0, 0]
    xmax = [50, 50]
    coeff = [[0.1], [0.2]]
    mask_shape = (2048, 2048)

    # Generate the mask with the mock trace function
    mask = trace(mock_trace_func, y0, xmin, xmax, coeff, mask_shape=mask_shape)

    # Assertions
    assert mask.shape == mask_shape, "Mask shape does not match the specified shape."
    assert not np.any(mask), "Mask should be entirely False."

if __name__ == '__main__':
    pytest.main([__file__])
