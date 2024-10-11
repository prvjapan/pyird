
import importlib
import numpy as np
from pyird.image.mask import trace_from_iraf_trace_file
import pytest

def test_read_trace_file():
    mask_shape = (2048, 2048)
    # image for calibration
    pathlist = (importlib.resources.files('pyird').joinpath('data/samples/aprefA'))
    mask = trace_from_iraf_trace_file(pathlist, mask_shape=mask_shape)
    
    # Assertions
    assert mask.shape == mask_shape, "Mask shape does not match the specified shape."
    assert mask.dtype == bool, "Mask data type should be boolean."

    # Check that the mask is not entirely False
    assert np.any(mask), "Mask should contain some True values."

if __name__ == '__main__':
    pytest.main([__file__])