import pytest
import numpy as np
from pyird.image.aptrace import aptrace

def test_aptrace():
    # Create synthetic data
    npix = 2048
    step = 100
    nap = npix//step
    dat_test = np.zeros((npix,npix))

    # Add vertical lines at regular intervals
    for x in range(0,npix,step):
        dat_test[:,x-2:x+2] = 1000

    # Run aptrace and get output parameters
    cutrow = 1000
    y0, xmin, xmax, coeff = aptrace(dat_test,cutrow,nap,plot=False)

    # Perform assertions
    assert y0 is not None, "y0 output is None"
    assert xmin is not None, "xmin output is None"
    assert xmax is not None, "xmax output is None"
    assert coeff is not None, "coeff output is None"
    assert np.sum(y0) == 20960, f"Sum of y0 is {np.sum(y0)}, expected 20960"


if __name__ == '__main__':
    test_aptrace()
