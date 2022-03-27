import pytest
from pyird.image.aptrace import aptrace

def test_aptrace():
    import numpy as np

    npix = 2048
    step = 100
    nap = npix//step

    dat_test = np.zeros((npix,npix))
    for x in range(0,npix,step):
        dat_test[:,x-2:x+2] = 1000

    cutrow = 1000
    y0, xmin, xmax, coeff = aptrace(dat_test,cutrow,nap)
    assert np.sum(y0) == 20960

if __name__ == '__main__':
    test_aptrace()
