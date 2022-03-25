import pytest
import pkg_resources
from pyird.io.iraf_trace import read_trace_file
from pyird.image.trace_function import trace_legendre
from pyird.image.mask import trace
import numpy as np


def test_read_trace_file():
    path = (pkg_resources.resource_filename('pyird', 'data/samples/aprefA'))
    y0, interp_function, xmin, xmax, coeff = read_trace_file(path)
    assert np.sum(y0) == 13517.2013

def test_mask():
    im=np.zeros((2048,2048))
    path = (pkg_resources.resource_filename('pyird', 'data/samples/aprefA'))
    y0, interp_function, xmin, xmax, coeff = read_trace_file(path)
    unique_interp=np.unique(interp_function)
    if len(unique_interp)==1 and unique_interp[0]==2:
        mask = trace(im, trace_legendre, y0, xmin, xmax, coeff)
    else:
        print("other interpolation function than legendre is not supported yet.")
    import matplotlib.pyplot as plt
    im[mask]=1.0
    assert int(np.sum(im))==246366
    
if __name__ == '__main__':
    test_read_trace_file()
    test_mask()
