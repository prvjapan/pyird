import pytest
import pkg_resources
from pyird.io.iraf_trace import read_trace_file
from pyird.image.trace_function import trace_legendre
from pyird.image.mask import trace_from_iraf_trace_file
import numpy as np


def test_read_trace_file():
    path = (pkg_resources.resource_filename('pyird', 'data/samples/aprefA'))
    y0, interp_function, xmin, xmax, coeff = read_trace_file(path)
    assert np.sum(y0) == 13517.2013

def test_mask_from_trace_file():
    path = (pkg_resources.resource_filename('pyird', 'data/samples/aprefA'))
    im=np.zeros((2048,2048))
    mask=trace_from_iraf_trace_file(path)
    im[mask]=1.0
    assert int(np.sum(im))==246366
    
if __name__ == '__main__':
    test_read_trace_file()
    test_mask_from_trace_file()
