import pytest
import importlib
from pyird.io.iraf_trace import read_trace_file
import numpy as np

def test_read_trace_file():
    path = (importlib.resources.files('pyird').joinpath('data/samples/aprefA'))
    y0, interp_function, xmin, xmax, coeff = read_trace_file(path)
    assert np.sum(y0) == 13517.2013

if __name__ == '__main__':
    test_read_trace_file()
