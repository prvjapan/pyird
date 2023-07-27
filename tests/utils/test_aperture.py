import pytest
from pyird.image.aptrace import aptrace
from pyird.image.trace_function import trace_legendre
from pyird.utils.aperture import TraceAperture

def test_aperture():
    import numpy as np
    #The following values are computed by test_aptrace.py
    y0=np.array([98, 198, 298, 398, 498, 598, 698, 798, 898, 998, 1098, 1198, 1298, 1398, 1498, 1598, 1698, 1798, 1898, 1998])
    xmin=np.array([2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])
    xmax=np.array([1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999])
    coeff=[np.array([ 1.00000000e+00, -4.24474857e-08,  1.43754930e-08]), np.array([ 1.00000000e+00, -1.50951790e-08, -1.47497638e-08]), np.array([1.00000000e+00, 9.85897478e-08, 1.10251354e-07]), np.array([1.00000000e+00, 9.85897478e-08, 1.10251354e-07]), np.array([1.00000000e+00, 9.85897478e-08, 1.10251354e-07]), np.array([ 1.00000000e+00,  1.94637642e-08, -1.63696642e-07]), np.array([ 1.00000000e+00,  1.94637642e-08, -1.63696642e-07]), np.array([ 1.00000000e+00,  1.94637642e-08, -1.63696642e-07]), np.array([ 1.00000000e+00,  1.94637642e-08, -1.63696642e-07]), np.array([ 1.00000000e+00,  1.94637642e-08, -1.63696642e-07]), np.array([ 1.00000000e+00, -6.44062833e-08,  1.82082310e-07]), np.array([ 1.00000000e+00, -6.44062833e-08,  1.82082310e-07]), np.array([ 1.00000000e+00, -6.44062833e-08,  1.82082310e-07]), np.array([ 1.00000000e+00, -6.44062833e-08,  1.82082310e-07]), np.array([ 1.00000000e+00, -6.44062833e-08,  1.82082310e-07]), np.array([ 1.00000000e+00, -6.44062833e-08,  1.82082310e-07]), np.array([ 1.00000000e+00, -6.44062833e-08,  1.82082310e-07]), np.array([ 1.00000000e+00, -6.44062833e-08,  1.82082310e-07]), np.array([ 1.00000000e+00, -6.44062833e-08,  1.82082310e-07]), np.array([ 1.00000000e+00, -6.44062833e-08,  1.82082310e-07])]

    ta = TraceAperture(trace_legendre, y0, xmin, xmax, coeff,'IRD')
#    ta.trace_function = trace_legendre
    #import matplotlib.pyplot as plt
    #plt.imshow(ta.mask())
    #plt.show()
    assert np.sum(ta.mask()) >= 239760

if __name__ == '__main__':
    test_aperture()
