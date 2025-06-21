import pytest
from pyird.image.trace_function import trace_legendre
from pyird.utils.aperture import TraceAperture

def test_aperture_real_case():
    import numpy as np
    #The following values are computed by test_aptrace.py
    y0=np.array([98, 198, 298, 398, 498, 598, 698, 798, 898, 998, 1098, 1198, 1298, 1398, 1498, 1598, 1698, 1798, 1898, 1998])
    xmin=np.array([2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])
    xmax=np.array([1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999, 1999])
    coeff=[np.array([ 1.00000000e+00, -4.24474857e-08,  1.43754930e-08]), np.array([ 1.00000000e+00, -1.50951790e-08, -1.47497638e-08]), np.array([1.00000000e+00, 9.85897478e-08, 1.10251354e-07]), np.array([1.00000000e+00, 9.85897478e-08, 1.10251354e-07]), np.array([1.00000000e+00, 9.85897478e-08, 1.10251354e-07]), np.array([ 1.00000000e+00,  1.94637642e-08, -1.63696642e-07]), np.array([ 1.00000000e+00,  1.94637642e-08, -1.63696642e-07]), np.array([ 1.00000000e+00,  1.94637642e-08, -1.63696642e-07]), np.array([ 1.00000000e+00,  1.94637642e-08, -1.63696642e-07]), np.array([ 1.00000000e+00,  1.94637642e-08, -1.63696642e-07]), np.array([ 1.00000000e+00, -6.44062833e-08,  1.82082310e-07]), np.array([ 1.00000000e+00, -6.44062833e-08,  1.82082310e-07]), np.array([ 1.00000000e+00, -6.44062833e-08,  1.82082310e-07]), np.array([ 1.00000000e+00, -6.44062833e-08,  1.82082310e-07]), np.array([ 1.00000000e+00, -6.44062833e-08,  1.82082310e-07]), np.array([ 1.00000000e+00, -6.44062833e-08,  1.82082310e-07]), np.array([ 1.00000000e+00, -6.44062833e-08,  1.82082310e-07]), np.array([ 1.00000000e+00, -6.44062833e-08,  1.82082310e-07]), np.array([ 1.00000000e+00, -6.44062833e-08,  1.82082310e-07]), np.array([ 1.00000000e+00, -6.44062833e-08,  1.82082310e-07])]

    ta = TraceAperture(trace_legendre, y0, xmin, xmax, coeff, inst='IRD', mask_shape=(2048, 2048))
#    ta.trace_function = trace_legendre
    #import matplotlib.pyplot as plt
    #plt.imshow(ta.mask())
    #plt.show()
    assert np.sum(ta.mask()) >= 239760

@pytest.fixture
def trace_aperture():
    return TraceAperture(
        trace_function=None, 
        y0=[i for i in range(100)],  
        xmin=[i for i in range(100)],
        xmax=[i for i in range(100)],
        coeff=[1 for i in range(100)],
        inst="REACH",
        mask_shape=(2048, 2048)
    )

def test_invalid_fiber_type(trace_aperture):
    with pytest.raises(ValueError):
        trace_aperture.choose_aperture("invalid_fiber")

def test_choose_aperture_updates_variables(trace_aperture):
    original_y0 = trace_aperture.y0.copy()
    original_xmin = trace_aperture.xmin.copy()
    original_xmax = trace_aperture.xmax.copy()
    original_coeff = trace_aperture.coeff.copy()
    
    trace_aperture.choose_aperture("mmf1")

    assert trace_aperture.y0 == original_y0[1::2]
    assert trace_aperture.xmin == original_xmin[1::2]
    assert trace_aperture.xmax == original_xmax[1::2]
    assert trace_aperture.coeff == original_coeff[1::2]
    assert trace_aperture.mmf == "m1"

    original_y0 = trace_aperture.y0.copy()
    original_xmin = trace_aperture.xmin.copy()
    original_xmax = trace_aperture.xmax.copy()
    original_coeff = trace_aperture.coeff.copy()

    trace_aperture.choose_aperture("mmf2")

    assert trace_aperture.y0 == original_y0[::2]
    assert trace_aperture.xmin == original_xmin[::2]
    assert trace_aperture.xmax == original_xmax[::2]
    assert trace_aperture.coeff == original_coeff[::2]
    assert trace_aperture.mmf == "m2"


if __name__ == '__main__':
    test_aperture_real_case()    
