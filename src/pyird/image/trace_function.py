import numpy as np


def trace_legendre(x, y0, xmin, xmax, coeff):
    """trace Legendre function.

    Args:
       x: x-array
       y0: y-offset
       xmin: xmin
       xmax: xmax
       coeff: Legendre polynomial coefficients
    """
    from numpy.polynomial.legendre import legval
    norder = len(y0)
    trace_lines = []
    for i in range(0, norder):
        x_ = np.array(x[i])
        xmax_ = xmax[i]
        xmin_ = xmin[i]
        nx = (2.*x_ - (xmax_+xmin_))/(xmax_-xmin_)
        f = legval(nx, coeff[i])+y0[i]-1
        trace_lines.append(f)

    return trace_lines


if __name__ == '__main__':
    import pkg_resources
    from pyird.io.iraf_trace import read_trace_file
    import numpy as np

    pathC = (pkg_resources.resource_filename('pyird', 'data/samples/aprefC'))
    path_c = (pkg_resources.resource_filename('pyird', 'data/samples/apref_c'))
    y0, interp_function, xmin, xmax, coeff = read_trace_file([pathC, path_c])
