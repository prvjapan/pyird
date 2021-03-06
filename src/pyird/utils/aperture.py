"""Class for aperture"""

import sys
from pyird.image.mask import trace

__all__ = ['TraceAperture']


class TraceAperture(object):
    """aperture instance for trace class

    """
    def __init__(self, trace_function, y0, xmin, xmax, coeff):
        """initialization

        Args:
            trace_function:  trace function used in interpolation
            y0:  y0
            xmin:  xmin
            xmax:  xmax
            coeff:  polinomial coeff

        """
        self.trace_function = trace_function
        self.y0 = y0
        self.xmin = xmin
        self.xmax = xmax
        self.coeff = coeff

    def mask(self):
        """mask image

        Returns:
            mask image

        """
        return trace(self.trace_function, self.y0, self.xmin, self.xmax, self.coeff)

    def mmf2(self):
        """choose apertures for mmf2 (star fiber)

        Returns:
            updated variables (y0, xmin, xmax, coeff)

        """
        self.y0 = self.y0[::2]
        self.xmin = self.xmin[::2]
        self.xmax = self.xmax[::2]
        self.coeff = self.coeff[::2]
