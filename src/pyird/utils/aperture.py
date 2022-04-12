"""Class for aperture"""

import sys
from pyird.image.mask import trace

__all__ = ['TraceAperture']


class TraceAperture(object):
    """aperture instance for trace class

    """    
    def __init__(self, trace_funcion, y0, xmin, xmax, coeff):
        """initialization
        
        Args:
            trace_funcion:  trace funcion used in interpolation
            y0:  y0
            xmin:  xmin
            xmax:  xmax
            coeff:  polinomial coeff
        
        """
        self.trace_funcion = trace_funcion
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

    
