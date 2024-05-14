"""Class for aperture"""

from pyird.image.mask import trace
import warnings

__all__ = ["TraceAperture"]


class TraceAperture(object):
    """aperture instance for trace class"""

    def __init__(self, trace_function, y0, xmin, xmax, coeff, inst):
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
        self.mmf = "mmf12"
        self.inst = inst
        self.info = False
        self.width = None

    def mask(self):
        """mask image

        Returns:
            mask image

        """
        if self.info:
            print("trace_mask used: ", self.mmf)
        return trace(
            self.trace_function,
            self.y0,
            self.xmin,
            self.xmax,
            self.coeff,
            inst=self.inst,
            width=self.width,
        )

    def mmf2(self):
        warn_msg = "Deprecated Use `choose_mmf2_aperture` instead"
        warnings.warn(warn_msg, FutureWarning)
        self.choose_mmf2_aperture()

    def choose_mmf2_aperture(self):
        """choose apertures for mmf2 (star fiber)

        Returns:
            updated variables (y0, xmin, xmax, coeff)

        """
        self.warn_single_mmf()
        self.y0 = self.y0[::2]
        self.xmin = self.xmin[::2]
        self.xmax = self.xmax[::2]
        self.coeff = self.coeff[::2]
        self.mmf = "m2"

    def mmf1(self):
        warn_msg = "Deprecated Use `choose_mmf1_aperture` instead"
        warnings.warn(warn_msg, FutureWarning)
        self.choose_mmf1_aperture()

    def choose_mmf1_aperture(self):
        """choose apertures for mmf1 (comb fiber)

        Returns:
            updated variables (y0, xmin, xmax, coeff)

        """
        self.warn_single_mmf()
        self.y0 = self.y0[1::2]
        self.xmin = self.xmin[1::2]
        self.xmax = self.xmax[1::2]
        self.coeff = self.coeff[1::2]
        self.mmf = "m1"

    def warn_single_mmf(self):
        if len(self.y0) in [21, 51, 52]:
            msg = "Looks a single MMF on the detector. Choosing a mmf aperture may have trouble."
            warnings.warn(msg, UserWarning)
