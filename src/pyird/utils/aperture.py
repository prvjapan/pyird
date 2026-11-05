"""Class for aperture"""

from pyird.image.mask import trace
from pyird.image.aperture_width import resolve_aperture_width
import numpy as np
import warnings

__all__ = ["TraceAperture"]


class TraceAperture(object):
    """aperture instance for trace class"""

    def __init__(
        self, trace_function, y0, xmin, xmax, coeff, inst, mask_shape, width=None
    ):
        """initialization

        Args:
            trace_function:  trace function used in interpolation
            y0:  y0
            xmin:  xmin
            xmax:  xmax
            coeff:  polinomial coeff
            inst: instrument (IRD or REACH or IRCS)
            mask_shape: shape of mask
            width: aperture widths ([width_start, width_end])

        """
        self.trace_function = trace_function
        self.y0 = y0
        self.xmin = xmin
        self.xmax = xmax
        self.coeff = coeff
        self.inst = inst
        if inst in ["IRD", "REACH"]:
            self.mmf = "mmf12"
        elif inst in ["IRCS"]:
            self.mmf = "ec" # echelle mode
        else:
            raise ValueError('Error: please set "inst" to either "IRD", "REACH", or "IRCS"')
        self.info = False
        self.width = width
        self.mask_shape = mask_shape

    @property
    def width(self):
        return self._width

    @width.setter
    def width(self, width):
        self._width = resolve_aperture_width(self.inst, width)

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
            mask_shape=self.mask_shape,
        )

    def mask_oneord(self, order=None, axis=-1):
        """mask image for a specific order

        Args:
            order: order number (1,2,...)
            axis: axis along which to identify orders (default: -1)
        """
        mask_all = self.mask()

        prev = np.roll(mask_all, 1, axis=axis)
        idx = [slice(None)] * mask_all.ndim
        idx[axis] = 0
        prev[tuple(idx)] = False

        starts = mask_all & ~prev

        # culculate the block id for each True pixel
        block_id = np.cumsum(starts, axis=axis)

        # set the mask for the specified order
        return mask_all & (block_id == order)

    def choose_aperture(self, fiber):
        """choose apertures of mmf1 (comb fiber) or mmf2 (star fiber)

        Args:
            fiber:  fiber type (mmf1 or mmf2)

        Returns:
            updated variables (y0, xmin, xmax, coeff)

        """
        if len(self.y0) in [21, 51, 52]:
            self.set_single_mmf(fiber)
        else:
            if fiber == "mmf1":
                self.choose_odd_aperture()
            elif fiber == "mmf2":
                self.choose_even_aperture()
            else:
                raise ValueError("fiber should be mmf1 or mmf2")
        print(f"Processing {self.mmf} fiber.")

    def mmf2(self):
        warn_msg = "Deprecated Use `choose_aperture(fiber='mmf2')` instead"
        warnings.warn(warn_msg, FutureWarning)
        self.choose_aperture(fiber="mmf2")

    def choose_mmf2_aperture(self):
        warn_msg = "Deprecated Use `choose_aperture(fiber='mmf2')` instead"
        warnings.warn(warn_msg, FutureWarning)
        self.choose_aperture(fiber="mmf2")

    def choose_even_aperture(self):
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
        warn_msg = "Deprecated Use `choose_aperture(fiber='mmf1')` instead"
        warnings.warn(warn_msg, FutureWarning)
        self.choose_aperture(fiber="mmf1")

    def choose_mmf1_aperture(self):
        warn_msg = "Deprecated Use `choose_aperture(fiber='mmf1')` instead"
        warnings.warn(warn_msg, FutureWarning)
        self.choose_aperture(fiber="mmf1")

    def choose_odd_aperture(self):
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

    def set_single_mmf(self, fiber):
        """set mmf1 or mmf2 for single mmf

        Returns:
            updated variables (mmf) while not reducing the number of apertures

        """
        if fiber == "mmf1":
            self.mmf = "m1"
        elif fiber == "mmf2":
            self.mmf = "m2"
        else:
            raise ValueError("fiber should be mmf1 or mmf2")

    def write_iraf(self, stream, background="none", rot180=False):
        """Write this aperture trace in IRAF APEXTRACT database format.

        The output is written below ``stream.anadir/database``. See
        ``pyird.io.iraf_trace.write_trace_file`` for details.
        """
        from pyird.io.iraf_trace import write_trace_file

        return write_trace_file(
            self,
            stream,
            background=background,
            rot180=rot180,
        )
