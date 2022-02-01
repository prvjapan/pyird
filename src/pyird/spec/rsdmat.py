"""Raw Spectral Detector matrix (RSD matrix)"""
import numpy as np


def multiorder_to_rsd(rawspec, pixcoord, npix=2048, fill_value=np.nan):
    """conversion multiorder rawspec+pixcoord to RSD matrix.

    Args:
        rawspec: multiorder rawspec
        pixcoord: multiorder pixel coordinate
        npix: number of detector pixels in y direction
        fill_value: filled value in empty elements

    Returns:
        RSD matrix (npix x norder)
    """
    norder = len(rawspec)
    assert norder == len(pixcoord)
    rsd = np.ones((npix, norder))*fill_value
    for i in range(0, len(pixcoord)):
        rsd[pixcoord[i], i] = rawspec[i]
    return rsd
