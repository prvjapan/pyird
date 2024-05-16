"""Raw Spectral Detector matrix (RSD matrix)"""
import numpy as np
from scipy.signal import medfilt

def multiorder_to_rsd(rawspec, pixcoord, npix=2048, fill_value=np.nan):
    """
    conversion multiorder rawspec+pixcoord to RSD matrix.

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

def rsd_order_medfilt(rsd, kernel_size=9):
    """
    median filtering for spectrum

    Args:
        rsd:    RSD matrix (npix x norder)
        kernel_size:    kernel size for median filter

    Returns:
        median filtered RSD matrix
    """
    rsd_filtered = []
    for i in range(len(rsd[0])):
        filtered_data = medfilt(rsd[:, i], kernel_size=kernel_size)
        rsd_filtered.append(filtered_data)
    rsd_filtered = np.array(rsd_filtered).T
    return rsd_filtered
