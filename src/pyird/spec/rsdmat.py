"""Raw Spectral Detector matrix (RSD matrix)


"""

def multiorder_to_rsd(rawspec, pixcoord, npix=2048):
    """conversion multiorder rawspec+pixcoord to RSD matrix

    Args: 
        rawspec: multiorder rawspec
        pixcoord: multiorder pixel coordinate
        npix: number of detector pixels in y direction

    Returns:
        RSD matrix (npix x norder)

    """
    norder=len(rawspec)
    assert norder == len(pixcoord)
    rsd=np.zeros((npix,norder))

