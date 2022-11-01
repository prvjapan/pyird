"""reading hotpixel mask."""
import astropy.io.fits as pyf


def read_hotpix(filename):
    """read hotpixel mask file.

    Args:
       path to a hotpixel mask file

    Returns:
       hotpixel mask

    Example:
        >>> # example to read a hotpixel mask used in IRD
        >>> import pkg_resources
        >>> if flat.band=='h':
        >>>     path=pkg_resources.resource_filename('pyird', 'data/hotpix_mask_h_202210_180s.fits')
        >>> elif flat.band=='y':
        >>>     path=pkg_resources.resource_filename('pyird', 'data/hotpix_mask_y_202210_180s.fits')
        >>> hotpix_mask=read_hotpix(path)
    """
    hdu = pyf.open(filename)
    im_hp = hdu[0].data
    if not 'bool' in str(type(im_hp[0][0])):
        im_hp = im_hp.astype(bool)
    return im_hp


if __name__ == '__main__':
    import pkg_resources
    path = (pkg_resources.resource_filename('pyird', 'data/hotpix_mask_h_202210_180s.fits'))
    hotpix_mask = read_hotpix(path)
