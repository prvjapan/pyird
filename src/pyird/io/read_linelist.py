"""reading linelist."""
import numpy as np


def read_linelist(filename):
    """read thar list file.

    Args:
       path to a linelist file

    Returns:
       wavref ndarray

    Example:
        >>> # example to read a Thorium Argon line list used in IRD
        >>> import pkg_resources
        >>> path=(pkg_resources.resource_filename('pyird', "data/thar_ird2.dat"))
        >>> wavref=read_linelist(path)
    """
    wavref = np.loadtxt(filename, comments='#')
    return wavref


if __name__ == '__main__':
    import pkg_resources
    path = (pkg_resources.resource_filename('pyird', 'data/thar_ird2.dat'))
    wavref = read_linelist(path)
