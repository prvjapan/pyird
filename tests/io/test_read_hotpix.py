import pytest
import pkg_resources
from pyird.io.read_hotpix import read_hotpix
import numpy as np


def test_read_hotpix():
    path = (pkg_resources.resource_filename('pyird', 'data/hotpix_mask_h_202210_180s.fits'))
    hotpix_mask = read_hotpix(path)
    assert np.sum(hotpix_mask.ravel()) == 14579


if __name__ == '__main__':
    test_read_hotpix()
