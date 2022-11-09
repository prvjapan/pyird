import pytest
from pyird.image.hotpix import identify_hotpix,identify_hotpix_sigclip
import numpy as np
import pathlib
import astropy.io.fits as pyf

def test_identify_hotpix():
    basedir = pathlib.Path(__file__).parent.parent.parent
    darkfile = basedir / 'data/dark/IRDA00041018.fits'
    im = pyf.open(str(darkfile))[0].data
    hotpix_mask, obj = identify_hotpix(im)
    assert np.sum(hotpix_mask.ravel())==1153

def test_identify_hotpix_sigclip():
    basedir = pathlib.Path(__file__).parent.parent.parent
    darkfile = basedir / 'data/dark/IRDA00041018.fits'
    im = pyf.open(str(darkfile))[0].data
    hotpix_mask = identify_hotpix_sigclip(im)
    assert np.sum(hotpix_mask.ravel())==16556

if __name__ == '__main__':
    test_identify_hotpix()
    test_identify_hotpix_sigclip()
