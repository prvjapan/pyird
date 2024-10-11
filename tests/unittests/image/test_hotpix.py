import pytest
from pyird.image.hotpix import identify_hotpix, identify_hotpix_sigclip, apply_hotpixel_mask
import numpy as np
import pathlib
import astropy.io.fits as pyf
import importlib

def test_identify_hotpix():
    basedir = pathlib.Path(__file__).parent.parent.parent.parent
    darkfile = basedir / 'data/dark/IRDA00041018.fits'
    im = pyf.open(str(darkfile))[0].data
    hotpix_mask, obj = identify_hotpix(im)
    assert np.sum(hotpix_mask.ravel())==1153

def test_identify_hotpix_sigclip():
    basedir = pathlib.Path(__file__).parent.parent.parent.parent
    darkfile = basedir / 'data/dark/IRDA00041018.fits'
    im = pyf.open(str(darkfile))[0].data
    hotpix_mask = identify_hotpix_sigclip(im)
    assert np.sum(hotpix_mask.ravel())==16556

def test_apply_hotpixel_mask():
    path=importlib.resources.files('pyird').joinpath('data/hotpix_mask_h_202210_180s.fits')

    npix = 2048
    norder = 21
    rsd = np.ones((npix,norder))
    rsd[1,:] = np.nan #including nan
    xmin = np.zeros(norder,dtype=int)
    xmax = np.ones(norder,dtype=int) * npix

    rsd_masked = apply_hotpixel_mask(None, rsd, None, xmin, xmax, None, save_path=path)
    #assert np.nansum(rsd_masked) == (npix-1)*norder  ## nans remains nans
    assert np.nansum(rsd_masked) == (npix)*norder ## when nans are interpolated

if __name__ == '__main__':
    test_identify_hotpix()
    test_identify_hotpix_sigclip()
    test_apply_hotpixel_mask()
