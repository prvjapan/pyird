import pytest
from pyird.image.hotpix import identify_hotpix, identify_hotpix_sigclip, apply_hotpixel_orderedge_mask
import numpy as np
import pathlib
import astropy.io.fits as pyf

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

@pytest.mark.parametrize("hotpix_mode, expected_nan", [("nan", True), ("interp", False)])
def test_apply_hotpixel_orderedge_mask_hotpix_handling(tmp_path, hotpix_mode, expected_nan):
    npix = 8
    norder = 1
    rsd = np.arange(npix, dtype=float).reshape(npix, norder)
    rsd[3, 0] = np.nan
    rsd[5, 0] = 100.0

    hotpix_mask_1d = np.zeros_like(rsd)
    hotpix_mask_1d[5, 0] = 1
    mask_path = tmp_path / "hotpix_mask_1d.fits"
    pyf.writeto(mask_path, hotpix_mask_1d, overwrite=True)

    xmin = np.array([1])
    xmax = np.array([6])

    rsd_masked = apply_hotpixel_orderedge_mask(
        np.ones_like(rsd),
        rsd,
        None,
        xmin,
        xmax,
        None,
        hotpix_mode=hotpix_mode,
        save_path=str(mask_path)
    )

    if expected_nan:
        assert np.isnan(rsd_masked[3, 0])
        assert np.isnan(rsd_masked[5, 0])
    else:
        assert rsd_masked[3, 0] == pytest.approx(3.0)
        assert rsd_masked[5, 0] == pytest.approx(5.0)

    assert rsd_masked[0, 0] == 0
    assert rsd_masked[7, 0] == 0


def test_apply_hotpixel_orderedge_mask_only_edges():
    npix = 6
    rsd = np.column_stack((np.arange(npix, dtype=float), np.arange(npix, dtype=float) + 10))
    rsd[2, 0] = np.nan

    xmin = np.array([1, 0])
    xmax = np.array([3, 4])

    rsd_masked = apply_hotpixel_orderedge_mask(None, rsd, None, xmin, xmax, None)

    assert rsd_masked[0, 0] == 0
    assert rsd_masked[4, 0] == 0
    assert rsd_masked[5, 0] == 0
    assert rsd_masked[5, 1] == 0
    assert rsd_masked[1, 0] == rsd[1, 0]
    assert np.isnan(rsd_masked[2, 0])
    assert rsd_masked[3, 0] == rsd[3, 0]
    assert np.all(rsd_masked[:5, 1] == rsd[:5, 1])

if __name__ == '__main__':
    test_identify_hotpix()
    test_identify_hotpix_sigclip()
    test_apply_hotpixel_orderedge_mask_hotpix_handling()
    test_apply_hotpixel_orderedge_mask_only_edges()
