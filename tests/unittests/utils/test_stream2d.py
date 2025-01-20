import pytest
from pyird.utils.irdstream import StreamCommon, Stream2D
from pyird.utils.fitsutils import load_fits_data_header
import pathlib
import numpy as np
from unittest.mock import MagicMock, patch

def test_path_1():
    basedir = pathlib.Path(__file__).parent.parent.parent.parent
    datadir = basedir / 'data/dark/'
    anadir = basedir / 'data/dark/'
    s2d = Stream2D('targets', datadir, anadir, fitsid = [41018])
    refs = basedir / 'data/dark/IRDA00041018.fits'
    assert s2d.path()[0]==refs

def test_path_2():
    basedir = pathlib.Path(__file__).parent.parent.parent.parent
    datadir = basedir / 'data/dark/'
    anadir = basedir / 'data/dark/'
    s2d = Stream2D('targets', datadir, anadir)
    s2d.fitsid=[41018]
    refs = basedir / 'data/dark/IRDA00041018.fits'
    assert s2d.path()[0]==refs

def setup_stream2d():
    basedir = pathlib.Path(__file__).parent.parent.parent.parent
    datadir = basedir / 'data/dark/'
    anadir = basedir / 'data/dark/'
    s2d = Stream2D('targets', datadir, anadir, fitsid = [41018])
    return s2d

## StreamCommon
@patch('os.remove')
def test_extclean(mock_remove):
    # Mock the paths
    stream = setup_stream2d()
    mock_path = MagicMock()
    mock_path.exists.return_value = True
    stream.path = MagicMock(return_value=[mock_path, mock_path])

    stream.extclean('.fits')
    assert mock_remove.call_count == 2

def test_check_existence():
    # Mock the paths
    stream = setup_stream2d()
    mock_path_exists = MagicMock()
    mock_path_exists.exists.return_value = True

    mock_path_not_exists = MagicMock()
    mock_path_not_exists.exists.return_value = False

    stream.extpath = MagicMock(side_effect=[
        [mock_path_not_exists, mock_path_exists],
        [mock_path_not_exists, mock_path_exists]
    ])

    ext_noexist, extf_noexist = stream.check_existence('in', 'out')
    assert len(ext_noexist) == 1
    assert len(extf_noexist) == 1

def test_remove_bias():
    stream = setup_stream2d()
    stream.rawpath = 'test_rawpath'
    stream.remove_bias()
    assert stream.extension == '_rb'

def test_print_if_info_is_true(capsys):
    stream = setup_stream2d()
    stream.print_if_info_is_true('test message')
    captured = capsys.readouterr()
    assert 'test message' in captured.out

## Stream2D
def test_detector_handling_rotate():
    stream2d = setup_stream2d()
    stream2d.rotate = True
    img, _ = load_fits_data_header(stream2d.rawpath[0])
    rotated_img = stream2d.detector_handling(img, mode='load')
    assert rotated_img[0][0] == img[0][int(img.shape[0])-1]
    assert rotated_img[int(img.shape[0])-1][int(img.shape[0])-1] == img[int(img.shape[0])-1][0]

def test_detector_handling_inverse():
    stream2d = setup_stream2d()
    stream2d.inverse = True
    img, _ = load_fits_data_header(stream2d.rawpath[0])
    inversed_img = stream2d.detector_handling(img, mode='load')
    assert inversed_img[0][0] == img[int(img.shape[0])-1][0]
    assert inversed_img[int(img.shape[0])-1][int(img.shape[0])-1] == img[0][int(img.shape[0])-1]

def test_immedian():
    stream2d = setup_stream2d()
    median_image = stream2d.immedian()
    img, _ = load_fits_data_header(stream2d.rawpath[0])
    np.testing.assert_array_equal(median_image, img)

def test_write_df_spec_wav(tmp_path):
    stream2d = setup_stream2d()
    spec_m2 = np.array([[1, 2], [3, 4]])
    reference = np.array([[10, 20], [30, 40]])
    save_path = tmp_path / "test.csv"
    wspec = stream2d.write_df_spec_wav(spec_m2, reference, save_path)
    assert save_path.exists()
    assert len(wspec) == 4
    assert list(wspec.columns) == ["wav", "order", "flux"]

def test_fitsid_increment():
    stream = setup_stream2d()
    stream.not_ignore_warning = False 
    val = stream.fitsid[0]
    stream.fitsid_increment()
    assert stream.fitsid[0] == val+1

def test_fitsid_decrement():
    stream = setup_stream2d()
    stream.not_ignore_warning = False 
    val = stream.fitsid[0]
    stream.fitsid_increment()
    stream.fitsid_decrement()
    assert stream.fitsid[0] == val
    
def test_fitsid_decrement_value_error_for_no_incremented():
    stream = setup_stream2d()
    stream.not_ignore_warning = True
    with pytest.raises(ValueError):
        stream.fitsid_decrement()

def test_fitsid_increment_value_error_for_double_increment():
    stream = setup_stream2d()
    stream.not_ignore_warning = False
    stream.fitsid_increment()
    with pytest.raises(ValueError):
        stream.fitsid_increment()

def setup_stream2d_band(band):
    basedir = pathlib.Path(__file__).parent.parent.parent.parent
    datadir = basedir / 'data/dark/'
    anadir = basedir / 'data/dark/'
    s2d = Stream2D('targets', datadir, anadir, fitsid = [41018], band=band, not_ignore_warning=False)
    return s2d


def test_fitsid_sets_band_at_stream2d_h():
    stream = setup_stream2d_band("h")
    assert stream.fitsid[0] == 41019
    assert stream.band == "h"

def test_fitsid_sets_band_at_stream2d_y():
    stream = setup_stream2d_band("y")
    assert stream.fitsid[0] == 41018
    assert stream.band == "y"

def test_fitsid_sets_band_at_stream2d_h_increment_error():
    stream = setup_stream2d_band("h")
    with pytest.raises(ValueError):
        stream.fitsid_increment()
    

    
        
if __name__ == '__main__':
    test_path_1()
    test_path_2()
    