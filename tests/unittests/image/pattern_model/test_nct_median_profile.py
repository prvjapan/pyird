import pytest
import numpy as np
from pyird.image.pattern_model import cal_nct, median_profile_xy, pad_nans_with_xy_median
from pyird.image.channel import eopixel_combine, channel_cube_to_image

@pytest.fixture
def test_data():
    channel_median_offset = np.random.random((2, 32))
    xy_profile = np.random.random((32,2048))
    cal_eotensor = xy_profile[np.newaxis, np.newaxis, :, :] + channel_median_offset[:, :, np.newaxis, np.newaxis]
    cal_eotensor_true = cal_eotensor.copy()
    # add outliers
    pixval_low, pixval_high = -100, 100
    Nout = 200
    for i in range(Nout):
        i0 = np.random.randint(0,2); i1 = np.random.randint(0,32); 
        i2 = np.random.randint(0,32); i3 = np.random.randint(0,2048); 
        if (cal_eotensor[i0,i1,i2,i3] == pixval_low) or (cal_eotensor[i0,i1,i2,i3] == pixval_high):
            continue
        else:
            if i<int(Nout/2):
                cal_eotensor[i0,i1,i2,i3] = pixval_low
            else:
                cal_eotensor[i0,i1,i2,i3] = pixval_high
    channel_cube = eopixel_combine(cal_eotensor)
    sample_image = channel_cube_to_image(channel_cube)
    return sample_image, cal_eotensor, cal_eotensor_true, channel_median_offset

@pytest.fixture
def sample_image_nan(test_data):
    sample_image, _, _, _ = test_data
    sample_image_nan = sample_image.copy()
    sample_image_nan[100:110, 100:110] = np.nan  # introduce NaNs for testing
    return sample_image_nan

def test_cal_nct(test_data):
    # Test cal_nct function output shape and type
    sample_image, _, _, _ = test_data
    nct_model = cal_nct(sample_image, margin_npixel=4, Ncor=64, sigma=0.1, xscale=32, yscale=64, cube=False)
    assert isinstance(nct_model, np.ndarray)
    assert nct_model.shape[0] == 2048
    assert nct_model.min() >= 0

def test_cal_nct_with_nan(sample_image_nan):
    # Test cal_nct with input containing NaNs
    nct_model = cal_nct(sample_image_nan, margin_npixel=4, Ncor=64, sigma=0.1, xscale=32, yscale=64, cube=False)
    assert isinstance(nct_model, np.ndarray)
    assert nct_model.shape[0] == 2048


# Tests for median_profile_xy
@pytest.mark.parametrize("channel_axis", [1])
def test_median_profile_xy(test_data, channel_axis):
    _, cal_eotensor_filtered, cal_eotensor_true, channel_median_offset = test_data
    offset_subtracted, xyprofile_offset_subtracted_model = median_profile_xy(
        cal_eotensor_filtered, channel_median_offset, channel_axis
    )

    # Verify that the shapes of the output arrays match the expected shape
    assert offset_subtracted.shape == cal_eotensor_filtered.shape
    assert xyprofile_offset_subtracted_model.shape == cal_eotensor_filtered.shape

    # Verify that the offset is correctly subtracted
    expected_offset_subtracted = cal_eotensor_filtered - channel_median_offset[:, :, np.newaxis, np.newaxis]
    np.testing.assert_almost_equal(offset_subtracted, expected_offset_subtracted, decimal=5)

    # Verify that the median profile is correctly recovered
    expected_eotensor = channel_median_offset[:, :, np.newaxis, np.newaxis] + xyprofile_offset_subtracted_model
    np.testing.assert_almost_equal(cal_eotensor_true, expected_eotensor, decimal=5)


# Tests for pad_nans_with_xy_median
@pytest.mark.parametrize("channel_axis, x_axis, y_axis", [(1, 2, 3)])
def test_pad_nans_with_xy_median(test_data, channel_axis, x_axis, y_axis):
    _, cal_eotensor_filtered, cal_eotensor_true, channel_median_offset = test_data
    # Generate input for pad_nans_with_xy_median
    offset_subtracted, xyprofile_offset_subtracted_model = median_profile_xy(
        cal_eotensor_filtered, channel_median_offset, channel_axis
    )
    image_pattern_model_eotensor = np.copy(xyprofile_offset_subtracted_model)
    image_pattern_model_eotensor[:, 5:-5, :, 1500:1505] = np.nan  # Introducing NaN values for testing

    # Run pad_nans_with_xy_median
    filled_model_eotensor = pad_nans_with_xy_median(
        offset_subtracted, channel_median_offset, image_pattern_model_eotensor, channel_axis, x_axis, y_axis
    )

    # Check that there are no NaNs left in the output array
    assert not np.isnan(filled_model_eotensor).any()

    # Verify that the shape remains consistent
    assert filled_model_eotensor.shape == xyprofile_offset_subtracted_model.shape