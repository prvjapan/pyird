import numpy as np
import pytest
from pyird.image.channel import image_to_channel_cube, channel_cube_to_image
from pyird.image.bias import bias_subtract, bias_subtract_image

@pytest.fixture
def synthetic_image():
    """
    Creates a synthetic 2D image for testing.
    The image contains different pixel values to simulate a bias pattern.
    """
    np.random.seed(0)  # For reproducibility
    image = np.random.normal(loc=1000, scale=50, size=(2048, 2048))
    return image

def test_bias_subtract():
    """
    Test the bias_subtract function using a synthetic channel cube.
    """
    # Create a synthetic channel cube
    np.random.seed(0)
    channel_cube = np.random.normal(loc=1000, scale=50, size=(4, 512, 512))

    # Run bias subtraction
    bias_subtracted_channel_cube, channel_bias = bias_subtract(channel_cube)

    # Assertions to check output shapes
    assert bias_subtracted_channel_cube.shape == channel_cube.shape, \
        "The shape of the bias-subtracted cube should match the input cube."
    assert channel_bias.shape == (4,), \
        "The shape of the channel bias array should match the number of channels."

    # Check if the mean bias is subtracted from the original data
    for ch_num in range(4):
        assert np.isclose(np.mean(channel_cube[ch_num] - channel_bias[ch_num]), 
                          np.mean(bias_subtracted_channel_cube[ch_num]), atol=1e-1), \
            f"Channel {ch_num} bias subtraction failed."

def test_bias_subtract_image(synthetic_image):
    """
    Test the bias_subtract_image function using a synthetic image.
    """
    # Run bias subtraction on the image
    bias_subtracted_image = bias_subtract_image(synthetic_image)

    # Assertions to check output shape
    assert bias_subtracted_image.shape == synthetic_image.shape, \
        "The shape of the bias-subtracted image should match the input image."

    # Verify that the bias is indeed subtracted
    original_mean = np.mean(synthetic_image)
    corrected_mean = np.mean(bias_subtracted_image)
    assert corrected_mean < original_mean, \
        "The mean of the bias-subtracted image should be less than the original mean."

if __name__ == '__main__':
    pytest.main([__file__])
