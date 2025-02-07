import numpy as np
import pytest
from pyird.image.channel import (
    image_to_channel_cube,
    channel_cube_to_image,
    revert_channel_cube,
    eopixel_split,
    eopixel_combine,
    bias_region,
)

@pytest.fixture
def synthetic_image():
    """
    Create a synthetic 2D image for testing.
    """
    np.random.seed(0)
    return np.random.normal(loc=1000, scale=50, size=(1024, 2048))

@pytest.fixture
def synthetic_channel_cube():
    """
    Create a synthetic 3D channel cube for testing.
    """
    np.random.seed(0)
    return np.random.normal(loc=1000, scale=50, size=(32, 32, 2048))

def test_image_to_channel_cube(synthetic_image):
    """
    Test the conversion from an image to a channel cube.
    """
    Nch = 32
    channel_cube = image_to_channel_cube(synthetic_image, Nch=Nch, revert=False)

    # Assertions
    assert channel_cube.shape == (Nch, 1024 // Nch, 2048), \
        "Channel cube shape mismatch."
    assert np.allclose(synthetic_image, channel_cube_to_image(channel_cube, revert=False)), \
        "Image reconstruction from channel cube failed."

def test_channel_cube_to_image(synthetic_channel_cube):
    """
    Test the conversion from a channel cube to an image.
    """
    # Convert the channel cube to an image
    image = channel_cube_to_image(synthetic_channel_cube, revert=False)

    # Assertions
    assert image.shape == (32 * 32, 2048), \
        "Reconstructed image shape mismatch."
    # Verify that the original cube can be retrieved
    reconstructed_cube = image_to_channel_cube(image, revert=False, Nch=32)
    assert np.allclose(reconstructed_cube, synthetic_channel_cube), \
        "Channel cube reconstruction from image failed."

def test_revert_channel_cube(synthetic_channel_cube):
    """
    Test reverting the y-direction of odd channels in the channel cube.
    """
    reverted_cube = revert_channel_cube(np.copy(synthetic_channel_cube))

    # Assertions
    assert reverted_cube.shape == synthetic_channel_cube.shape, \
        "Reverted cube shape mismatch."
    # Verify that odd channels have been reversed
    assert np.allclose(reverted_cube[1::2, :, :], synthetic_channel_cube[1::2, ::-1, :]), \
        "Odd channels were not correctly reverted."

def test_eopixel_split_and_combine(synthetic_channel_cube):
    """
    Test splitting the even and odd pixels and then recombining them.
    """
    eo_tensor = eopixel_split(synthetic_channel_cube)

    # Assertions for split
    assert eo_tensor.shape == (2, 32, 16, 2048), \
        "EO tensor shape mismatch after splitting."
    
    # Combine the split data
    combined_cube = eopixel_combine(eo_tensor)

    # Assertions for combined cube
    assert combined_cube.shape == synthetic_channel_cube.shape, \
        "Combined cube shape mismatch."
    assert np.allclose(combined_cube, synthetic_channel_cube), \
        "EO pixel recombination did not match the original data."

def test_bias_region(synthetic_channel_cube):
    """
    Test extracting the bias region from the channel cube.
    """
    margin = 4
    bias = bias_region(synthetic_channel_cube, margin=margin)

    # Assertions
    assert bias.shape == (32, 32, 2 * margin), \
        "Bias region shape mismatch."
    # Verify the bias region was extracted from both ends
    assert np.allclose(bias[:, :, :margin], synthetic_channel_cube[:, :, :margin]), \
        "Front bias region extraction failed."
    assert np.allclose(bias[:, :, margin:], synthetic_channel_cube[:, :, -margin:]), \
        "Back bias region extraction failed."

if __name__ == '__main__':
    pytest.main([__file__])
