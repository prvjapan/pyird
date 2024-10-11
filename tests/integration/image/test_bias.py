import pytest
import numpy as np
from pyird.image.bias import bias_subtract
from pyird.image.channel import image_to_channel_cube, channel_cube_to_image


def test_rmbias():
    # Create synthetic data
    np.random.seed(1)
    a = np.random.normal(0.0, 1.0, (2048, 2048))

    # bias subtract in channel cube
    channel_cube = image_to_channel_cube(a)
    c, b = bias_subtract(channel_cube)
    image_rmbias = channel_cube_to_image(c)
    assert np.std(image_rmbias) < np.std(a)


if __name__ == '__main__':
    test_rmbias()
