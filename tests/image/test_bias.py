import pytest
from pyird.image.bias import bias_subtract
from pyird.image.channel import image_to_channel_cube, channel_cube_to_image

def test_rmbias():
    import numpy as np
    np.random.seed(1)
    a=np.random.normal(0.0,1.0,(2048,2048))    
    channel_cube=image_to_channel_cube(a)    
    c=bias_subtract(channel_cube)
    image_rmbias=channel_cube_to_image(c)
    assert (np.abs(np.sum(image_rmbias)+10613.724679286175)) < 1.e-8

if __name__=="__main__":
    test_rmbias()
