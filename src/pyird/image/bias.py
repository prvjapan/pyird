""" Bias subtraction
   
    - originally developed by M. Kuzuhara


"""
import numpy as np
import scipy.stats as stats
from pyird.image.channel import bias_region

def bias_subtract(channel_cube):
    """
    Args:
       channel_cube: channel cube

    Returns:
       unbiased_channel_cube: channel cube

    """
    snclip = 3.0 #clipping value
    ref=bias_region(channel_cube)    
    meancp=[]
    for ch_num in range(np.shape(channel_cube)[0]):
        cp=stats.sigmaclip(ref[ch_num,:,:], snclip, snclip)[0]
        meancp.append(np.nanmedian(cp))
    meancp=np.array(meancp)
    unbias_channel_cube=channel_cube-meancp[:,np.newaxis,np.newaxis]
    return unbias_channel_cube

if __name__=="__main__":
    import numpy as np
    from pyird.image.channel import image_to_channel_cube, channel_cube_to_image

    np.random.seed(1)
    a=np.random.normal(0.0,1.0,(2048,2048))    
    channel_cube=image_to_channel_cube(a)    
    c=bias_subtract(channel_cube)
    image_rmbias=channel_cube_to_image(c)

    print(np.abs(np.sum(image_rmbias)-33139.3454845))
