import matplotlib.pyplot as plt
def identify_hotpix(im):
    """Identification of hotpixels

    Args:
       im: image

    """
    import sep
    try:
        bkg = sep.Background(im.byteswap().newbyteorder())
    except:
        bkg = sep.Background(im)

    bkg_im = bkg.back()
    im_sub=im-bkg_im
    obj = sep.extract(im_sub, 1.5, err=bkg.globalrms)
    mask = np.zeros_like(im, dtype=bool)
    sep.mask_ellipse(mask,x=obj['x'],y=obj['y'],a=obj['a'],b=obj['b'],theta=obj['theta'])
    from pyird.plot.showmask import show_hotpix
    show_hotpix(obj,mask)
    
if __name__=="__main__":
    import numpy as np
    from pyird.image.bias import bias_subtract_image
    from pyird.image.channel import image_to_channel_cube, channel_cube_to_image
    from pyird.utils import fitsset,irdstream
    import astropy.io.fits as pyf
    import tqdm
    import pathlib

    
    mode="YJ"
    datadir=pathlib.Path("/home/kawahara/pyird/data/dark/")
    anadir=pathlib.Path("/home/kawahara/pyird/data/dark/") 
    dark=irdstream.Stream2D("targets",datadir,anadir)
    dark.fitsid=[41018]    
    if mode=="H":
        dark.fitsid_increment()        

    print(dark.rawpath)
    for data in tqdm.tqdm(dark.rawpath):
        im = pyf.open(str(data))[0].data
        hd = pyf.open(str(data))[0].header

    im_subbias=bias_subtract_image(im)
    identify_hotpix(im_subbias)
    
    import matplotlib.pyplot as plt

