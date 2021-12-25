import matplotlib.pyplot as plt

def identify_hotpix(im):
    """Identification of hotpixels using sep

    Args:
       im: image

    Returns:
       hotpix_mask: hot pixel mask
       obj: sep object

    """
    from sep import Background,extract,mask_ellipse
    try:
        bkg = Background(im.byteswap().newbyteorder())
    except:
        bkg = Background(im)

    bkg_im = bkg.back()
    im_sub=im-bkg_im
    obj = extract(im_sub, 1.5, err=bkg.globalrms)
    hotpix_mask = np.zeros_like(im, dtype=bool)
    mask_ellipse(hotpix_mask,x=obj['x'],y=obj['y'],a=obj['a'],b=obj['b'],theta=obj['theta'])
    
    return hotpix_mask,obj

if __name__=="__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    from pyird.image.bias import bias_subtract_image
    from pyird.utils import fitsset,irdstream
    from pyird.plot.showmask import show_hotpix
    import astropy.io.fits as pyf
    import pathlib

    
    mode="YJ"
    datadir=pathlib.Path("/home/kawahara/pyird/data/dark/")
    anadir=pathlib.Path("/home/kawahara/pyird/data/dark/") 
    dark=irdstream.Stream2D("targets",datadir,anadir)
    dark.fitsid=[41018]    
    if mode=="H":
        dark.fitsid_increment()        

    print(dark.rawpath)
    for data in dark.rawpath:
        im = pyf.open(str(data))[0].data
        hd = pyf.open(str(data))[0].header

    im_subbias=bias_subtract_image(im)
    mask,obj=identify_hotpix(im_subbias)
    show_hotpix(obj,mask)
    

