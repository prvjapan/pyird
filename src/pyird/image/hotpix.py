import numpy as np
import pandas as pd
import astropy.io.fits as pyf


def identify_hotpix(im, threshold=10.0):
    """Identification of hotpixels using sep.

    Args:
       im: image
       threshold: threshold of sep

    Returns:
       hotpix_mask: hot pixel mask
       obj: sep object
    """
    from sep import Background, extract, mask_ellipse
    try:
        bkg = Background(im.byteswap().newbyteorder())
    except:
        bkg = Background(im)

    bkg_im = bkg.back()
    im_sub = im-bkg_im
    obj = extract(im_sub, threshold, err=bkg.globalrms)
    hotpix_mask = np.zeros_like(im, dtype=bool)
    mask_ellipse(hotpix_mask, x=obj['x'], y=obj['y'],
                 a=obj['a'], b=obj['b'], theta=obj['theta'])

    return hotpix_mask, obj

def hotpix_fits_to_dat(fitsfile,save_path):
    """Convert .fits to .dat so that it can be used in the RV pipeline.

    Args:
       fitsfile: .fits file of hot pixel mask
       save_path: save path for .dat file

    Returns:
        hot pixel mask in .dat file {pixels, orders, intensity}
    """
    hdu = pyf.open(fitsfile)
    data = hdu[0].data
    wspec = pd.DataFrame([],columns=['wav','order','flux'])
    for i in range(len(data[0,:])):
        pixel = range(1,2049)
        order = np.ones(len(pixel))
        order[:] = i+1
        flux_ord = data[:,i]
        data_order = [pixel,order,flux_ord]
        df_order = pd.DataFrame(data_order,index=['wav','order','flux']).T
        wspec = pd.concat([wspec,df_order])
    wspec = wspec.fillna(0)
    wspec.to_csv(save_path,header=False,index=False,sep=' ')
    return

if __name__ == '__main__':
    import numpy as np
    from pyird.image.bias import bias_subtract_image
    from pyird.utils import irdstream
    from pyird.plot.showmask import show_hotpix
    import astropy.io.fits as pyf
    import pathlib
    mode = 'H'
    datadir = pathlib.Path('/home/kawahara/pyird/data/dark/')
    anadir = pathlib.Path('/home/kawahara/pyird/data/dark/')
    dark = irdstream.Stream2D('targets', datadir, anadir)
    dark.fitsid = [41018]
    if mode == 'YJ':
        dark.fitsid_increment()

    print(dark.rawpath)
    for data in dark.rawpath:
        im = pyf.open(str(data))[0].data
        hd = pyf.open(str(data))[0].header

    im_subbias = bias_subtract_image(im)
    mask, obj = identify_hotpix(im_subbias)
    show_hotpix(obj, im_subbias)
