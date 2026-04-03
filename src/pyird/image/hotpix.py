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

    Notes:
        See also issue #121 for a comparison among `identify_hotpix`, `identify_hotpix_sigclip`, and another method.
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

def identify_hotpix_sigclip(im,sigma=4,maxiters=5,frac=0.005):
    """Identification of hotpixels using sigma clip.

    Args:
       im: image
       sigma: maximum value of sigma exploration range
       maxiters: maximum value of iterations
       frac: mask percentage target value (empirically 0.5%)

    Returns:
       hotpix_mask: hot pixel mask

    Notes:
        See also issue #121 for a comparison among `identify_hotpix`, `identify_hotpix_sigclip`, and another method.
    """
    from astropy.stats import SigmaClip
    import itertools
    import functools

    def sigmaclip(im,*params):
        sigma, count = params[0][0], params[0][1]
        sigclip = SigmaClip(sigma=sigma,maxiters=count)
        filtered_data = sigclip(im)
        hotpix_mask = filtered_data.mask
        frac_masked = len(im[hotpix_mask].ravel())/len(im.ravel())
        return hotpix_mask,frac_masked

    ## search for best parameters
    sigmas = np.arange(1,sigma+1,1)
    iters = np.arange(1,maxiters+1,1)
    si_list = list(itertools.product(sigmas,iters))
    results = list(map(functools.partial(sigmaclip,im),si_list))
    results = np.array(results,dtype='object')
    ## hotpix mask
    ind_frac = np.argmin(np.abs(results[:,1]-frac))
    si_use = si_list[ind_frac]
    hotpix_mask,frac_masked = sigmaclip(im,si_use)
    print(f"{frac_masked * 100:.2f}% of pixels flagged in hot-pixel mask.")

    return hotpix_mask

def apply_hotpixel_orderedge_mask(
        hotpix_mask, 
        rsd, 
        y0, xmin, xmax, coeff, 
        hotpix_mode="nan", 
        save_path=None
        ):
    """ mask hotpixel and/or order edge. 

    Args:
        hotpix_mask: mask made from dark. if None, only apply order edge mask.
        rsd: extracted spectrum to be masked
        y0, xmin, xmax, coeff: trace infomation
        hotpix_mode: raplace bad pixels by np.nan if hotpix_mode="nan", while replace them based on spline interpolation if hotpix_mode="interp". 
        save_path: path to save hotpixel mask

    Returns:
        masked and/or interpolated spectrum

    """
    from scipy.interpolate import InterpolatedUnivariateSpline as IUS
    from pyird.image.oned_extract import flatten
    from pyird.image.trace_function import trace_legendre
    from pyird.spec.rsdmat import multiorder_to_rsd

    if hotpix_mask is not None:
        try:
            hdu = pyf.open(save_path)[0]
            hotpix1D = hdu.data
        except:
            mask = hotpix_mask.astype(int)
            rawspec, pixcoord, _, _, _, _ = flatten(mask, trace_legendre, y0, xmin, xmax, coeff)
            hotpix1D = multiorder_to_rsd(rawspec, pixcoord)
            hdux = pyf.PrimaryHDU(hotpix1D)
            hdulist = pyf.HDUList([hdux])
            if save_path is not None:
                hdulist.writeto(save_path, overwrite=True)
        mask_ind = hotpix1D>0

    flux = rsd
    pixels = np.array([np.arange(0,flux.shape[0],1)]*flux.shape[1]).T
    flux_new = []
    for j in range(flux.shape[1]):
        mask_edge = np.ones(len(flux[:,j]),dtype=bool)
        mask_edge[xmin[j]:xmax[j]+1] = False
        mask_nan = np.isnan(flux[:,j])
        flux_tmp = flux[:,j].copy()
        if hotpix_mask is not None:
            mask_ord = mask_ind[:,j]
            pix_masked = pixels[:,j][~(mask_ord | mask_edge | mask_nan)]
            flux_masked = flux[:,j][~(mask_ord | mask_edge | mask_nan)]

            ### raplace bad pixels based on spline interpolation ###
            spline_func = IUS(pix_masked, flux_masked, check_finite=True)

            interp_flux = spline_func(pixels[:,j])
            if hotpix_mode=="interp":
                flux_tmp[mask_ord | mask_nan] = interp_flux[mask_ord | mask_nan]
            elif hotpix_mode=="nan":
                flux_tmp[mask_ord | mask_nan] = np.nan
        flux_tmp[mask_edge] = 0 #np.nan
        flux_new.append(flux_tmp)
    flux_new = np.array(flux_new).T
    return flux_new

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
