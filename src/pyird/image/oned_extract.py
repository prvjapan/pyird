import tqdm
import numpy as np
import pandas as pd


def flatten(im, trace_func, y0, xmin, xmax, coeff, inst='IRD',onepix=False,npix=2048,width=None):
    """make mask for trace parameters for multiorder.

    Args:
       im: image
       trace_func: trace function
       x: x-array
       y0: y-offset
       xmin: xmin
       xmax: xmax
       coeff: coefficients
       inst: instrument (IRD or REACH)
       onepix: extract the spectrum pixel by pixel in an aperture
       npix: number of pixels
       width: list of aperture widths ([width_start,width_end])

    Returns:
       raw multiorder spectra and multiorder pixel coordinate
       if onepix is True, return pandas DataFrame of spectra in each pixel
    """

    if width is None:
        if inst=='IRD':
            width_str = 2
            width_end = 4
        elif inst=='REACH':
            width_str = 2
            width_end = 3
    else:
        width_str = width[0]
        width_end = width[1]

    if len(y0)==21: #h band
        rotim = np.copy(im[::-1, ::-1])
    elif len(y0)==51: # yj band
        rotim = np.copy(im)

    x = []
    for i in range(len(y0)):
        x.append(list(range(xmin[i], xmax[i]+1)))
    tl = trace_func(x, y0, xmin, xmax, coeff)
    nx, ny = np.shape(im)

    spec, pixcoord, iys_all, iye_all = [], [], [], []
    orders=range(1,len(y0)+1); pixels=range(1,npix+1)
    df_zero = pd.DataFrame([],columns=pixels,index=orders)
    if onepix:
        df_onepix = {}
        for i in range(-width_str,width_end):
            df_onepix['ec%d'%(i)] = df_zero.copy()  ## initialize (need .copy())
    for i in tqdm.tqdm(range(len(y0))):
        tl_int = tl[i].astype(int)
        tl_decimal = tl[i] - tl_int
        eachspec = []
        eachpixcoord = []
        iys_tmp,iye_tmp = [], []
        for j, ix in enumerate(x[i]):
            if onepix:
                for k in range(-width_str,width_end):
                    if k<0:
                        iys = np.max([0, tl_int[j]+k])
                        iye = np.max([0, tl_int[j]+k+1])
                    else:
                        iys = np.min([ny-1, tl_int[j]+k])
                        iye = np.min([ny-1, tl_int[j]+k+1])
                    # At the ends of the aperture partial pixels are used. (cf. IRAF apall)
                    apsum = rotim[ix,iys]*(1-tl_decimal[j]) + rotim[ix,iye]*tl_decimal[j]
                    df_onepix['ec%d'%(k)].loc[i+1,ix+1] = apsum
            else:
                iys = np.max([0, tl_int[j]-width_str])
                iye = np.min([ny-1, tl_int[j]+width_end])
                # At the ends of the aperture partial pixels are used. (cf. IRAF apall)
                apsum = np.sum(rotim[ix, iys+1:iye]) + rotim[ix,iys]*(1-tl_decimal[j]) + rotim[ix,iye]*tl_decimal[j]
                eachspec.append(apsum)
                eachpixcoord.append(ix)
                iys_tmp.append(iys)
                iye_tmp.append(iye)
        spec.append(eachspec)
        pixcoord.append(eachpixcoord)
        iys_all.append(iys_tmp)
        iye_all.append(iye_tmp)
    if onepix:
        return df_onepix
    else:
        return spec, pixcoord, rotim, tl, iys_all, iye_all

if __name__ == '__main__':
    import numpy as np
    import pkg_resources
    from pyird.utils import irdstream
    import pathlib
    import matplotlib.pyplot as plt
    from pyird.image.mask import trace
    from pyird.image.trace_function import trace_legendre
    from pyird.image.pattern_model import median_XY_profile
    from pyird.io.iraf_trace import read_trace_file
    import astropy.io.fits as pyf

    datadir = pathlib.Path('/home/kawahara/pyird/data/samples/REACH/')
    anadir = pathlib.Path('/home/kawahara/pyird/data/samples/REACH/')
    target = irdstream.Stream2D('targets', datadir, anadir)
#    target.fitsid=[47077]
    target.fitsid = [47103]
    # Load an image
    for datapath in target.rawpath:
        im = pyf.open(str(datapath))[0].data

    # image for calibration
    pathC = (pkg_resources.resource_filename('pyird', 'data/samples/aprefC'))
    path_c = (pkg_resources.resource_filename('pyird', 'data/samples/apref_c'))
    y0, interp_function, xmin, xmax, coeff = read_trace_file([pathC, path_c])

    mask = trace(im, trace_legendre, y0, xmin, xmax, coeff)
    calim = np.copy(im)
    calim[mask] = np.nan

    model_im = median_XY_profile(calim, show=False)
    corrected_im = im-model_im

    # trace
    pathA = (pkg_resources.resource_filename('pyird', 'data/samples/aprefA'))
    y0, interp_function, xmin, xmax, coeff = read_trace_file(pathA)

    # flatten
    spec, pixcoord = flatten(
        corrected_im, trace_legendre, y0, xmin, xmax, coeff)

    fig = plt.figure(figsize=(12, 7))
    for i, esp in enumerate(spec):
        plt.plot(pixcoord[i], esp, label='order '+str(i))
    plt.xlabel('pixel coordinate')
    plt.ylabel('raw spectra')
    plt.legend()
    plt.savefig('fig_flatten.png')
    plt.show()
