import tqdm
import numpy as np


def flatten(im, trace_func, y0, xmin, xmax, coeff, inst='IRD'):
    """make mask for trace parameters for multiorder.

    Args:
       im: image
       trace_func: trace function
       x: x-array
       y0: y-offset
       xmin: xmin
       xmax: xmax
       coeff: coefficients

    Returns:
       raw multiorder spectra
       multiorder pixel coordinate
    """

    if len(y0)==21: #h band
        rotim = np.copy(im[::-1, ::-1])
        width_str = 3
        if inst=='REACH':
            width_end = 3
        elif inst=='IRD':
            width_end = 4
    elif len(y0)==51: # yj band
        rotim = np.copy(im)
        if inst=='REACH':
            width_str = 4
            width_end = 5
        elif inst=='IRD':
            width_str = 5
            width_end = 6

    x = []
    for i in range(len(y0)):
        x.append(list(range(xmin[i], xmax[i]+1)))
    tl = trace_func(x, y0, xmin, xmax, coeff)
    nx, ny = np.shape(im)
    spec = []
    pixcoord = []
    iys_all,iye_all = [], []
    for i in tqdm.tqdm(range(len(y0))):
        tl_tmp = np.array(tl[i], dtype=int)
        eachspec = []
        eachpixcoord = []
        iys_tmp,iye_tmp = [], []
        for j, ix in enumerate(x[i]):
            iys = np.max([0, tl_tmp[j]-width_str])
            iye = np.min([ny, tl_tmp[j]+width_end])
            eachspec.append(np.sum(rotim[ix, iys:iye]))
            eachpixcoord.append(ix)
            iys_tmp.append(iys)
            iye_tmp.append(iye)
        spec.append(eachspec)
        pixcoord.append(eachpixcoord)
        iys_all.append(iys_tmp)
        iye_all.append(iye_tmp)
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
