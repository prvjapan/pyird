import numpy as np
import tqdm


def trace_from_iraf_trace_file(pathlist, mask_shape=None):
    """make mask from the iraf trace file.

    Args:
       pathlist: trace files path list
       mask_shape: (optional) shape of mask, c.f. np.shape(image)

    Returns:
       mask image (same shape as im)

    Note:
       Currently, Only legendre is supported.
    """
    from pyird.io.iraf_trace import read_trace_file
    from pyird.image.trace_function import trace_legendre

    y0, interp_function, xmin, xmax, coeff = read_trace_file(pathlist)
    unique_interp = np.unique(interp_function)
    if len(unique_interp) == 1 and unique_interp[0] == 2:
        mask = trace(trace_legendre, y0, xmin, xmax, coeff, mask_shape=mask_shape)
    else:
        print('Error: other interpolation function than legendre is not supported yet.')

    return mask


def trace(trace_func, y0, xmin, xmax, coeff, mask_shape=None, inst='IRD'):
    """make mask for trace parameters for multiorder.

    Args:
       trace_func: trace function
       x: x-array
       y0: y-offset
       xmin: xmin
       xmax: xmax
       coeff: coefficients
       mask_shape: (optional) shape of mask, c.f. np.shape(image)
       inst: IRD or REACH

    Returns:
       mask image (same shape as im)

    Examples:
        >>> from pyird.image.trace_function import trace_legendre
        >>> mask=trace(im, trace_legendre, y0, xmin, xmax, coeff)
    """
    if mask_shape is None:
        mask_shape=(2048,2048)

    x = []
    for i in range(len(y0)):
        x.append(list(range(xmin[i], xmax[i]+1)))
    tl = trace_func(x, y0, xmin, xmax, coeff)
#    mask = np.zeros_like(im, dtype=bool)
    mask = np.zeros(mask_shape, dtype=bool)
    if len(y0)==21 or len(y0)==42: #h band
        width_str = 3
        if inst=='REACH':
            width_end = 3
        elif inst=='IRD':
            width_end = 4
    else: #if len(y0)==51 or len(y0)==102: # yj band
        if inst=='REACH':
            width_str = 4
            width_end = 5
        elif inst=='IRD':
            width_str = 5
            width_end = 6
    nx, ny = mask_shape
    for i in tqdm.tqdm(range(len(y0))):
        tl_tmp = np.array(tl[i], dtype=int)
        for j, ix in enumerate(x[i]):
            iys = np.max([0, tl_tmp[j]-width_str])
            iye = np.min([ny, tl_tmp[j]+width_end])
            mask[ix, iys:iye] = True
    return mask[::-1, ::-1]


if __name__ == '__main__':
    import numpy as np
    import pkg_resources
    from pyird.utils import irdstream
    import pathlib
    import matplotlib.pyplot as plt
    from pyird.image.trace_function import trace_legendre
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
    pathA = (pkg_resources.resource_filename('pyird', 'data/samples/aprefA'))
    y0, interp_function, xmin, xmax, coeff = read_trace_file(pathA)
    mask = trace(im, trace_legendre, y0, xmin, xmax, coeff)

    immasked = np.copy(im)
    immasked[mask] = None

    fig = plt.figure()
    ax = fig.add_subplot(121)
    ax.imshow(im, vmin=-3.0, vmax=40.0)
    ax.set_title('raw')
    ax2 = fig.add_subplot(122)
    ax2.imshow(immasked, vmin=-3.0, vmax=40.0)
    ax2.set_title('masked')
    plt.show()
