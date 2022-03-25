import numpy as np
import tqdm


def trace_from_iraf_trace_file(im, pathlist):
    """make mask from the iraf trace file.

    Args:
       im: image
       pathlist: trace files path list

    Returns:
       mask image (same shape as im)

    Note:
       Currently, Only legendre is supported.

    """
    from pyird.io.iraf_trace import read_trace_file
    from pyird.image.trace_function import trace_legendre
    
    y0, interp_function, xmin, xmax, coeff = read_trace_file(pathlist)
    unique_interp=np.unique(interp_function)
    if len(unique_interp)==1 and unique_interp[0]==2:
        mask = trace(im, trace_legendre, y0, xmin, xmax, coeff)
    else:
        print("Error: other interpolation function than legendre is not supported yet.")
        
    return mask
    
def trace(im, trace_func, y0, xmin, xmax, coeff):
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
       mask image (same shape as im)

    Examples:        
        >>> from pyird.image.trace_function import trace_legendre
        >>> mask=trace(im, trace_legendre, y0, xmin, xmax, coeff)
    """

    x = []
    for i in range(len(y0)):
        x.append(list(range(xmin[i], xmax[i]+1)))
    tl = trace_func(x, y0, xmin, xmax, coeff)
    mask = np.zeros_like(im, dtype=bool)
    width = 2
    nx, ny = np.shape(im)
    for i in tqdm.tqdm(range(len(y0))):
        tl_tmp = np.array(tl[i], dtype=int)
        for j, ix in enumerate(x[i]):
            iys = np.max([0, tl_tmp[j]-width])
            iye = np.min([ny, tl_tmp[j]+width+2])
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
