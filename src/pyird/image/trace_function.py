import numpy as np

def trace_legendre(x, y0, xmin, xmax, coeff):
    """ trace Legendre function 

    Args:
       x: x-array
       y0: y-offset
       xmin: xmin
       xmax: xmax
       coeff: Legendre polynomial coefficients

    """
    from numpy.polynomial.legendre import legval
    norder=len(y0)
    trace_lines=[]
    for i in range(0,norder):
        x_=np.array(x[i])
        xmax_=xmax[i]
        xmin_=xmin[i]        
        nx=(2.*x_ - (xmax_+xmin_))/(xmax_-xmin_)
        f=legval(nx, coeff[i])+y0[i]-1
        trace_lines.append(f)

    return trace_lines
    
if __name__=="__main__":
    import pkg_resources
    from pyird.utils import irdstream
    from pyird.io.iraf_trace import read_trace_file
    from pyird.image.channel import image_to_channel_cube, channel_cube_to_image, eopixel_split, eopixel_combine
    from pyird.image.bias import bias_subtract

    import numpy as np
    import matplotlib.pyplot as plt
    import pathlib

    path=(pkg_resources.resource_filename('pyird', "data/samples/aprefA"))
    y0, interp_function, xmin, xmax, coeff=read_trace_file(path)
    x=[]
    for i in range(len(y0)):
        x.append(list(range(xmin[i],xmax[i]+1)))
    #x=np.linspace(xmin,xmax,1000)
    
    tl=trace_legendre(x, y0, xmin, xmax, coeff)

    datadir=pathlib.Path("/home/kawahara/pyird/data/samples/REACH/")
    anadir=pathlib.Path("/home/kawahara/pyird/data/samples/REACH/")
    reach=irdstream.Stream2D("targets",datadir,anadir)
    reach.fitsid=[47103]

    import astropy.io.fits as pyf
    for datapath in reach.rawpath:
        im = pyf.open(str(datapath))[0].data
    channel_cube=image_to_channel_cube(im,revert=True)
    bs_channel_cube, bias=bias_subtract(channel_cube)
    im=channel_cube_to_image(bs_channel_cube)

    if True:
        c=plt.imshow(im[::-1,::-1],vmin=-3.0,vmax=50.0,cmap="bone")
        plt.colorbar(c)
        for i in range(0,len(y0)):
            plt.plot(tl[i],x[i],alpha=1.0,color="red",ls="dotted",lw=2)
        plt.savefig("ontrace.png")
        plt.show()

