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

    pathC=(pkg_resources.resource_filename('pyird', "data/samples/aprefC"))
    path_c=(pkg_resources.resource_filename('pyird', "data/samples/apref_c"))

    y0, interp_function, xmin, xmax, coeff=read_trace_file([pathC,path_c])

    #
    x=[]
    for i in range(len(y0)):
        x.append(list(range(xmin[i],xmax[i]+1)))
    #
    tl=trace_legendre(x, y0, xmin, xmax, coeff)

    datadir=pathlib.Path("/home/kawahara/pyird/data/samples/REACH/")
    anadir=pathlib.Path("/home/kawahara/pyird/data/samples/REACH/")
    reach=irdstream.Stream2D("targets",datadir,anadir)
    #reach.fitsid=[47103]
    reach.fitsid=[47077]

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

    import tqdm
    mask=np.zeros_like(im,dtype=bool)    
    width=2
    nx,ny=np.shape(im)
    for i in tqdm.tqdm(range(len(y0))):
        tl_tmp=np.array(tl[i],dtype=int)
        for j,ix in enumerate(x[i]):
            iys=np.max([0,tl_tmp[j]-width])
            iye=np.min([ny,tl_tmp[j]+width+2])
            mask[ix,iys:iye]=True 
            
    rotim=im[::-1,::-1]
    rotim[mask]=None

    if True:
        c=plt.imshow(rotim,vmin=-3.0,vmax=50.0,cmap="bone_r")
        plt.colorbar(c)
        for i in range(0,len(y0)):
            plt.plot(tl[i],x[i],alpha=1.0,color="red",ls="dotted",lw=2)
        plt.savefig("ontrace.png")
        plt.show()
