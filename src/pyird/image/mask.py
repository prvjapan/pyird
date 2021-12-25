import numpy as np
import tqdm

def trace(im, trace_func, y0, xmin, xmax, coeff):
    """make mask for trace parameters for multiorder
    
    Args:
       im: image
       trace_func: trace function
       x: x-array
       y0: y-offset
       xmin: xmin
       xmax: xmax
       coeff: coefficients

    Returns:
       mask

    Examples:        
        >>> from pyird.image.trace_function import trace_legendre
        >>> mask=trace(im, trace_legendre, y0, xmin, xmax, coeff)


    """
    x=[]
    for i in range(len(y0)):
        x.append(list(range(xmin[i],xmax[i]+1)))
    tl=trace_func(x, y0, xmin, xmax, coeff)
    mask=np.zeros_like(im,dtype=bool)    
    width=2
    nx,ny=np.shape(im)
    for i in tqdm.tqdm(range(len(y0))):
        tl_tmp=np.array(tl[i],dtype=int)
        for j,ix in enumerate(x[i]):
            iys=np.max([0,tl_tmp[j]-width])
            iye=np.min([ny,tl_tmp[j]+width+2])
            mask[ix,iys:iye]=True
    return mask[::-1,::-1]

