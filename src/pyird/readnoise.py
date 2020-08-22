import numpy as np
from pyird import gp2d

def RNestimate_OEGP(img,xscale=10,yscale=5,sigma=0.01):
    #odd/even GP method by H.Kawahara
    from numpy import linalg as LA
    import scipy
    import numpy as np
    import time
    #####################
    #GP model for o/e pix and o/e chan
    nchan=32 # number od channel
    noechan=int(nchan/2) #number of odd or even channel
    nRN=64 # number of yaxis
    npix=2048 # number of xaxis (pixel)
    nSRN=int(nRN/2)
    ##fold
    subcube=[]
    for jchan in range(0,nchan):
        arr=img[jchan*nRN:(jchan+1)*nRN,:]
        if np.mod(jchan,2)==0:
            subcube.append(arr[0::2,:])
            subcube.append(arr[1::2,:])
        else:
            subcube.append(arr[-1::-2,:])
            subcube.append(arr[-2::-2,:])
        
    #GP
    subcube_median=np.nanmedian(subcube,axis=(1,2))
    subcubex=subcube-subcube_median[:,np.newaxis,np.newaxis]
    subarray=np.nanmedian(subcubex,axis=0)
    subarray=subarray[:,4:-4]#remove edges

    rectmp=gp2d.GP2D(subarray,sigma,xscale,yscale)
    rec=np.zeros((nSRN,npix))
    rec[:,4:-4]=rectmp

    ### RECONSTRUCT
    recimg=np.zeros(np.shape(img))
    for jchan in range(0,nchan):
        arr=np.zeros((nRN,npix))
        arr[0::2,:]=rec[:,:]+subcube_median[2*jchan]#+3
        arr[1::2,:]=rec[:,:]+subcube_median[2*jchan+1]
        if np.mod(jchan,2)==0:
            recimg[jchan*nRN:(jchan+1)*nRN,:]=arr
        else:
            recimg[jchan*nRN:(jchan+1)*nRN,:]=arr[::-1,:]
            
    return recimg
#    cimg=imgorig-recimg
