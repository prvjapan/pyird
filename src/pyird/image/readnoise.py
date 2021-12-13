import numpy as np
from pyird import gp2d

def coarse_gp(array,subarray,xscale,yscale,Ncor):
    # coaese graining and 2D GP recovery
    ## remove outlier
    each=(array-np.median(array))
    mask=np.abs(each)>5.0*np.std(each)
    marray=np.copy(array)
    marray[mask]=None    
    stacked_array=[]
    for j in range(0,Ncor):
        stacked_array.append(marray[:,j::Ncor])
    stacked_array=np.array(stacked_array)
    coarsed_array=np.nanmedian(stacked_array,axis=0)
    coarsed_array[coarsed_array!=coarsed_array]=np.nanmedian(coarsed_array)    
    sigma=0.001
    GLP=gp2d.GP2Dcross(coarsed_array,subarray,sigma,xscale,yscale)
    return GLP


def RNestimate_OEGP(img,xscale=10,yscale=5,sigma=0.01,cgdmode="gp"):
    #OEGP (Odd/Even+Gaussian Process) method by H.Kawahara
    
    #cgdmode="gp": channel global distribtuion infered by GP2d
    #cgdmode="median": channel global distribtuion infered by by median
    
    from numpy import linalg as LA
    import scipy
    import numpy as np
    import time

    #####################
    # infer common subprofile
    nchan=32 # number od channel
    noechan=int(nchan/2) #number of odd or even channel
    nRN=64 # number of yaxis
    npix=2048 # number of xaxis (pixel)
    nSRN=int(nRN/2)

    ##folding
    subcube=[]
    for jchan in range(0,nchan):
        arr=img[jchan*nRN:(jchan+1)*nRN,:]
        if np.mod(jchan,2)==0:
            subcube.append(arr[0::2,:])
            subcube.append(arr[1::2,:])
        else:
            subcube.append(arr[-1::-2,:])
            subcube.append(arr[-2::-2,:])
        
    subcube_median=np.nanmedian(subcube,axis=(1,2))
    subcubex=subcube-subcube_median[:,np.newaxis,np.newaxis]
    subarray=np.nanmedian(subcubex,axis=0)
    subarray=subarray[:,4:-4]#remove edges
    rec=np.zeros((nSRN,npix)) #recovered common subprofile
    rec[:,4:-4]=subarray
    #######################

    #######################
    if cgdmode=="gp":        
        # Channel Global Distribution
        CGDa=[]
        Ncor=64
        xs=256
        ys=64
        for jchan in range(0,nchan):
            arr=img[jchan*nRN:(jchan+1)*nRN,:]
            if np.mod(jchan,2)==0:
                cgda=arr[0::2,:]-rec    
                CGD=coarse_gp(cgda,subarray,xs,ys,Ncor)
                CGDa.append(CGD)
                
                cgda=arr[1::2,:]-rec    
                CGD=coarse_gp(cgda,subarray,xs,ys,Ncor)
                CGDa.append(CGD)
            else:
                cgda=arr[-1::-2,:]-rec
                CGD=coarse_gp(cgda,subarray,xs,ys,Ncor)
                CGDa.append(CGD)
                
                cgda=arr[-2::-2,:]-rec
                CGD=coarse_gp(cgda,subarray,xs,ys,Ncor)
                CGDa.append(CGD)
    #######################
    
    ######################
    #Recovering an image
    recimg=np.zeros(np.shape(img))
    iap=0
    for jchan in range(0,nchan):
        arr=np.zeros((nRN,npix))
        if cgdmode=="gp":
            CGD=np.zeros(np.shape(rec))
            CGD[:,4:-4]=CGDa[iap]
            arr[0::2,:]=rec[:,:]+CGD
        elif cgdmode=="median":
            arr[0::2,:]=rec[:,:]+subcube_median[2*jchan]
        else:
            sys.exit("No availbale cgdmode.")
        iap=iap+1

        if cgdmode=="gp":
            CGD=np.zeros(np.shape(rec))
            CGD[:,4:-4]=CGDa[iap]
            arr[1::2,:]=rec[:,:]+CGD
        elif cgdmode=="median":
            arr[1::2,:]=rec[:,:]+subcube_median[2*jchan+1]
        
        iap=iap+1
        if np.mod(jchan,2)==0:
            recimg[jchan*nRN:(jchan+1)*nRN,:]=arr
        else:
            recimg[jchan*nRN:(jchan+1)*nRN,:]=arr[::-1,:] #recovered image
        
    return recimg
