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



def extract_pattern_OEGP(img,channel_cube,xscale=10,yscale=5,sigma=0.01,cgdmode="gp"):
    #OEGP (Odd/Even+Gaussian Process) method by H.Kawahara
    
    #cgdmode="gp": channel global distribtuion infered by GP2d
    #cgdmode="median": channel global distribtuion infered by by median
    
    from numpy import linalg as LA
    import scipy
    import numpy as np
    import time

    Nch, ch_pix_num, xsize=np.shape(channel_cube)
    noechan=int(Nch/2) #number of odd or even channel
    nSRN=int(ch_pix_num/2)

    a=channel_cube[0::2,0::2,:]
    b=channel_cube[0::2,1::2,:]
    c=channel_cube[1::2,-1::-2,:]
    d=channel_cube[1::2,-2::-2,:]

    
    ##folding
    subcubea=[]
    subcubeb=[]
    subcubec=[]
    subcubed=[]

    for jchan in range(0,Nch):
        arr=img[jchan*ch_pix_num:(jchan+1)*ch_pix_num,:]
        if np.mod(jchan,2)==0:
            subcubea.append(arr[0::2,:])
            subcubeb.append(arr[1::2,:])
        else:
            subcubec.append(arr[-1::-2,:])
            subcubed.append(arr[-2::-2,:])
    subcubea=np.array(subcubea)
    print(subcubea[0,0,:])
    print(a[0,0,:])
    print(np.shape(a))
    import sys
    sys.exit()
            
    #here
    subcube_median=np.nanmedian(subcube,axis=(1,2))
    subcubex=subcube-subcube_median[:,np.newaxis,np.newaxis]
    subarray=np.nanmedian(subcubex,axis=0)
    subarray=subarray[:,4:-4]#remove margin
    rec=np.zeros((nSRN,xsize)) #recovered common subprofile
    rec[:,4:-4]=subarray
    #######################

    #######################
    if cgdmode=="gp":        
        # Channel Global Distribution
        CGDa=[]
        Ncor=64
        xs=256
        ys=64
        for jchan in range(0,Nch):
            arr=img[jchan*ch_pix_num:(jchan+1)*ch_pix_num,:]
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
    for jchan in range(0,Nch):
        arr=np.zeros((ch_pix_num,xsize))
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
            recimg[jchan*ch_pix_num:(jchan+1)*ch_pix_num,:]=arr
        else:
            recimg[jchan*ch_pix_num:(jchan+1)*ch_pix_num,:]=arr[::-1,:] #recovered image
        
    return recimg



if __name__=="__main__":
    import numpy as np
    from pyird.image.channel import image_to_channel_cube, channel_cube_to_image
    from pyird.image.bias import bias_subtract
    from pyird.utils import fitsset,irdstream
    import astropy.io.fits as pyf
    import tqdm
    import pathlib
    import matplotlib.pyplot as plt
    from scipy.stats import median_absolute_deviation as mad
    
    mode="YJ"
    datadir=pathlib.Path("/home/kawahara/pyird/data/dark/")
    anadir=pathlib.Path("/home/kawahara/pyird/data/dark/") 
    dark=irdstream.Stream2D("targets",datadir,anadir)
    dark.fitsid=[41018]    
    if mode=="H":
        dark.fitsid_increment()        

    print(dark.rawpath)
    for data in tqdm.tqdm(dark.rawpath):
        im = pyf.open(str(data))[0].data
        hd = pyf.open(str(data))[0].header

    
    #np.random.seed(1)
    #a=np.random.normal(0.0,1.0,(2048,2048))    
    channel_cube=image_to_channel_cube(im)
    c=bias_subtract(channel_cube)

    if False:
        image_rmbias=channel_cube_to_image(c)
        fig=plt.figure()
        ax1=fig.add_subplot(121)
        cc=ax1.imshow(im,vmin=-15.0,vmax=-8.0)
        plt.colorbar(cc)
        ax2=fig.add_subplot(122)
        cc=ax2.imshow(image_rmbias,vmin=-3.0,vmax=4.0)
        plt.colorbar(cc)    
        print(mad(im.flatten()),"->",mad(image_rmbias.flatten()))
        plt.show()

    extract_pattern_OEGP(a.T,channel_cube)
