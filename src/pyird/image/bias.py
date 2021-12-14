""" Bias subtraction
   
    - originally developed by M. Kuzuhara


"""
import numpy as np
import scipy.stats as stats

def bias_subtract_reference(image):
    """
    Args:
       image: 2D image

    """
    Nch = 32 # number of the channels of IRD detector
    ysize, xsize = np.shape(image)
    image_rmbias = np.ones((xsize,ysize))    
    ch_pix_num = xsize/Nch

    ###### bias removal for each stripe ######
    for ch_num in range(Nch):
        x_start = int(ch_pix_num*ch_num)
        x_end = int(ch_pix_num*(ch_num+1))    
        stripe = image[0:ysize, x_start:x_end]
        
        
        ref1 = image[0:4, x_start:x_end]
        ref2 = image[2044:ysize, x_start:x_end]
        ref_all = np.concatenate( (ref1, ref2) )        
        clipped = stats.sigmaclip(ref_all, 3.0, 3.0)[0]
        mean = np.nanmedian(clipped) #by H.K.
        unbias_stripe = stripe - mean
        image_rmbias[0:ysize,x_start: x_end] = unbias_stripe

    return image_rmbias

def bias_subtract(image):
    """
    Args:
       image: 2D image

    """
    Nch = 32 # number of the channels of IRD detector
    margin = 4 # margin width of the bias estimate zone
    snclip = 3.0 #clipping value
    ysize, xsize = np.shape(image)
    ch_pix_num = int(xsize/Nch)
    image_split=np.array(np.split(image.T,Nch)) # (Nch, ch_pix_num, 2048)
    #image=(image_split.reshape(Nch*ch_pix_num, 2048)).T # recover
    ref1x=(image_split[:,:,0:margin])
    ref2x=(image_split[:,:,ysize-margin:ysize])
    ref=np.concatenate([ref1x,ref2x],axis=2) # both margin image (Nch, ch_pix_num, 2*margin)
    meancp=[]
    for ch_num in range(Nch):
        cp=stats.sigmaclip(ref[ch_num,:,:], snclip, snclip)[0]
        meancp.append(np.nanmedian(cp))
    meancp=np.array(meancp)
    unbias_stripe=image_split-meancp[:,np.newaxis,np.newaxis]
    image_rmbias=(unbias_stripe.reshape(Nch*ch_pix_num, 2048)).T
    return image_rmbias
    

if __name__=="__main__":
    import numpy as np
    import time
    a=np.random.normal(0.0,1.0,(2048,2048))
    #a=np.array(range(0,2048*2048)).reshape(2048,2048)
    ts=time.time()
    b=bias_subtract_reference(a)
    te=time.time()
    print(ts-te)
    
    ts=time.time()
    c=bias_subtract(a)
    te=time.time()
    print(ts-te)
