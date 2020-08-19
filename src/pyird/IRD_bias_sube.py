#!/usr/bin/env python
############################

"""
This module can be used to correct stripe-pattern noises in IRD data. To subtract the biases (vertical stripes), the reference pixel or light-insensitive (botom part of image) regions are used. Then, median count of each channel are subtracted.  Meanwhile, the horizontal stripe patters are subtracted using the unilluminated pixels. 

Hot pixel mask can be used.  This code identify pixels with counts of 1.

"""
# Author: Masayuki Kuzuhara (Astrobiology Center/NAOJ; m.kuzuhara@nao.ac.jp)
# Last update: 2017-January-5

############################

import sys
import numpy as np
#import pyfits as pyf
import astropy.io.fits as pyf
import scipy.stats as stats
from datetime import datetime
import pathlib

########################
#### header edit #######

def header_edit(hd, method):

    current_date = str(datetime.now())[0:-5]
    hd.set("rbdate", current_date, \
                     "Fits modification date for removing bias") 
    hd.set("rbmeth", method, \
                     "Method for removing bias") 

    return hd

### subtract the biases in each readout channel; #####
### reference pixels (+/- 4 edge pixels) are used #### 
### This is recommended for IRD currently         ####

### im: image
### hd: header

def bias_subtract_reference(im, hd, Nch = 32):
 
    xsize = np.shape(im)[1]; ysize = np.shape(im)[0]
    dummy = np.ones((xsize,ysize))    
    ch_pix_num = xsize/Nch

    ###### bias removal for each stripe ######
    for ch_num in range(Nch):

        ysize = int(ysize)
        x_start = int(ch_pix_num*ch_num)
        x_end = int(ch_pix_num*(ch_num+1))    

        stripe = im[0:ysize, x_start:x_end]
        ref1 = im[0:4, x_start:x_end]
        ref2 = im[2044:ysize, x_start:x_end]

        ref_all = np.concatenate( (ref1, ref2) )
        
        clipped = stats.sigmaclip(ref_all, 3.0, 3.0)[0]

#        mean = np.nanmean(clipped)
        mean = np.nanmedian(clipped) #by H.K.

        unbias_stripe = stripe - mean

        dummy[0:ysize,x_start: x_end] = unbias_stripe

    hd = header_edit(hd, 'reference')

    return dummy, hd


#######################################################################
### subtract bias of each channel;                                  ###
### the bias is estimated by using the light-insensitive regions    ###
#######################################################################

### Yfaint_pix is the y range where scatterd light less emarged ####

def bias_subtract_block(im, Nch = 32, Yfaint_pix = 40, hotpixmask = None):
 
    xsize = np.shape(im)[1]; ysize = np.shape(im)[0]

    im[0:ysize, 0:4] = np.nan
    im[0:ysize, xsize-4:xsize] = np.nan

    ############ hot pixel masking ########################

    if hotpixmask is not None:

       mask_index = np.where(hotpixmask == 1.0)
       im[mask_index] = np.nan

    ch_pix_num = xsize/Nch

    ###### bias subtraction begins #######
    
    for i in range(Nch):

        stripe = im[0:ysize,ch_pix_num*i:ch_pix_num*(i+1)] 
        lessiluminated = im[4:Yfaint_pix, ch_pix_num*i:ch_pix_num*(i+1)]
        notnanindex = np.where(lessiluminated == lessiluminated)

        lessiluminated = lessiluminated[notnanindex]

        clip = stats.sigmaclip(lessiluminated, 3.0, 3.0)[0] 

        im[0:ysize, ch_pix_num*i:ch_pix_num*(i+1)] = \
                stripe - np.mean(clip)

    ###################################

    return im
 
#######################################################################
### subtract bias of each channel;                                  ###
### the bias is estimated using whole pixels except reference       ###
### pixels in each channel                                          ### 
#######################################################################

def bias_subtract_whole(im, Nch = 32, hotpixmask = None):

    xsize = np.shape(im)[1]; ysize = np.shape(im)[0]

    im[0:ysize, 0:4] = np.nan
    im[0:ysize, xsize-4:xsize] = np.nan

    ############ hot pixel masking ########################

    if hotpixmask is not None:

       mask_index = np.where(hotpixmask == 1.0)
       im[mask_index] = np.nan

    ############ bias subtraction begins  ###############  
    ch_pix_num = xsize/Nch

    for i in range(Nch):

        stripe = im[0:ysize,ch_pix_num*i:ch_pix_num*(i+1)] 
        notnanindex = np.where(stripe == stripe)

        tmp = stripe[notnanindex]

        clip = stats.sigmaclip(tmp, 5.0, 5.0)[0] 

        im[0:ysize, ch_pix_num*i:ch_pix_num*(i+1)] = \
                stripe - np.mean(clip)

    


    return im
   
#####################################################################

def main(outdir, data_list, method, rot=None, hotpix_im = None, stripe_direction = 'horizon'):
    import tqdm
    ######### open hot-pixel image ############


    if hotpix_im is not None:
       hotpix_im = pyf.open(hotpix_im)[0].data
       
       ###### transpose image if the stripe directioin is horizontal ######    

       if stripe_direction == 'horizon':
                hotpix_im = np.transpose(hotpix_im)
    
    if method == 'block':
 
       print("----------------------------------------------------------")
       print(" ")
       print("Bias subtraction using science pixels near edge" )
       print(" ")
       print("----------------------------------------------------------")
    data_list_noexist=[]
    skip=0
    for data in data_list:
        outdata_name = data.name.replace(".fits","") + '_rb.fits'
        if not (outdir/outdata_name).exists():
            data_list_noexist.append(data)
        else:
            skip=skip+1
            
    if skip>1:
        print("Bias Correction: Skipped "+str(skip)+" files because they already exist.")
    
    for data in tqdm.tqdm(data_list_noexist):
        #        print "processing for " + str(data)
        im = pyf.open(str(data))[0].data
        hd = pyf.open(str(data))[0].header

        ###### transpose image if the stripe directioin is horizontal ######    

        if stripe_direction == 'horizon':

           im = im.transpose()  

        ####################################################################
        
        if method == 'block': 
    
           im = bias_subtract_block(im, hotpixmask = hotpix_im)

        elif method == 'whole':

           im = bias_subtract_whole(im, hotpixmask = hotpix_im)

        elif method == 'reference':

           im, hd = bias_subtract_reference(im,hd)
 
        #### output ##### 
        outdata_name = data.name.replace(".fits","") + '_rb.fits'
     
        ###### transpose image if the stripe directioin is horizontal ######    

        if stripe_direction == 'horizon':

           im = im.transpose()  

        ####################################################################
        #180 deg rotation, usually for H band
        if rot == "r":
            im=im[::-1,::-1]
        
        try:
            pyf.writeto(outdir/outdata_name, im, hd)
        except:
            print("Output file = "+str(outdir/outdata_name)+" already Exists.")

if __name__ == "__main__":
 
    #method = 'whole'
    #method = 'block'
    method = 'reference'

    input_data = sys.argv[1]
    

    if method == 'reference':

       hotpix_img = None

    main(input_data, method, hotpix_im = hotpix_img)

