############################
## Code Name:  mask_maskimg
## Author(s): M. Kuzuhara (ABC)
## Email: m.kuzuhara@nao.ac.jp
## Last Update: 2019-01-08
## included to pyird by H. Kawahara
"""

This script makes a mask image, where spectra regions have counts 1000, using iraf "apmask" and a previously-generated iraf configration database for spectrum extraction such as the one made via "apall". 

To use this script, you need to put "database" directory on your working directory, and the database directory needs to have the correcting configuration database. 

usage>

python make_maskimg.py "NAME_of_configration" (like white_mmf1_18Oct_yj)

 
"""

import numpy as np
import astropy.io.fits as pyf
import pyraf.iraf as iraf
import sys
import time
import glob
import os

def makemsk(anadir,apfile,reffits):
    #mv to anadir
    currentdir=os.getcwd() 
    os.chdir(anadir)
    outputf="mask_from_" + reffits
    
    ## check which your target data, YJ or H ##
    ## and determine where the blocked region id ##
    
    line_num = len(open(apfile, "r").readlines())

    #if line_num > 1000:

    block_region_y0 = 2000
    block_region_y1 = 2048

    #else:

    #   block_region_y0 = 0
    #   block_region_y1 = 40
      
    ########################################### 
    ## make dummy image ##
    dummy = np.ones((2048,2048))*1000.
    dummy[block_region_y0:block_region_y1,:] = 0.
    dummy_fits_name = "dummy_formsk_" + str(time.time()) + ".fits"
    pyf.writeto(dummy_fits_name, dummy) 
    ######################

    iraf.imred(Stdout=1)
    iraf.echelle(Stdout=1)
    temporary_mask_pl = "tmp_" + str(time.time())
    iraf.apmask(input = dummy_fits_name, output = temporary_mask_pl, 
                apertures=apfile,references = reffits,
                interactive='no', find='no', recenter='no', resize='no',
                edit='no', trace='no', fittrace='no', 
                mask='yes', line='INDEF', nsum=10, buffer=0.0)
    iraf.imar(dummy_fits_name, "*", temporary_mask_pl + ".pl", \
              outputf + ".fits")
    iraf.imdel(dummy_fits_name + "," + temporary_mask_pl + ".pl")

    ###### handle for H-band data #########
    if line_num < 1000:
       print("here")
       msk_tmp = pyf.open(outputf + ".fits")[0].data
       msk_tmp = msk_tmp[::-1,::-1]
       iraf.imdel(outputf + ".fits")
       pyf.writeto(outputf + ".fits", msk_tmp)
       
    os.chdir(currentdir)
    return outputf

#if __name__ == "__main__":
#
#   makemsk(sys.argv[1])
