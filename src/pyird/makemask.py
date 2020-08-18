############################
## Code Name:  mask_maskimg
## Author(s): M. Kuzuhara (ABC)
## Email: m.kuzuhara@nao.ac.jp
## Last Update: 2019-01-08
## included and modified for pyird by H. Kawahara
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
import pathlib
def makemsk(anadir,reffitlist=[]):
    #mv to anadir
    currentdir=os.getcwd() 
    os.chdir(anadir)
    outputf="mask_from"
    for reffits in reffitlist:
       outputf=outputf+"_"+reffits
       
    if pathlib.Path(outputf).exists():
        iraf.imdel(outputf)
    block_region_y0 = 2000
    block_region_y1 = 2048
    ########################################### 
    iraf.imred(Stdout=1)
    iraf.echelle(Stdout=1)
    dumout=[]
    for l,reff in enumerate(reffitlist):
        dummy = np.ones((2048,2048))*1000.
        dummy[block_region_y0:block_region_y1,:] = 0.
        dummy_fits_name = "dummy_formsk_" + str(time.time()) + ".fits"
        pyf.writeto(dummy_fits_name, dummy) 

        temporary_mask_pl = "tmp_" + str(time.time())
        iraf.apmask(input = dummy_fits_name, output = temporary_mask_pl, 
                    references = reff,
                    interactive='no', find='no', recenter='no', resize='no',
                    edit='no', trace='no', fittrace='no', 
                    mask='yes', line='INDEF', nsum=10, buffer=0.0)
        dummy_fits_name_out = "dummy_formsk_out" + str(time.time()) + ".fits"
        iraf.imar(dummy_fits_name, "*", temporary_mask_pl + ".pl", \
                  dummy_fits_name_out)
        dumout.append(dummy_fits_name_out)

    if len(reffitlist)>1:
        dummy_fits_name=dumout[0]
        for l in range(len(reffitlist)-1):
            if l==len(reffitlist)-2:
                dummy_fits_name_swap=outputf + ".fits"
            else:
                dummy_fits_name_swap="dummy_foradd_" + str(time.time()) + ".fits"
            iraf.imar(dummy_fits_name, "+", dumout[l+1], \
                      dummy_fits_name_swap)
            dummy_fits_name=dummy_fits_name_swap
        

    ###### handle for H-band data #########
    msk_tmp = pyf.open(outputf + ".fits")[0].data
    msk_tmp=np.array(msk_tmp)
    msk_tmp[msk_tmp>1000]=1000.0
    iraf.imdel(outputf + ".fits")
    pyf.writeto(outputf + ".fits", msk_tmp)
       
    os.chdir(currentdir)
    return outputf

#if __name__ == "__main__":
#
#   makemsk(sys.argv[1])
