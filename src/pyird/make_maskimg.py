#!/usr/bin/env python
############################
## Code Name:  mask_maskimg
## Author(s): M. Kuzuhara (ABC)
## Email: m.kuzuhara@nao.ac.jp
## Last Update: 2019-01-08
"""

This script makes a mask image, where spectra regions have counts 1000, using iraf "apmask" and a previously-generated iraf configration database for spectrum extraction such as the one made via "apall". 

To use this script, you need to put "database" directory on your working directory, and the database directory needs to have the correcting configuration database. 

usage>

python make_maskimg.py "NAME_of_configration" (like white_mmf1_18Oct_yj)

 
"""

import numpy as np
import pyfits as pyf
import pyraf.iraf as iraf
import sys
import time
import glob

## rawimg is fits image, iraf_config is IRAF database file made through apperture tracing (e.g., apall)
def makemsk(iraf_config):

    ## check which your target data, YJ or H ##
    ## and determine where the blocked region id ##
    
    config_file = glob.glob("./database/*"+iraf_config)[0]
    line_num = len(open(config_file, "r").readlines())

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
                references = iraf_config,
                interactive='no', find='no', recenter='no', resize='no',
                edit='no', trace='no', fittrace='no', 
                mask='yes', line='INDEF', nsum=10, buffer=0.0)
    iraf.imar(dummy_fits_name, "*", temporary_mask_pl + ".pl", \
              "mask_from_" + iraf_config + ".fits")
    iraf.imdel(dummy_fits_name + "," + temporary_mask_pl + ".pl")

    ###### handle for H-band data #########

    if line_num < 1000:

       print "here"
       msk_tmp = pyf.open("mask_from_" + iraf_config + ".fits")[0].data
       msk_tmp = msk_tmp[::-1,::-1]
       iraf.imdel("mask_from_" + iraf_config + ".fits")
       pyf.writeto("mask_from_" + iraf_config + ".fits", msk_tmp)


if __name__ == "__main__":

   makemsk(sys.argv[1])
