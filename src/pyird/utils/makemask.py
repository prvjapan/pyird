############################
"""

Note:
    This script makes a mask image, where spectra regions have counts 1000, using iraf "apmask" and a previously-generated iraf configration database for spectrum extraction such as the one made via "apall". 

- Code Name:  mask_maskimg
- Author(s): M. Kuzuhara (ABC)
- Email: m.kuzuhara@nao.ac.jp
- Last Update: 2019-01-08
- included and modified for pyird by H. Kawahara

 
"""

import numpy as np
import pyraf.iraf as iraf
import sys
import time
import pathlib
from astropy.io import fits
import sys


def makemsk(anadir, reffitlist=[], directfitslist=[], directfits_criterion=[], directrot=[]):
    """make mask.

    Args:
        anadir: directory for analysis
        reffitlist: mask fits made by apall
        directfitslist: direct fits
        directfits_criterion: criterion for directfits
        directrot: if True, rotate image in 90 deg.

    Examples:
        makemask.makemsk(str(anadir),["refC"],["hotpixel_200709_h","MMF"],[0.99,20.0],[True,False])
    """
    # Moving analysis directory
    iraf.cd(str(anadir))
    print('Working at ', iraf.pwd())

    ### uniform mask ###
    if len(reffitlist) == 0 and len(directfitslist) == 0:
        dummy = np.zeros((2048, 2048))
        fits.writeto('mask_uniform.fits', dummy)
        return

    outputf = 'mask_from'
    for reffits in reffitlist:
        outputf = outputf+'_'+reffits

    if len(directfitslist) > 0:
        for reffits in directfitslist:
            outputf = outputf+'_D'+reffits

    if pathlib.Path(outputf).exists():
        iraf.imdel(outputf)
    block_region_y0 = 2000
    block_region_y1 = 2048
    ###########################################
    iraf.imred(Stdout=1)
    iraf.echelle(Stdout=1)
    dumout = []
    for l, reff in enumerate(reffitlist):
        dummy = np.ones((2048, 2048))*1000.
        dummy[block_region_y0:block_region_y1, :] = 0.
        dummy_fits_name = 'dummy_formsk_' + str(time.time()) + '.fits'
        fits.writeto(dummy_fits_name, dummy)

        temporary_mask_pl = 'tmp_' + str(time.time())
        iraf.apmask(input=dummy_fits_name, output=temporary_mask_pl,
                    references=reff,
                    interactive='no', find='no', recenter='no', resize='no',
                    edit='no', trace='no', fittrace='no',
                    mask='yes', line='INDEF', nsum=10, buffer=0.0)
        dummy_fits_name_out = 'dummy_formsk_out' + str(time.time()) + '.fits'
        if len(reffitlist) == 1:
            iraf.imar(dummy_fits_name, '*', temporary_mask_pl + '.pl',
                      outputf + '.fits')
        else:
            iraf.imar(dummy_fits_name, '*', temporary_mask_pl + '.pl',
                      dummy_fits_name_out)

            dumout.append(dummy_fits_name_out)

    if len(reffitlist) > 1:
        dummy_fits_name = dumout[0]
        for l in range(len(reffitlist)-1):
            if l == len(reffitlist)-2:
                dummy_fits_name_swap = outputf + '.fits'
            else:
                dummy_fits_name_swap = 'dummy_foradd_' + \
                    str(time.time()) + '.fits'
            iraf.imar(dummy_fits_name, '+', dumout[l+1],
                      dummy_fits_name_swap)
            dummy_fits_name = dummy_fits_name_swap

    dumdirect = []
    if len(directfitslist) > 0:
        #        from scipy.signal import medfilt
        for l, fitsimg in enumerate(directfitslist):
            crit = directfits_criterion[l]
            fits_file = str(fitsimg)
            hdulist = fits.open(fits_file+'.fits')
            hdulist[0].header
            img = hdulist[0].data
            if len(directrot) == len((directfitslist)):
                if directrot[l]:
                    print('Rotation applied: '+str(fits_file))
                    img = img[::-1, ::-1]
            dmimg = np.zeros(np.shape(img))
            dmimg[img > crit] = 1000.0
            dumdirect.append(dmimg)
    try:
        msk_tmp = fits.open(outputf + '.fits')[0].data
    except:
        print('Failed. Probably no ref file or directory mismatch.')
        sys.exit()
    msk_tmp = np.array(msk_tmp)
    if len(directfitslist) > 0:
        for l, fitsimg in enumerate(directfitslist):
            msk_tmp = msk_tmp+dumdirect[l]
    msk_tmp[msk_tmp > 1000] = 1000.0
    iraf.imdel(outputf + '.fits')
    fits.writeto(outputf + '.fits', msk_tmp)

    return outputf

# if __name__ == "__main__":
#
#   makemsk(sys.argv[1])
