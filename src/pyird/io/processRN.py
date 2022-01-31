
def wrap_kawahara_processRN(filen, filemask, fitsout):
    # OEGP method
    import astropy.io.fits as pyf
    import numpy as np
    from pyird import readnoise as rn
    from pyird import large_scaled

    hdulist = pyf.open(filemask)
    hdu = hdulist[0]
    header = hdu.header
    mask = np.array(hdu.data)
    mask[mask == 0] = 0
    mask[mask > 0] = 1
    mask = np.array(mask, dtype=np.bool)
    # input img
    hdulist = pyf.open(filen)
    hdu = hdulist[0]
    header = hdu.header
    img = np.array(hdu.data)
    imgorig = np.copy(img)
    # print("# of Masked pix is "+str(len(img[mask])))
    img[mask] = None

    # Large Scale Distribution
    check = large_scaled.check_mask_filling(mask, Ncor=64)
    if check:
        LSD = large_scaled.get_LSD(img, gpscale=1024, Ncor=64, sigma=0.001)
    else:
        sys.exit('exit')
    ###############################
    img = img-LSD
    recimg = rn.RNestimate_OEGP(img)
    ###############################

    cimg = imgorig-recimg-LSD

    hdu = pyf.PrimaryHDU(cimg, header)
    hdulist = pyf.HDUList([hdu])
    hdulist.writeto(fitsout, overwrite=True)
