
def wrap_kawahara_processRN(filen,filemask,fitsout):
    #OEGP method
    import astropy.io.fits as pyf
    import numpy as np
    from pyird import readnoise as rn
    from pyird import large_scaled

    hdulist=pyf.open(filemask)
    hdu=hdulist[0]
    header=hdu.header
    mask=np.array(hdu.data) 
    mask[mask==0]=0
    mask[mask>0]=1
    mask=np.array(mask,dtype=np.bool)
    #input img
    hdulist=pyf.open(filen)
    hdu=hdulist[0]
    header=hdu.header
    img=np.array(hdu.data) 
    imgorig=np.copy(img)
    #print("# of Masked pix is "+str(len(img[mask])))
    img[mask]=None
    
    check=large_scaled.check_mask_filling(mask,Ncor=64)
    if check:
        LSD=large_scaled.get_LSD(img,gpscale=64,Ncor=64,sigma=0.001)
    else:
        sys.exit("exit")
    ###############################
    img=img-LSD
    recimg=rn.RNestimate_OEGP(img)
    ###############################

    cimg=imgorig-recimg-LSD
    
    hdu = pyf.PrimaryHDU(cimg, header)
    hdulist = pyf.HDUList([hdu])
    hdulist.writeto(fitsout,overwrite=True)


def wrap_hirano_processRN(filen="raw/IRDA00003502_rb.fits",filemmf="../../destripe_codes/combined/mask_from_white_mmf12_201906_yj.fits",fitsout="obj/r03502.fits"):
    
    import astropy.io.fits as pyf
    import numpy as np
    import process_RNe
    
    hdulist=pyf.open(filen)
    hdu=hdulist[0]
    adata=np.array(hdu.data, dtype=np.float64) 
    header=hdu.header
    
    hdulist=pyf.open(filemmf)
    hdu=hdulist[0]
    bdata=np.array(hdu.data, dtype=np.float64) 
    
    outa=process_RNe.processRN((adata.flatten()),bdata.flatten())
    
    hdu = pyf.PrimaryHDU(outa.T, header)
    hdulist = pyf.HDUList([hdu])
    hdulist.writeto(fitsout,overwrite=True)
    
if __name__ == "__main__":
    import pandas as pd
    import argparse
    import tqdm
    parser = argparse.ArgumentParser(description='wrap hirano-kun cpp code')
    parser.add_argument('-f', nargs=1, default=["raw/all_rb.lst"], help='list', type=str)
    print("If ended as core dump, set ulimit -s unlimited.")
    args = parser.parse_args()    
    numlist=pd.read_csv(args.f[0],names=["NUM"])
    for num in tqdm.tqdm(numlist["NUM"]):
        wrap_hirano_processRN(filen="raw/IRDA0000"+str(num)+"_rb.fits",\
                       filemmf="../../destripe_codes/combined/mask_from_white_mmf12_201906_yj.fits",\
                       fitsout="obj/r0"+str(num)+".fits")

        
