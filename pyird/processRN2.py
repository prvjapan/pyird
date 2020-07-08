def wrap_processRN(filen="raw/IRDA00003502_rb.fits",filemmf="../../destripe_codes/combined/mask_from_white_mmf12_201906_yj.fits",fitsout="obj/r03502.fits"):
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
    parser.add_argument('-f', nargs=1, default=["raw/IRDBD_rb.lst"], help='list', type=str)

    args = parser.parse_args()    
    numlist=pd.read_csv(args.f[0],names=["NUM"])
    for num in tqdm.tqdm(numlist["NUM"]):
        wrap_processRN(filen="raw/IRDBD0000"+str(num)+"_rb.fits",\
                       filemmf="../../destripe_codes/combined/mask_from_white_mmf12_201906_yj.fits",\
                       fitsout="obj/r90"+str(num)+".fits")

        
