import processRN
    
if __name__ == "__main__":
    import pandas as pd
    import argparse
    import tqdm
    parser = argparse.ArgumentParser(description='wrap hirano-kun cpp code')
    parser.add_argument('-f', nargs=1, default=["raw/IRDBD_rb.lst"], help='list', type=str)

    args = parser.parse_args()    
    numlist=pd.read_csv(args.f[0],names=["NUM"])
    for num in tqdm.tqdm(numlist["NUM"]):
        processRN.wrap_processRN(filen="raw/IRDBD0000"+str(num)+"_rb.fits",\
                       filemmf="../../destripe_codes/combined/mask_from_white_mmf12_201906_yj.fits",\
                       fitsout="obj/r90"+str(num)+".fits")

        
