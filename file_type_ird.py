#!/usr/bin/env python
from astropy.io import fits
import sys
import argparse
import numpy as np

def checktype(filelist):
    flat=[]
    expflat=[]
    comp=[]
    bias=[]
    expbias=[]
    obs=[]
    obsname=[]
    flatt=[]
    for i in range(0,len(filelist)):
#        file=args.f[i]
        file=filelist[i]
        hduread=fits.open(file)
        name=hduread[0].header["Object"]
        expt=hduread[0].header["EXPTIME"]
        slt=hduread[0].header["SLT-LEN"]

        if name=="FLAT":
            flat.append(file)
            expflat.append(expt)
            flatt.append(slt)
        elif name=="COMPARISON":
            comp.append(file)        
        elif name=="BIAS":
            bias.append(file)
            expbias.append(expt)
        else:
            obs.append(file)
            obsname.append(name)
        
    return     obs,bias,comp,flat,expflat,expbias,obsname,flatt

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='check file type')
    parser.add_argument('-f', nargs="+", required=True, help='fits')
    args = parser.parse_args()    

    obs,bias,comp,flat,expflat,expbias,obsname,flatt=checktype(args.f)

    print((str(len(obs))+" OBJECTS files found:\n", obs, "\n"))
    print(("OBJECT NAME=\n",obsname, "\n"))
    print((str(len(bias))+" BIAS files found:\n", bias, "\n"))
#print "EXPOSURE TIME=", expbias, "\n"
    print((str(len(comp))+" COMPARISON files found:\n", comp, "\n"))
    print((str(len(flat))+" FLAT files found:", flat, "\n"))
#print "EXPOSURE TIME=",expflat

    if min(flatt) < max(flatt):
        mask = (np.array(flatt)==min(flatt))
        flatt=np.array(flatt)
        flat=np.array(flat)
        print("########################################")
        print(("THESE "+str(len(flat[mask]))+" FILES ARE ORDER TRACE FLAT ?"))
        print((flat[mask]))
        print(("SLT-LEN=",flatt[mask],"mm"))
        print("These has shorter slit-length (SLT-LEN) than others.") 
        print("Check them.")
        print("########################################")

