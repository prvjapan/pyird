from pyird import fitsset 
from pyird import irdstream
from pyird import makemask
import sys
import fitsids_IRDBD

targets, flat_mmf, thars1, thars2=fitsids_IRDBD.IRDBD_Stream2D(mode="YJ")

datadir,anadir=fitsids_IRDBD.directories(mode="YJ")

#Bias Correction 
targets.remove_bias()
flat_mmf.remove_bias()
thars1.remove_bias()
thars2.remove_bias()

#Combined Flat
#combined_flat=flat_mmf.imcombine("combined_flat")
combined_flat=fitsset.FitsSet("combined_flat",anadir)

#Read-Noise Correction
maskfits=fitsset.FitsSet("mask_from_refMMF1_refMMF2_Dhotpixel_20190412_YJ",anadir)

targets.rm_readnoise(maskfits)
flat_mmf.rm_readnoise(maskfits)

#RBN combine flat
#flat_mmf.extension="_rb"
#combined_flat_rb=flat_mmf.imcombine("combined_flat_rb")

#Combined Thorium-Argon
#combined_thars1=thars1.imcombine("combined_thars1")
#combined_thars2=thars2.imcombine("combined_thars2")


Mext=""
extin="rb"
apref_target=fitsset.FitsSet("refA"+str(Mext),anadir) #aperture reference
wavref_target=fitsset.FitsSet("apthars1",anadir) #wave reference
combined_flat.apall(apref_target,wavref_target,extout=extin+"_1d"+str(Mext),extwout="flat_center_"+extin+"_1dw"+str(Mext))

Mext=""
extin="_rbo"
apref_target=fitsset.FitsSet("refA"+str(Mext),anadir) #aperture reference
wavref_target=fitsset.FitsSet("apthars1",anadir) #wave reference
targets.extract1D(apref_target,wavref_target,extin=extin,extout="_"+extin+"_1d"+str(Mext),extwout="_"+extin+"_1dw"+str(Mext))

