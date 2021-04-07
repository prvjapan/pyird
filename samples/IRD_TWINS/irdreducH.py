from pyird import fitsset 
from pyird import irdstream
from pyird import makemask
from pyird import maskcheck as mc
import sys
import fitsids


targets, flat_mmf, thars1=fitsids.IRDBD_Stream2D(mode="H")
datadir,anadir,fdatadir=fitsids.directories(mode="H")

#Bias Correction with 180 deg rotation
targets.remove_bias("r")
thars1.remove_bias("r")
flat_mmf.remove_bias("r")

#Combined Flat
combined_flat=flat_mmf.imcombine("combined_flat")
combined_flat=fitsset.FitsSet("combined_flat",anadir)

#Notes:
#refA should be made using apall for combined_flat
#refC = all orders A+B with wider width (apall)

#MASK
makemask.makemsk(anadir,["refC"],["hotpixel_20190412_H"],[0.99],[True])
maskfits=fitsset.FitsSet("mask_from_refC_Dhotpixel_20190412_H",anadir)

mc.plotmask(maskfits, combined_flat,vmin=0.0,vmax=2000)

#Removing read noise
targets.rm_readnoise(maskfits)
flat_mmf.rm_readnoise(maskfits)

#RBN combine flat
combined_flat_rb=flat_mmf.imcombine("combined_flat_rb")

#Combined Thorium-Argon
combined_thars1=thars1.imcombine("combined_thars1")

#Apall 1D Flat
combined_flat=fitsset.FitsSet("combined_flat_rb",anadir)
apref_target=fitsset.FitsSet("refA",anadir) #aperture reference
wavref_target=fitsset.FitsSet("combined_thars1.ec",anadir) #wave reference
combined_flat.apall(apref_target,wavref_target)

#Apall 1D
extin="_rbo"
apref_target=fitsset.FitsSet("refA",anadir) #aperture reference
wavref_target=fitsset.FitsSet("combined_thars1.ec",anadir) #wave reference
targets.extract1D(apref_target,wavref_target,extin=extin,extout="_"+extin+"_1d",extwout="_"+extin+"_1dw")

