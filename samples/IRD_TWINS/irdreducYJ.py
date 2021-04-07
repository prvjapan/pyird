from pyird import fitsset 
from pyird import irdstream
from pyird import makemask
import sys
import fitsids


targets, flat_mmf, thars1=fitsids.IRDBD_Stream2D(mode="YJ")
datadir,anadir,fdatadir=fitsids.directories(mode="YJ")

#Bias Correction with 180 deg rotation
targets.remove_bias()
thars1.remove_bias()
flat_mmf.remove_bias()

#Combined Flat
#combined_flat=flat_mmf.imcombine("combined_flat")
combined_flat=fitsset.FitsSet("combined_flat",anadir)

#Read-Noise Correction
#mask should be prepared by mask.py 
#maskfits=fitsset.FitsSet("mask_from_refC_refMMF1",anadir) #aperture reference

#maskfits=fitsset.FitsSet("mask_uniform",anadir)
#dark.extclean("_rbn")
#dark.rm_readnoise(maskfits)

#refMMF1,2 are mask for readnoise removal, do not use them to extract spectra
maskfits=fitsset.FitsSet("mask_from_refMMF1_refMMF2_Dhotpixel_20190412_YJ",anadir)

targets.rm_readnoise(maskfits)
flat_mmf.rm_readnoise(maskfits)

#RBN combine flat
#flat_mmf.extension="_rb"
#combined_flat_rb=flat_mmf.imcombine("combined_flat_rb")


#Combined Thorium-Argon
#combined_thars1=thars1.imcombine("combined_thars1")

#combined_thars2=thars2.imcombine("combined_thars2")

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

