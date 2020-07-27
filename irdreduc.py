from pyird import fitsset 
from pyird import irdstream
import pathlib
import sys

datadir=pathlib.Path("/media/kawahara/kingyo/ird/data/ird/IRDseminar/YJ")
anadir=pathlib.Path("/media/kawahara/kingyo/ird/ana/IRDSeminar/yband")

#Define Stream2D from raw data
targets=irdstream.Stream2D("targets",datadir/"star",anadir)
targets.fitsid=list(range(3126,3186,2))+list(range(3194,3200,2))
flat_mmf1=irdstream.Stream2D("flat_mmf1",datadir/"flat"/"mmf1",anadir)
flat_mmf1.fitsid=list(range(3502,3562,2))
flat_mmf2=irdstream.Stream2D("flat_mmf2",datadir/"flat"/"mmf2",anadir)
flat_mmf2.fitsid=list(range(3436,3496,2))
comb=irdstream.Stream2D("comp",datadir/"comp"/"20180627",anadir,"IRDBD0000")
comb.fitsid=list(range(3686,3693,1))+list(range(3714,3716,1))
thar=irdstream.Stream2D("comp",datadir/"comp"/"20180627",anadir,"IRDBD0000")
thar.fitsid=list(range(3693,3714,1))

#Bias Correction
targets.remove_bias()
flat_mmf1.remove_bias()
flat_mmf2.remove_bias()
comb.remove_bias()
thar.remove_bias()

#Combine for making mask
combined_flat_for_mask=flat_mmf1.imcombine("combined_flat_for_mask")
maskfits=combined_flat_for_mask.apall("combined_flat_for_mask_ap")
    
#Define Mask file
maskpath=pathlib.Path("/media/kawahara/kingyo/ird/data/destripe_codes/combined")
maskfits=fitsset.FitsSet("mask_from_white_mmf12_201906_yj",maskpath)

#Read-Noise Correction
targets.rm_readnoise(maskfits)
flat_mmf1.rm_readnoise(maskfits)
flat_mmf2.rm_readnoise(maskfits)
comb.rm_readnoise(maskfits)
thar.rm_readnoise(maskfits)

#Combined Flat
combined_flat_mmf1=flat_mmf1.imcombine("combined_flat_mmf1")
combined_flat_mmf2=flat_mmf2.imcombine("combined_flat_mmf2")

#Combined Thorium-Argon
combined_thar=thar.imcombine("combined_thar")
