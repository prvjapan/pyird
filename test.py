from pyird import fitsset 
import pathlib

datadir=pathlib.Path("/media/kawahara/kingyo/ird/data/ird/IRDseminar/YJ")
anadir=pathlib.Path("/media/kawahara/kingyo/ird/ana/IRDSeminar/yband")

# READ FitsStream
targets=fitsset.FitsStream("targets",datadir/"star",anadir)
targets.fitsid=list(range(3126,3186,2))+list(range(3194,3200,2))

flat_mmf1=fitsset.FitsStream("flat_mmf1",datadir/"flat"/"mmf1",anadir)
flat_mmf1.fitsid=list(range(3502,3562,2))

flat_mmf2=fitsset.FitsStream("flat_mmf2",datadir/"flat"/"mmf2",anadir)
flat_mmf2.fitsid=list(range(3436,3496,2))

comparison=fitsset.FitsStream("comp",datadir/"comp"/"20180627",anadir,"IRDBD0000")
comparison.fitsid=list(range(3686,3718,2))

if False:
    #Bias correction
    targets.remove_bias()
    flat_mmf1.remove_bias()
    flat_mmf2.remove_bias()
    comparison.remove_bias()
    
if False:
    #Combine for making mask
    flat_mmf1.combine_rb()
    flat_mmf2.combine_rb()

print("Use ",flat_mmf1.rbcombined.path(string=True)," to make mask file.")
    
#Define Mask file
maskpath=pathlib.Path("/media/kawahara/kingyo/ird/data/destripe_codes/combined")
maskfits=fitsset.FitsSet("mask_from_white_mmf12_201906_yj",maskpath)

if False:
    #Read-Noise correction
    targets.rm_readnoise(maskfits)
    flat_mmf1.rm_readnoise(maskfits)
    flat_mmf2.rm_readnoise(maskfits)
    comparison.rm_readnoise(maskfits)

if False:
    #combined flat
    flat_mmf1.combine_rbn()
    flat_mmf2.combine_rbn()
