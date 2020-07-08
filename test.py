from pyird import fitsset 
from pyird import IRD_bias_sube
import pathlib
from pyraf import iraf

fset=fitsset.FitsSet("IRDSeminar")
datadir=pathlib.Path("/media/kawahara/kingyo/ird/data/ird/IRDseminar/YJ")
fset.anadir=pathlib.Path("/media/kawahara/kingyo/ird/ana/IRDSeminar/yband")

# set identifiers
fset.targets={"id":list(range(3126,3186,2))+list(range(3194,3200,2)),"dir":datadir/"star","tag":"IRDA0000"}
fset.flat_mmf1={"id":list(range(3502,3562,2)),"dir":datadir/"flat"/"mmf1","tag":"IRDA0000"}
fset.flat_mmf2={"id":list(range(3436,3496,2)),"dir":datadir/"flat"/"mmf2","tag":"IRDA0000"}
fset.comparison={"id":list(range(3686,3718,2)),"dir":datadir/"comp"/"20180627","tag":"IRDBD0000"}

# bias
if False:
    method = 'reference'
    hotpix_img = None
    IRD_bias_sube.main(fset.anadir,fset.rawfitslist(), method, hotpix_im = hotpix_img)

# combine mmf1 and mmf2
iraf.imcombine(input=fitsset.at_rblist(fset.flat_mmf1,fset.anadir),output='mmf1.fits',combine="median")
iraf.imcombine(input=fitsset.at_rblist(fset.flat_mmf2,fset.anadir),output='mmf2.fits',combine="median")

