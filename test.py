from pyird import fitsset 
import pathlib

fset=fitsset.FitsSet("IRDSeminar")
fset.anadir=pathlib.Path("/media/kawahara/kingyo/ird/ana/yband")

datadir=pathlib.Path("/media/kawahara/kingyo/ird/data/ird/IRDseminar/YJ")

# set identifiers
fset.targets={"num":list(range(3126,3186,2))+list(range(3194,3200,2)),"dir":datadir/"star","tag":"IRDA0000"}
fset.flat_mmf1={"num":list(range(3502,3562,2)),"dir":datadir/"flat"/"mmf1","tag":"IRDA0000"}
fset.flat_mmf2={"num":list(range(3436,3496,2)),"dir":datadir/"flat"/"mmf2","tag":"IRDA0000"}
fset.comparison={"num":list(range(3686,3718,2)),"dir":datadir/"comp"/"20180627","tag":"IRDBD0000"}
fset.check_ready()

