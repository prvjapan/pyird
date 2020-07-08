from pyird import fitsset 
import pathlib

datadir=pathlib.Path("/media/kawahara/kingyo/ird/data/ird/IRDseminar/YJ")
anadir=pathlib.Path("/media/kawahara/kingyo/ird/ana/IRDSeminar/yband")

# READ FitsStream
targets=fitsset.FitsStream("targets",datadir/"star",anadir)
targets.fitsid=list(range(3126,3186,2))+list(range(3194,3200,2))
targets.remove_bias()
targets.rb_combine()

flat_mmf1=fitsset.FitsStream("flat_mmf1",datadir/"flat"/"mmf1",anadir)
flat_mmf1.fitsid=list(range(3502,3562,2))
flat_mmf1.remove_bias()
flat_mmf1.rb_combine()

flat_mmf2=fitsset.FitsStream("flat_mmf2",datadir/"flat"/"mmf2",anadir)
flat_mmf2.fitsid=list(range(3436,3496,2))
flat_mmf2.remove_bias()
flat_mmf2.rb_combine()

comparison=fitsset.FitsStream("comp",datadir/"comp"/"20180627",anadir,"IRDBD0000")
comparison.fitsid=list(range(3686,3718,2))
comparison.remove_bias()

