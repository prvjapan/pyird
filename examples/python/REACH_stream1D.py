from pyird.utils import irdstream
import pathlib

basedir = pathlib.Path('~/pyird/data/20211110_REACH/').expanduser()

datadir=basedir/'reduc'
anadir=basedir/'reduc'
target_dat = irdstream.Stream1D("HR37635",datadir,anadir,prefix='nw', extension='_m2')
target_dat.fitsid=[53205]
target_dat.band='h'
#target_dat.fitsid_increment()
target_dat.date='20211110'
target_dat.remove_fringe()
