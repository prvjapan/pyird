import pytest
from pyird.utils.irdstream import Stream2D
import pathlib

def test_path_1():
    datadir = pathlib.Path('~/pyird/data/dark/').expanduser()
    anadir = pathlib.Path('~/pyird/data/dark/').expanduser()
    s2d = Stream2D('targets', datadir, anadir, fitsid = [41018])
    refs=pathlib.Path('~/pyird/data/dark/').expanduser()/pathlib.Path("IRDA00041018.fits")
    assert s2d.path()[0]==refs

def test_path_2():
    datadir = pathlib.Path('~/pyird/data/dark/').expanduser()
    anadir = pathlib.Path('~/pyird/data/dark/').expanduser()
    s2d = Stream2D('targets', datadir, anadir)
    s2d.fitsid=[41018]
    refs=pathlib.Path('~/pyird/data/dark/').expanduser()/pathlib.Path("IRDA00041018.fits")
    assert s2d.path()[0]==refs

if __name__ == '__main__':
    test_path_1()
    test_path_2()
