import pytest
from pyird.utils.irdstream import Stream2D
import pathlib

def test_path_1():
    basedir = pathlib.Path(__file__).parent.parent.parent
    datadir = basedir / 'data/dark/'
    anadir = basedir / 'data/dark/'
    s2d = Stream2D('targets', datadir, anadir, fitsid = [41018])
    refs = basedir / 'data/dark/IRDA00041018.fits'
    assert s2d.path()[0]==refs

def test_path_2():
    basedir = pathlib.Path(__file__).parent.parent.parent
    datadir = basedir / 'data/dark/'
    anadir = basedir / 'data/dark/'
    s2d = Stream2D('targets', datadir, anadir)
    s2d.fitsid=[41018]
    refs = basedir / 'data/dark/IRDA00041018.fits'
    assert s2d.path()[0]==refs

if __name__ == '__main__':
    test_path_1()
    test_path_2()
