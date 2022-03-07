import pytest
from pyird.utils.irdstream import Stream2D
import pathlib

def test_path():
    datadir = pathlib.Path('/home/kawahara/pyird/data/dark/')
    anadir = pathlib.Path('/home/kawahara/pyird/data/dark/')
    s2d = Stream2D('targets', datadir, anadir)
    s2d.fitsid = [41018]
    refs=pathlib.Path('/home/kawahara/pyird/data/dark/')/pathlib.Path("IRDA00041018.fits")
    assert s2d.path()[0]==refs


if __name__ == '__main__':
    test_path()
