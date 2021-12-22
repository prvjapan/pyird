import pytest
from pyird.utils.getinfo import get_radec


def test_get_radec():
    ra,dec=get_radec("HR8799")
    print(ra,dec)
    delta=((ra - 346.86964833333326)**2+(dec - 21.134250555555553)**2)
    assert delta==0.0

if __name__=="__main__":
    test_get_radec()
