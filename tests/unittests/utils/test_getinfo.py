import pytest
from pyird.utils.getinfo import get_radec


def test_get_radec():
    ra, dec = get_radec('HR8799 A')
    print(ra, dec)
    delta = ((ra - 346.86964875)**2+(dec - 21.134252777777778)**2)
    assert pytest.approx(delta) == 0.0


if __name__ == '__main__':
    test_get_radec()
