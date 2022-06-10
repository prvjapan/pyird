import pytest
import pkg_resources
from pyird.spec.continuum import comb_norm

def test_normalize1D():
    import numpy as np
    flat_path = (pkg_resources.resource_filename('pyird', 'data/samples/wflat_m2_20210317.dat'))
    df_interp = comb_norm(flat_path,flat_path)
    nflux = df_interp['nflux']
    assert round(np.mean(nflux))==1

if __name__ == '__main__':
    test_normalize1D()
