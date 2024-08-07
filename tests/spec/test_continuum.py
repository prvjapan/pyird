import pytest
import pkg_resources
from pyird.spec.normalize import SpectrumNormalizer

def test_normalize1D():
    import numpy as np
    spectrum_normalizer = SpectrumNormalizer()
    flat_path = (pkg_resources.resource_filename('pyird', 'data/samples/wflat_m2_20210317.dat'))
    _, df_interp = spectrum_normalizer.combine_normalize(flat_path,flat_path,blaze=False)
    nflux = df_interp['nflux']
    assert round(np.mean(nflux))==1

def test_normalize1D_blaze():
    import numpy as np
    spectrum_normalizer = SpectrumNormalizer()
    flat_path = (pkg_resources.resource_filename('pyird', 'data/samples/wblaze_h_m2_20210317.dat'))
    _, df_interp = spectrum_normalizer.combine_normalize(flat_path,flat_path,blaze=True)
    nflux = df_interp['nflux']
    assert round(np.mean(nflux))==1

if __name__ == '__main__':
    test_normalize1D()
