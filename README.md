<img src=documents/figures/pyird_logo.png width=70%>

PyIRD
===========
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/prvjapan/pyird/blob/develop/LICENSE)
[![Version](https://img.shields.io/badge/version-v1.0.0-blue?logo=github)](https://github.com/prvjapan/pyird/releases/tag/v1.0.0)
[![PyPI Version](https://img.shields.io/pypi/v/pyird)](https://pypi.org/project/pyird/)
[![Run pytest](https://github.com/prvjapan/pyird/actions/workflows/main.yml/badge.svg)](https://github.com/prvjapan/pyird/actions/workflows/main.yml)

`PyIRD` is a Python-based pipeline for reducing spectroscopic data obtained with [IRD](https://ird.mtk.nao.ac.jp/IRDpub/index_tmp.html) and [REACH](https://secondearths.sakura.ne.jp/reach/) on the Subaru Telescope. 
It is designed to process raw images into one-dimensional spectra in a semi-automatic manner. 
Unlike traditional methods, it does not rely on IRAF (Tody et al. [1986](https://ui.adsabs.harvard.edu/abs/1986SPIE..627..733T/abstract), [1993](https://ui.adsabs.harvard.edu/abs/1993ASPC...52..173T/abstract)), a software traditionally used for astronomical data reduction. This approach simplifies the workflow while maintaining efficiency and accuracy.
Additionally, the pipeline includes an updated method for removing readout noise patterns from the detector, enabling efficient extraction of spectra even for faint targets such as brown dwarfs.


Install
------------------

```
pip install pyird
```

or

```
git clone https://github.com/prvjapan/pyird.git
cd pyird
python setup.py install
```

Examples
------------------

`PyIRD` is designed to perform data reduction semi-automatically by following a general workflow for high-dispersion spectroscopic data reduction (e.g., readout noise subtraction, flat fielding, aperture extraction, wavelength calibration, and normalization).

See the following examples how to create 1D spectra from raw data.
- [pyird/examples/python/IRD_stream.py](https://github.com/prvjapan/pyird/blob/master/examples/python/IRD_stream.py) -- for IRD data; read [the docs](https://secondearths.sakura.ne.jp/pyird/tutorials/IRD_stream.html) for the detailed explanation.
- [pyird/examples/python/REACH_stream.py](https://github.com/prvjapan/pyird/blob/master/examples/python/REACH_stream.py) -- for REACH data; basically the same as `IRD_stream.py`, but with the variable `inst` set to `REACH` instead.

The raw data for `IRD_stream.py` can be downloaded from here. 


License
------------------------------
`PyIRD` is publicly available under the MIT license. For developers, please read [CONTRIBUTING.md](https://github.com/prvjapan/pyird/blob/develop/CONTRIBUTING.md).