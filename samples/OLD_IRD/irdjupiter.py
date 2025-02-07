"""Raw data analysis of Jupiter/miniIRD using pyird
"""
import matplotlib.pyplot as plt
import pathlib
import numpy as np
from pyird.utils import irdstream

basedir = pathlib.Path("~/jovispec/analysis/data/H").expanduser()

# Specifies directories
datadir = basedir / "raw/"
anadir = basedir / "reduc/"

# flat and tracing apertures
flat = irdstream.Stream2D("flat", datadir, anadir, rawtag="2017011", fitsid=[8140833, 8140945, 8141330], rotate=True, inverse=True, detector_artifact=True)
trace_mmf = flat.aptrace(cutrow=580, nap=21)
trace_mask = trace_mmf.mask()
## note that we do not apply trace_mmf.choose_mmf1_aperture nor trace_mmf.choose_mmf1_aperture
## because there is only a single mmf on the detector

# dark and bias subtraction
from pyird.image.bias import bias_subtract_image
from pyird.image.hotpix import identify_hotpix_sigclip
dark = irdstream.Stream2D("dark", datadir, anadir, rawtag="2017011", fitsid=[8133836], rotate=True, inverse=True)
median_image = dark.immedian()
im_subbias = bias_subtract_image(median_image)
hotpix_mask = identify_hotpix_sigclip(im_subbias)

# generates the normalized flat
flat.trace = trace_mmf
flat.clean_pattern(trace_mask=trace_mask,extin='', extout='_cp', hotpix_mask=hotpix_mask)
flat.imcomb = True # median combine
flat.flatten(hotpix_mask=hotpix_mask)
df_flatn = flat.apnormalize()

# Th-Ar wavelength calibration
thar=irdstream.Stream2D("thar",datadir,anadir,rawtag="2017011",fitsid=[8143952, 8144142, 8144237, 8144415],rotate=True, inverse=True)
thar.trace = trace_mmf
thar.clean_pattern(trace_mask=trace_mask,extin='', extout='_cp', hotpix_mask=hotpix_mask)
thar.calibrate_wavelength()

# target reduction
target = irdstream.Stream2D("target", datadir, anadir, rawtag="2017011", fitsid=[8034325], rotate=True, inverse=True)
target.info = True
target.trace = trace_mmf
target.clean_pattern(trace_mask=trace_mask,extin='', extout='_cp', hotpix_mask=hotpix_mask)

# flat fielding
target.apext_flatfield(df_flatn,hotpix_mask=hotpix_mask)

# simliar to the IRAF task "dispcor"
target.dispcor(master_path=thar.anadir,extin='_flnhp')
