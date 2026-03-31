from pyird.utils import irdstream
import pathlib

#--------SETTINGS--------#
basedir = pathlib.Path('~/pyird/data/20210317/').expanduser()

band = 'h' #'h' or 'y'
mmf = 'mmf2' #'mmf1' (comb fiber) or 'mmf2' (star fiber)

rawdir = basedir/'reduc/'
anadir_1d = basedir/'reduc_1d/'

# last 5 digits of FITS file numbers: [start, end file number]
target_id = [41510, 41511] # target image
target_prefix = "nw"
#-------------------------#


target_1d = irdstream.Stream1D("target_1d", 
                               rawdir=rawdir, 
                               anadir=anadir_1d, 
                               fitsid=list(range(target_id[0], target_id[-1], 2)),
                               prefix=target_prefix, # w or nw
                               extension=f"_m{mmf[-1]}", # fiber
                               inst="IRD",
                               band=band)

comb_master_path = rawdir / f"wcomb_master_{band}_m{mmf[-1]}.dat"
df_recalib = target_1d.recalibrate_wavelength_with_comb(comb_master_path, fiber=mmf, n_poly=6)

# df_recalib.to_csv(anadir_1d/f"pcomb_master_{band}_m{mmf[-1]}.dat", index=False, sep=' ', header=None)