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
target_name = "G196-3B"

combine_flg = False
#-------------------------#

target_1d = irdstream.Stream1D(streamid=target_name, 
                               rawdir=rawdir, 
                               anadir=anadir_1d, 
                               fitsid=list(range(target_id[0], target_id[-1], 2)),
                               prefix=target_prefix, # w or nw
                               extension=f"_m{mmf[-1]}", # fiber
                               inst="IRD",
                               band=band)

# wavelength recalibration with comb
comb_master_path = rawdir / f"wcomb_master_{band}_m{mmf[-1]}.dat"
df_recalib = target_1d.recalibrate_wavelength_with_comb(comb_master_path, fiber=mmf, n_poly=6)

# airglow mask
if mmf == "mmf2":
    target_1d.mask_airglow()

# combine spectra
if combine_flg:
    df_merge = target_1d.spec_combine(method="mean")

    savedir_merge = anadir_1d/f"{target_1d.prefix}{target_name}_{band}{target_1d.extension}.dat"
    df_merge.to_csv(savedir_merge, **target_1d.tocsvargs)

# concatenate spectra (YJ + H)
if band == "h":
    target_1d.spec_concat_yjh()