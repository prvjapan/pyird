from pyird.utils import irdstream
from pyird.spec.wavcal import wavcal_thar
from pyird.image.bias import bias_subtract_image
from pyird.image.hotpix import identify_hotpix
from pyird.spec.rsdmat import multiorder_to_rsd
from pyird.image.oned_extract import flatten
from pyird.plot.detector import corrected_detector
from pyird.image.pattern_model import median_XY_profile
from pyird.image.trace_function import trace_legendre
from pyird.image.mask import trace
from pyird.io.iraf_trace import read_trace_file
from pyird.plot.order import plot_refthar
from pyird.io.read_linelist import read_linelist
import astropy.io.fits as pyf
from scipy.signal import medfilt
import pathlib
import numpy as np
import pandas as pd
import pkg_resources

import argparse
import matplotlib.pyplot as plt

def imcorrect(date,im):
    """read pattern removal.
    """
    # image for calibration
    calim = np.copy(im)
    # read pattern removal
    pathC = 'path/to/apwhite_mask'  #(pkg_resources.resource_filename('pyird', 'data/samples/aprefC'))
    #pathC = pathlib.Path('/Volumes/IRD/IRD/observation_BD/data/hband/%d/database/apwhite_mask'%(date))
    y0, interp_function, xmin, xmax, coeff = read_trace_file(pathC) #([pathC, path_c])
    mask = trace(im, trace_legendre, y0, xmin, xmax, coeff)
    calim[mask] = np.nan
    calim[hotpix_mask] = np.nan
    model_im = median_XY_profile(calim,show=plot)
    corrected_im = im-model_im
    #corrected_detector(im, model_im, corrected_im)
    return corrected_im

def im_to_rsd(date,corrected_im,mode='mmf2'):
    """One Dimensionalization.
    """
    # trace
    if mode=='mmf2':
        pathA = 'path/to/apwhite_mmf2_r'  #(pkg_resources.resource_filename('pyird', 'data/samples/aprefB'))
        #pathA = pathlib.Path('/Volumes/IRD/IRD/observation_BD/data/hband/%d/database/apwhite_mmf2_r'%(date))
    elif mode=='mmf1':
        pathA = 'path/to/apwhite_mmf1_r'
        #pathA = pathlib.Path('/Volumes/IRD/IRD/observation_BD/data/hband/%d/database/apwhite_mmf1_r'%(date))
    y0, interp_function, xmin, xmax, coeff = read_trace_file(pathA)
    # flatten
    rawspec, pixcoord = flatten(
        corrected_im, trace_legendre, y0, xmin, xmax, coeff)
    # rsd
    rsd = multiorder_to_rsd(rawspec, pixcoord)
    return rsd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='copy of REARCH.py')
    parser.add_argument('-d', type=int, nargs=1, default=[20210321], help='date')
    parser.add_argument('--dark', type=int, nargs=1,default=[41955], help='dark frame')
    parser.add_argument('--thar', type=int, nargs=2,default=[14893, 14942], help='thar_star frame; start, stop')
    parser.add_argument('--step', type=int, nargs=1,default=[1],help='thar_star step of frame number')
    parser.add_argument('--mode', nargs=1,default=['mmf2'], help='mmf1 (lfc fiber) or mmf2 (star fiber)')
    parser.add_argument('--plot',default=False, help='plot')

    args = parser.parse_args()
    date = args.d[0]
    fdark = args.dark
    fthar = list(range(args.thar[0],args.thar[1],args.step[0]))
    mode = arge.mode[0]
    plot = args.plot

    # hotpixel mask
    datadir = pathlib.Path('/Users/yuikasagi/IRD/PhDwork/pyird/data/dark/')
    anadir = pathlib.Path('/Users/yuikasagi/IRD/PhDwork/pyird/data/dark/')
    dark = irdstream.Stream2D('targets', datadir, anadir)
    dark.fitsid = fdark
    for data in dark.rawpath:
        im = pyf.open(str(data))[0].data
    im_subbias = bias_subtract_image(im)
    hotpix_mask, obj = identify_hotpix(im_subbias)

    # thar images
    datadir = pathlib.Path('/Volumes/IRD/IRD/observation_BD/data/hband/%d/raw/'%(date))
    anadir = pathlib.Path('/Users/yuikasagi/IRD/PhDwork/pyird/data/%d/'%(date))
    if args.step[0]==1:
        rawtag='IRDAD000'
    elif args.step[0]==2:
        rawtag='IRDA000'
    thars1 = irdstream.Stream2D('thars1', datadir, anadir, rawtag=rawtag)
    thars1.fitsid = fthar  # THARSTAR

    # Load an images
    corrected_im_all = []
    for datapath in thars1.rawpath:
        im = pyf.open(str(datapath))[0].data
        corrected_im = imcorrect(date,im)
        corrected_im_all.append(corrected_im)

    # median combine of thar spectra
    thar1_comb = np.nanmedian(corrected_im_all,axis=0) #np.nansum(corrected_im_all,axis=0)#
    rsd = im_to_rsd(date,thar1_comb,mode=mode)

    dat = rsd.T

    # wavlength calibration
    wavsol2, data2 = wavcal_thar(dat,maxiter=30,stdlim=0.001)
    print('features = ', len(data2.ravel()[data2.ravel()!=0]))
    plot_refthar(wavsol2, data2, 21)
    np.savetxt('./wavsol2.txt',wavsol2)
