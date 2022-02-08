from astropy.io import fits
import numpy as np
from scipy.signal import medfilt
import pandas as pd

if __name__ == '__main__':
    # Load a Th-Ar spectrum
    tharfile = '/Users/kasagiyui/IRD/PhDwork/pyird/data/20210317/thar_mmf2_a_h.fits'
    hd = fits.open(tharfile)
    dat = (hd[0].data)

    # Load a channel list
    if np.shape(dat)[0] == 21:
        chanfile = '../../../data/channel_H.list'
        norder = 21
        print('H band')
    elif np.shape(dat)[0] > 45:
        chanfile = '../../../data/channel_YJ.list'
        norder = np.shape(dat)[0]
        print('YJ band')

    pdat = pd.read_csv(chanfile, delimiter=',')
    j = 0
    l = norder
    offset = 0

    # allocate the line positions in pixcoord
    pdat_new = pd.DataFrame([],columns=['ORDER','CHANNEL','WAVELENGTH'])
    for i in range(j, l):
        porder = pdat['ORDER']
        pchan = pdat['CHANNEL']
        wav = pdat['WAVELENGTH']
        #st = pdat['ST']

        med = medfilt(dat[i, :], kernel_size=3)

        for k, po in enumerate(porder):
            if po == i+1-offset:
                paround = 10  ## search area
                plocal_med = max(med[pchan[k]-paround:pchan[k]+paround])
                pchanlocal_med = np.where(med==plocal_med)[0]
                pchanlocal_dat = np.where(dat[i,:]==max(dat[i,:][pchanlocal_med]))[0][0]
                data = [i+1-offset,pchanlocal_dat,wav[k]]
                df_tmp = pd.DataFrame([data],columns=pdat_new.columns)
                pdat_new = pd.concat([pdat_new,df_tmp],ignore_index=True)
    #print(pdat_new)
