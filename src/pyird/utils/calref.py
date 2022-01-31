from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.signal import medfilt
#from specutils import Spectrum1D
import pandas as pd
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='wav cal reference')
    parser.add_argument('-f', nargs=1, help='1D multispec fits', type=str)
    parser.add_argument('-o', nargs=1, help='order', type=int)
    parser.add_argument('-c', nargs=1, help='channel list', type=str)

    parser.add_argument(
        '-s', default=[0], nargs=1, help='offset value, for SMF/YJ probably set 1. Depending on position. 1-51 or 1-52', type=int)

    args = parser.parse_args()
    offset = args.s[0]
    hd = fits.open(args.f[0])
    dat = (hd[0].data)
    ssw = False
    if np.shape(dat)[0] == 2048:
        print('WARNING')
        print('The Input is probably 2D image. Use 1D fits.')
        ssw = True

    if args.c:
        chanfile = args.c[0]
    else:
        if np.shape(dat)[0] == 21:
            chanfile = '../../data/channel_H.list'
            norder = 21
            print('H band')
        elif np.shape(dat)[0] > 45:
            chanfile = '../../data/channel_YJ.list'
            norder = np.shape(dat)[0]
            print('YJ band')

        else:
            sys.exit('1D multispec is abnormal. Set chanfile (-c) by yourself.')

    pdat = pd.read_csv(chanfile, delimiter=',')

    if args.o:
        j = int(args.o[0])-1
        l = j+1
    else:
        j = 0
        l = norder

    line = ['dashed', 'dotted']
    for i in range(j, l):
        porder = pdat['ORDER']
        pchan = pdat['CHANNEL']
        wav = pdat['WAVELENGTH']
        st = pdat['ST']

        med = medfilt(dat[i, :], kernel_size=3)

        fig = plt.figure(figsize=(15, 12))
        ax = fig.add_subplot(211)
        if ssw:
            pos = (np.min(med)+np.max(med)*1.05)/2.0
            plt.text(0.0, pos, 'WARNING: INPUT IS PROBABLY 2D FITS.',
                     fontsize=50, color='red')
        ax.plot(med, '.')
        ax.plot(dat[i, :], alpha=0.3)
        plt.ylim(np.min(med), np.max(med)*1.05)
        s = 0
        for k, po in enumerate(porder):
            if po == i+1-offset:
                s = s+1
                plt.axvline(pchan[k], color='green',
                            ls=line[int(st[k])], alpha=0.3)
                plt.text(pchan[k], np.max(med)*(0.8-0.05*s),
                         str(wav[k]), color='green')

        plt.title('Order='+str(i+1))
        plt.ylabel('median filtered')

        ax2 = fig.add_subplot(212)
        ax2.plot(dat[i, :])
        s = 0
        for k, po in enumerate(porder):
            if po == i+1-offset:
                s = s+1
                plt.axvline(pchan[k], color='green',
                            ls=line[int(st[k])], alpha=0.3)
                plt.text(pchan[k], np.max(dat[i, :]) *
                         (0.8-0.05*s), str(wav[k]), color='green')
        plt.ylabel('spectrum')
        plt.xlabel('CHANNEL')
        plt.show()
        fig.clear()
