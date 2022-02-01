"""Aperture."""

if __name__ == '__main__':
    import numpy as np
    from pyird.image.operator import imcombine
    from pyird.utils import detect_peaks, irdstream
    import pathlib
    import matplotlib.pyplot as plt

    mode = 'H'
    datadir = pathlib.Path(
        '/media/kawahara/kingyo_bkup/kingyo/IRD_G196_3B/20210317_cal/')
    anadir = pathlib.Path('/home/kawahara/pyird/data/flat/')

    flat_mmf = irdstream.Stream2D('flat_mmf', datadir, anadir)
    flat_mmf.fitsid =\
        list(range(41704, 41904, 2))
    if mode == 'H':
        flat_mmf.fitsid_increment()
        Norder = 20*2+1
    elif mode == 'YJ':
        Norder = 52*2

    imcube = flat_mmf.load_fitsset()
    cflat = imcombine(imcube)

    pk = detect_peaks.detect_peaks(cflat[1024, :])
    ind = np.argsort(cflat[1024, pk])[::-1][:Norder]
    pk = (pk[ind])
    x = np.array(range(0, 2048))
    plt.plot(x, cflat[1024, :])
    plt.plot(x[pk], cflat[1024, pk], 'o')

    plt.show()
