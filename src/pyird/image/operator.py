"""image operator for imcube.

* imcube has a shape of [N x Ny x Nx], where N is the number of images.
"""
import numpy as np


def imcombine(imcube, mode='median'):
    """combine images to a image.

    Args:
       imcube: imcube
       mode: median or mean
    Return:
       image
    """
    if mode == 'median':
        mask = imcube != imcube
        if np.shape(imcube[mask])[0] > 0:  # nan check
            im = np.nanmedian(imcube, axis=0)
        else:
            im = np.median(imcube, axis=0)
    elif mode == 'mean':
        mask = imcube != imcube
        if np.shape(imcube[mask])[0] > 0:  # nan check
            im = np.nanmean(imcube, axis=0)
        else:
            im = np.mean(imcube, axis=0)
    else:
        print('No mode in imcombine.')
    return im


if __name__ == '__main__':
    import numpy as np
    from pyird.utils import irdstream
    import pathlib
    import matplotlib.pyplot as plt

    mode = 'YJ'
    datadir = pathlib.Path(
        '/media/kawahara/kingyo_bkup/kingyo/IRD_G196_3B/20210317_cal/')
    anadir = pathlib.Path('/home/kawahara/pyird/data/flat/')

    flat_mmf = irdstream.Stream2D('flat_mmf', datadir, anadir)
    flat_mmf.fitsid =\
        list(range(41704, 41904, 2))
    if mode == 'H':
        flat_mmf.fitsid_increment()
    imcube = flat_mmf.load_fitsset()
    cflat = imcombine(imcube)
    plt.plot(cflat[1024, :])
    plt.show()
#    def imcombine(self,tag,combine="median"):
#        combined_fitsset=FitsSet(tag,self.fitsdir)
#        combined_fitsset.clean()
#        iraf.imcombine(input=self.at_list(listname=tag),output=combined_fitsset.path(check=False)[0],combine=combine) #
#        return combined_fitsset

#    combined_flat=flat_mmf.imcombine("combined_flat")
