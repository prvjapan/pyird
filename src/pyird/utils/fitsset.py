"""Set of fits files."""

import sys
__all__ = ['FitsSet']


class FitsSet(object):

    def __init__(self, tag, fitsdir, extension=''):
        self.tag = tag
        self.extension = extension
        self.fitsdir = fitsdir
        self.fitsid = None
        self.unlock = True  # unlock for cleanning (i.e. remove)
        self.not_ignore_warning = True

    def data(self, indices=None):
        """get header and data in an array form
        Args:
            indices: indices of fitsset, if None, 0 is used

        Returns:
            array of data
        """
        import astropy.io.fits as pyf
        import numpy as np
        if indices == None:
            hdulist = pyf.open(np.array(self.path(string=True))[0])
            hdu = hdulist[0]
            return np.array(hdu.data)

        data_arr = []
        for filen in np.array(self.path(string=True))[np.array(indices)]:
            hdulist = pyf.open(filen)
            hdu = hdulist[0]
            data_arr.append(np.array(hdu.data))
        return data_arr

    def path(self, string=False, check=True):
        """get path array.

        Returns:
            array of paths
        """
        fitsarray = []
        if self.fitsid is not None:
            for num in self.fitsid:
                d = self.fitsdir/(self.tag+str(num)+self.extension+'.fits')
                if check and not d.exists():
                    print(d, 'does not exist.')
                    if self.not_ignore_warning:
                        print(
                            'If you wanna proceed despite this warning, set self.not_ignore_warning = False in your fitsset.')
                        sys.exit('Error.')
                else:
                    if string:
                        fitsarray.append(str(d))
                    else:
                        fitsarray.append(d)
        return fitsarray

    def clean(self):
        """Clean i.e. remove fits files if exists."""
        import os
        if self.unlock:
            for i in range(len(self.path())):
                if self.path(check=False)[i].exists():
                    os.remove(self.path(check=False, string=True)[i])
                    print('rm old '+self.path(check=False, string=True)[i])

    def at_list(self, listname='tmp'):
        """make at list used in IRAF."""

        listname = listname+'.list'
        atlistname = str(self.fitsdir/listname)
        f = open(atlistname, mode='w')
        for fits in self.path(string=True):
            f.write(fits+'\n')
        f.close()
        atlistname = '@'+atlistname
        return atlistname
