"""File stream for IRD analysis."""

from pyird.utils.fitsset import FitsSet
from pyird.utils import directory_util
import astropy.io.fits as pyf
import numpy as np
import tqdm
import os
__all__ = ['Stream1D', 'Stream2D']


class Stream1D(FitsSet):
    def __init__(self, streamid, rawdir, anadir, fitsid=None, rawtag='IRDA000', extension=''):
        """initialization
        Args:
           streamid: ID for stream
           rawdir: directory where the raw data are
           anadir: directory in which the processed file will put
        """
        super(Stream1D, self).__init__(rawtag, rawdir, extension=extension)
        self.streamid = streamid
        self.anadir = anadir
        self.unlock = False
        if fitsid is not None:
            print("fitsid:",fitsid)
            self.fitsid=fitsid
        else:
            print("No fitsid yet.")

    @property
    def fitsid(self):
        return self._fitsid

    @fitsid.setter
    def fitsid(self, fitsid):
        self._fitsid = fitsid
        self.rawpath = self.path(string=False, check=True)

    def fitsid_increment(self):
        """increase fitsid +1
        """
        for i in range(0, len(self.fitsid)):
            self.fitsid[i] = self.fitsid[i]+1
        self.rawpath = self.path(string=False, check=True)

    def fitsid_decrement(self):
        """decrease fitsid +1
        """
        for i in range(0, len(self.fitsid)):
            self.fitsid[i] = self.fitsid[i]-1
        self.rawpath = self.path(string=False, check=True)

    def extpath(self, extension, string=False, check=True):
        """decrease fitsid +1

        Args:
           extension: extension
        
        Returns:
           path array of fits files w/ extension

        """
        f = self.fitsdir
        e = self.extension
        self.fitsdir = self.anadir
        self.extension = extension
        path_ = self.path(string, check)
        self.fitsdir = f
        self.extension = e
        return path_


class Stream2D(FitsSet):
    def __init__(self, streamid, rawdir, anadir, fitsid=None, rawtag='IRDA000', extension=''):
        """initialization
        Args:
           streamid: ID for stream
           rawdir: directory where the raw data are
           anadir: directory in which the processed file will put
           fitsid: fitsid

        """
        super(Stream2D, self).__init__(rawtag, rawdir, extension='')        
        self.streamid = streamid
        self.rawdir = rawdir
        self.anadir = anadir
        self.unlock = False
        if fitsid is not None:
            print("fitsid:",fitsid)
            self.fitsid=fitsid
        else:
            print("No fitsid yet.")
            

    @property
    def fitsid(self):
        return self._fitsid

    @fitsid.setter
    def fitsid(self, fitsid):
        self._fitsid = fitsid
        self.rawpath = self.path(string=False, check=True)

    def fitsid_increment(self):
        """Increase fits id +1."""
        for i in range(0, len(self.fitsid)):
            self.fitsid[i] = self.fitsid[i]+1
        self.rawpath = self.path(string=False, check=True)

    def fitsid_decrement(self):
        """Decrease fits id +1."""
        for i in range(0, len(self.fitsid)):
            self.fitsid[i] = self.fitsid[i]-1
        self.rawpath = self.path(string=False, check=True)

    def load_fitsset(self):
        """Load fitsset and make imcube.

        Returns:
           imcube
        """
        imcube = []
        for data in tqdm.tqdm(self.rawpath):
            im = pyf.open(str(data))[0].data
            imcube.append(im)
        return np.array(imcube)

    ############################################################################################
    def extpath(self, extension, string=False, check=True):
        """decrease fitsid +1

        Args:
           extension: extension
        
        Returns:
           path array of fits files w/ extension

        """
        f = self.fitsdir
        e = self.extension
        self.fitsdir = self.anadir
        self.extension = extension
        path_ = self.path(string, check)
        self.fitsdir = f
        self.extension = e
        return path_

    def extclean(self, extension):
        """Clean i.e. remove fits files if exists.
        
        Args:
            extension: extension of which files to be removed

        """
        import os
        if extension == '':
            print('extclean cannot clean w/o extension fits.')
            return

        for i in range(len(self.path())):
            if self.path(check=False)[i].exists():
                os.remove(self.extpath(extension, check=False, string=True)[i])
                print('rm old '+self.extpath(extension,
                      check=False, string=True)[i])

    def remove_bias(self, rot=None, method='reference', hotpix_img=None, info=False):
        if info:
            print('remove_bias: files=')
            print(self.rawpath)
        if rot == 'r':
            print('180 degree rotation applied.')
        #IRD_bias_sube.main(self.anadir,self.rawpath, method,rot,hotpix_im = hotpix_img)
        self.fitsdir = self.anadir
        self.extension = '_rb'

    def rm_readnoise(self, maskfits, extin='_rb', extout='_rbo', info=False):
        if info:
            print('rm_readnoise: fits=')
            print(maskfits)
        currentdir = os.getcwd()
        os.chdir(str(self.anadir))

        ext_noexist, extf_noexist = self.check_existence(extin, extout)
        for i, fitsid in enumerate(tqdm.tqdm(ext_noexist)):
            maskfits.path()[0]
            # processRN.wrap_kawahara_processRN(filen=ext_noexist[i],filemask=filemask,fitsout=extf_noexist[i])
        self.fitsdir = self.anadir
        self.extension = extout

        os.chdir(currentdir)

    def flatfielding1D(self):
        return 

    def extract1D(self):
        """extract 1D specctra
        """
        return
    
    def check_existence(self, extin, extout):
        extf = self.extpath(extout, string=False, check=False)
        ext = self.extpath(extin, string=False, check=False)
        extf_noexist = []
        ext_noexist = []
        skip = 0
        for i, extfi in enumerate(extf):
            if not extfi.exists():
                extf_noexist.append(str(extfi.name))
                ext_noexist.append(str(ext[i].name))
            else:
                print(str(extfi.name), str(ext[i].name))
                skip = skip+1

        if skip > 1:
            print('Skipped '+str(skip)+' files because they already exists.')

        return ext_noexist, extf_noexist
