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
        self.info=False
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

    def remove_bias(self, rot=None, method='reference', hotpix_img=None):
        if self.info:
            print('remove_bias: files=')
            print(self.rawpath)
        if rot == 'r':
            print('180 degree rotation applied.')
        #IRD_bias_sube.main(self.anadir,self.rawpath, method,rot,hotpix_im = hotpix_img)
        self.fitsdir = self.anadir
        self.extension = '_rb'

    def clean_pattern(self, extout, trace_path_list, hotpix_mask=None, extin=""):
        """

        Args:
           extout: output extension
           path list of trace files
           extin: input extension
           hotpix_mask: hot pixel mask

        """
        from pyird.io.iraf_trace import read_trace_file
        from pyird.image.pattern_model import median_XY_profile
        from pyird.image.trace_function import trace_legendre
        from pyird.image.mask import trace
        
        currentdir = os.getcwd()
        os.chdir(str(self.anadir))
        if self.info:
            print('clean_pattern: output extension=',extout)

        extin_noexist, extout_noexist = self.check_existence(extin, extout)
        for i, fitsid in enumerate(tqdm.tqdm(extin_noexist)):
            filen=extin_noexist[i]
            hdu=pyf.open(filen)[0]
            im = hdu.data
            header = hdu.header
            calim = np.copy(im) # image for calibration
            y0, interp_function, xmin, xmax, coeff = read_trace_file(trace_path_list)
            mask = trace(im, trace_legendre, y0, xmin, xmax, coeff)
            calim[mask] = np.nan
            if hotpix_mask is not None:
                calim[hotpix_mask] = np.nan
            model_im = median_XY_profile(calim,show=False)
            corrected_im = im-model_im
            hdu = pyf.PrimaryHDU(corrected_im, header)
            hdulist = pyf.HDUList([hdu])
            hdulist.writeto(extout_noexist[i], overwrite=True)

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
        """check files do not exist or not

        Args:
           extin: extension of input files
           extout: extension of output files

        Returns:
           input file name w/ no exist
           output file name w/ no exist

        """
        
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
                if self.info:
                    print("Ignore ",str(ext[i].name),"->",str(extfi.name))
                skip = skip+1

        if skip > 1:
            print('Skipped '+str(skip)+' files because they already exists.')

        return ext_noexist, extf_noexist
