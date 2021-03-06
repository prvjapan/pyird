"""File stream for IRD analysis."""

from pyird.utils.fitsset import FitsSet
from pyird.utils import directory_util
from pyird.image.trace_function import trace_legendre
from pyird.image.aptrace import aptrace
from pyird.utils.aperture import TraceAperture
import astropy.io.fits as pyf
import numpy as np
import tqdm
import os
__all__ = ['Stream1D', 'Stream2D']


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
        self.info = False
        self.imcomb = False
        if fitsid is not None:
            print('fitsid:', fitsid)
            self.fitsid = fitsid
        else:
            print('No fitsid yet.')

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
        """path to file with extension.

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

    def clean_pattern(self, trace_path_list=None, hotpix_mask=None, extout='_cp', extin=None):
        """

        Args:
           trace_path_list: path list of trace files
           hotpix_mask: hot pixel mask
           extout: output extension
           extin: input extension

        """
        from pyird.image.pattern_model import median_XY_profile
        from pyird.image.mask import trace_from_iraf_trace_file

        currentdir = os.getcwd()
        os.chdir(str(self.anadir))
        if self.info:
            print('clean_pattern: output extension=', extout)

        if trace_path_list is None:
            mask = self.trace.mask()
        else:
            mask = trace_from_iraf_trace_file(trace_path_list, mask_shape=np.shape(im))

        if extin is None:
            extin = self.extension

        extin_noexist, extout_noexist = self.check_existence(extin, extout)

        for i, fitsid in enumerate(tqdm.tqdm(extin_noexist)):
            filen = extin_noexist[i]
            hdu = pyf.open(filen)[0]
            im = hdu.data
            header = hdu.header
            calim = np.copy(im)  # image for calibration
            calim[mask] = np.nan
            if hotpix_mask is not None:
                calim[hotpix_mask] = np.nan
            model_im = median_XY_profile(calim, show=False)
            corrected_im = im-model_im
            hdu = pyf.PrimaryHDU(corrected_im, header)
            hdulist = pyf.HDUList([hdu])
            hdulist.writeto(extout_noexist[i], overwrite=True)

        self.fitsdir = self.anadir
        self.extension = extout
        os.chdir(currentdir)

    def flatten(self, trace_path=None, extout='_fl', extin=None):
        """
        Args:
           trace_path: trace file to be used in flatten
           extout: output extension
           extin: input extension

        """
        from pyird.image.oned_extract import flatten
        from pyird.image.trace_function import trace_legendre
        from pyird.image.mask import trace
        from pyird.io.iraf_trace import read_trace_file
        from pyird.spec.rsdmat import multiorder_to_rsd

        currentdir = os.getcwd()
        os.chdir(str(self.anadir))
        if self.info:
            print('flatten: output extension=', extout)

        if trace_path is None:
            y0, xmin, xmax, coeff = self.trace.y0, self.trace.xmin, self.trace.xmax, self.trace.coeff
        else:
            y0, interp_function, xmin, xmax, coeff = read_trace_file(trace_path)

        if extin is None:
            extin = self.extension

        if not self.imcomb:
            extin_noexist, extout_noexist = self.check_existence(extin, extout)
            for i, fitsid in enumerate(tqdm.tqdm(extin_noexist)):
                filen = extin_noexist[i]
                hdu = pyf.open(filen)[0]
                im = hdu.data
                header = hdu.header
                rawspec, pixcoord = flatten(
                    im, trace_legendre, y0, xmin, xmax, coeff)
                rsd = multiorder_to_rsd(rawspec, pixcoord)
                hdux = pyf.PrimaryHDU(rsd, header)
                hdulist = pyf.HDUList([hdux])
                hdulist.writeto(extout_noexist[i], overwrite=True)
        else:
            median_image=self.immedian()
            hdu = pyf.open(self.path()[0])[0]
            header = hdu.header
            rawspec, pixcoord = flatten(
                median_image, trace_legendre, y0, xmin, xmax, coeff)
            rsd = multiorder_to_rsd(rawspec, pixcoord)
            hdux = pyf.PrimaryHDU(rsd, header)
            hdulist = pyf.HDUList([hdux])
            save_path = self.anadir/('%s_mmf12.fits'%(self.streamid))
            hdulist.writeto(save_path, overwrite=True)
        self.fitsdir = self.anadir
        self.extension = extout
        os.chdir(currentdir)

    def immedian(self):
        """take image median.

        Return:
           median image
        """
        imall = []
        for path in tqdm.tqdm(self.path(string=True, check=True)):
            imall.append(pyf.open(path)[0].data)
        imall = np.array(imall)
        # np.nansum(corrected_im_all,axis=0)#
        median_image = np.nanmedian(imall, axis=0)
        return median_image

    def calibrate_wavlength(self, trace_file_path=None, master='thar_master.fits' ,maxiter=30, stdlim=0.001, npix=2048):
        """wavelength calibration usgin Th-Ar.

        Args:
           trace_file_path: path to the trace file
           master: master file for the wavelength calibrated ThAr file
           maxiter: maximum number of iterations
           stdlim: When the std of fitting residuals reaches this value, the iteration is terminated.
           npix: number of pixels

        """
        from pyird.spec.wavcal import wavcal_thar
        from pyird.io.iraf_trace import read_trace_file
        from pyird.image.oned_extract import flatten
        from pyird.image.trace_function import trace_legendre
        from pyird.spec.rsdmat import multiorder_to_rsd

        currentdir = os.getcwd()
        os.chdir(str(self.anadir))

        median_image = self.immedian()

        if trace_file_path is None:
            y0, xmin, xmax, coeff = self.trace.y0, self.trace.xmin, self.trace.xmax, self.trace.coeff
        else:
            y0, interp_function, xmin, xmax, coeff = read_trace_file(trace_file_path)

        master_path = master
        if not os.path.exists(master_path):
            filen = self.path()[0] #header of the first file
            hdu = pyf.open(filen)[0]
            im = hdu.data
            header = hdu.header
            nord = len(y0)
            rawspec, pixcoord = flatten(
                median_image, trace_legendre, y0, xmin, xmax, coeff)
            rsd = multiorder_to_rsd(rawspec, pixcoord)
            print(np.shape(rsd.T)[0])
            wavsol, data = wavcal_thar(rsd.T, maxiter=maxiter, stdlim=stdlim)
            wavsol_2d = wavsol.reshape((npix,nord))
            hdux = pyf.PrimaryHDU(wavsol_2d, header)
            hdulist = pyf.HDUList([hdux])
            hdulist.writeto(master_path, overwrite=True)

        os.chdir(currentdir)

    def flatfielding1D(self):
        return

    def extract1D(self):
        """extract 1D specctra."""
        return

    def check_existence(self, extin, extout):
        """check files do not exist or not.

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
                    print('Ignore ', str(ext[i].name), '->', str(extfi.name))
                skip = skip+1

        if skip > 1:
            print('Skipped '+str(skip)+' files because they already exists.')

        return ext_noexist, extf_noexist

    def aptrace(self,cutrow = 1000,nap=42):
        """extract aperture of trace from a median image of current fitsset

        Args:
            cutrow: cutting criterion
            nap: number of apertures, nap = 42 ## 42 for H band, 102 for YJ band

        Returns:
            TraceAperture instance

        """
        flatmedian=self.immedian()
        if nap==42:
            flatmedian = flatmedian[::-1,::-1]
            y0, xmin, xmax, coeff = aptrace(flatmedian,cutrow,nap)

        return TraceAperture(trace_legendre, y0, xmin, xmax, coeff)

    def dispcor(self, input_ext='_fl', prefix='w', master='thar_master.fits',master_path=None):
        """dispersion correct and resample spectra

        Args:
            input_ext: extension of input files
            prefix: prefix for output files
            master: master file for the wavelength calibrated ThAr file
            master_path: path of the directory containing the master ThAr file

        """
        from pyird.plot.showspec import show_wavcal_spectrum
        def mkwspec(spec_m2,reference,save_path):
            import pandas as pd
            wspec = pd.DataFrame([],columns=['wav','order','flux'])
            for i in range(len(reference[0])):
                wav = reference[:,i]
                order = np.ones(len(wav))
                order[:] = i+1
                data_order = [wav,order,spec_m2[:,i]]
                df_order = pd.DataFrame(data_order,index=['wav','order','flux']).T
                wspec = pd.concat([wspec,df_order])
            wspec.to_csv(save_path,header=False,index=False,sep=' ')
            return wspec

        if master_path==None:
            master_path = self.anadir.joinpath('..','thar').resolve()/master

        if not self.imcomb:
            inputs = self.extpath(input_ext,string=False, check=True)
            for j,input in enumerate(inputs):
                hdu = pyf.open(input)[0]
                spec_m12 = hdu.data
                spec_m2 = spec_m12#[:,::2] # choose mmf2 (star fiber)

                hdu = pyf.open(master_path)[0]
                reference = hdu.data

                id = self.fitsid[j]
                save_path = self.anadir/('%s%d_m2.dat'%(prefix,id)) ##for mmf2
                wspec = mkwspec(spec_m2,reference,save_path)
                if self.info:
                    print('dispcor: output spectrum= %s%d_m2.dat'%(prefix,id))
                #plot
                show_wavcal_spectrum(wspec,alpha=0.5)
        else:
            hdu = pyf.open(self.anadir/('%s_mmf12.fits'%(self.streamid)))[0]
            spec_m12 = hdu.data
            spec_m2 = spec_m12#[:,::2] # choose mmf2 (star fiber)

            hdu = pyf.open(master_path)[0]
            reference = hdu.data

            save_path = self.anadir/('%s%s_m2.dat'%(prefix,self.streamid))
            wspec = mkwspec(spec_m2,reference,save_path)
            if self.info:
                print('dispcor: output spectrum= %s%s_m2.dat'%(prefix,self.streamid))
            #plot
            show_wavcal_spectrum(wspec,alpha=0.5)

    def normalize1D(self,flatid='flat'):
        """combine orders and normalize spectrum

        Args:
            flatid: streamid for flat data

        """
        from pyird.spec.continuum import comb_norm
        from pyird.plot.showspec import show_wavcal_spectrum
        for id in self.fitsid:
            wfile = self.anadir/('w%d_m2.dat'%(id))
            flatfile = self.anadir.joinpath('..','flat').resolve()/('w%s_m2.dat'%(flatid))
            df_interp = comb_norm(wfile,flatfile)
            df_save = df_interp[['wav','nflux']]
            save_path = self.anadir/('ncw%d_m2.dat'%(id))
            df_save.to_csv(save_path,header=False,index=False,sep=' ')
            if self.info:
                print('normalize1D: output normalized 1D spectrum= ncw%d_m2.dat'%(id))
            #plot
            show_wavcal_spectrum(df_save,alpha=0.5)
