"""File stream for IRD analysis."""

from pyird.utils.fitsset import FitsSet
from pyird.utils.datset import DatSet
from pyird.utils import directory_util
from pyird.image.trace_function import trace_legendre
from pyird.image.aptrace import aptrace
from pyird.utils.aperture import TraceAperture
import astropy.io.fits as pyf
import numpy as np
import pandas as pd
import tqdm
import os

import pkg_resources

__all__ = ['Stream1D', 'Stream2D']


class Stream2D(FitsSet):
    def __init__(self, streamid, rawdir, anadir, fitsid=None, rawtag='IRDA000', extension='', inst='IRD'):
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
        self.inst = inst
        self.tocsvargs = {'header':False,'index':False,'sep':' '}
        if fitsid is not None:
            print('fitsid:', fitsid)
            self.fitsid = fitsid
            if (rawtag=='IRDA000' and fitsid[0]%2==0) or (rawtag=='IRDBD000'):
                self.band = 'y'
            else:
                self.band = 'h'
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
        self.band = 'h'

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

    def clean_pattern(self, trace_mask=None, hotpix_mask=None, extout='_cp', extin=None):
        """

        Args:
           trace_mask: trace mask (from iraf, use image.mask.trace_from_iraf_trace_file)
           hotpix_mask: hot pixel mask
           extout: output extension
           extin: input extension

        """
        from pyird.image.pattern_model import median_XY_profile
        #from pyird.image.mask import trace_from_iraf_trace_file

        currentdir = os.getcwd()
        os.chdir(str(self.anadir))
        if self.info:
            print('clean_pattern: output extension=', extout)

        if trace_mask is None:
            trace_mask = self.trace.mask()

        if extin is None:
            extin = self.extension

        extin_noexist, extout_noexist = self.check_existence(extin, extout)
        # print('extin_noexist: ',extin_noexist)
        # print('extout_noexist: ',extout_noexist)
        for i, fitsid in enumerate(tqdm.tqdm(extin_noexist)):
            filen = self.rawdir/extin_noexist[i]
            hdu = pyf.open(filen)[0]
            im = hdu.data
            header = hdu.header
            if self.band == 'h':
                calim = np.copy(im)  # image for calibration
                calim[trace_mask] = np.nan
                if hotpix_mask is not None:
                    calim[hotpix_mask] = np.nan
            elif self.band == 'y':
                calim = np.copy(im)  # image for calibration
                trace_mask_r = trace_mask.reshape(1, int(trace_mask.size))
                trace_mask_r = trace_mask_r[0][::-1].reshape(im.shape[0], im.shape[1]) # rotate mask matrix 180 degree
                calim[trace_mask_r] = np.nan
                if hotpix_mask is not None:
                    calim[hotpix_mask] = np.nan

            model_im = median_XY_profile(calim)
            corrected_im = im - model_im
            hdu = pyf.PrimaryHDU(corrected_im, header)
            hdulist = pyf.HDUList([hdu])
            hdulist.writeto(extout_noexist[i], overwrite=True)
            # np.savetxt('im.txt',im,delimiter=',')
            # np.savetxt('model_im.txt',model_im,delimiter=',')
            # np.savetxt('calim.txt', calim, delimiter=',')

        self.fitsdir = self.anadir
        self.extension = extout
        os.chdir(currentdir)

    def flatten(self, trace_path=None, extout='_fl', extin=None, hotpix_mask=None, width=None):
        """
        Args:
           trace_path: trace file to be used in flatten
           extout: output extension
           extin: input extension
           hotpix_mask: hotpix masked spectrum ('extout' to be automatically '_hp')
           width: list of aperture widths ([width_start,width_end])

        """
        from pyird.image.oned_extract import flatten
        from pyird.image.trace_function import trace_legendre
        from pyird.image.mask import trace
        from pyird.io.iraf_trace import read_trace_file
        from pyird.spec.rsdmat import multiorder_to_rsd
        from pyird.image.hotpix import apply_hotpixel_mask

        currentdir = os.getcwd()
        os.chdir(str(self.anadir))

        if trace_path is None:
            y0, xmin, xmax, coeff = self.trace.y0, self.trace.xmin, self.trace.xmax, self.trace.coeff
            mmf = self.trace.mmf
            width = self.trace.width
        else:
            y0, interp_function, xmin, xmax, coeff = read_trace_file(trace_path)

        if extin is None:
            extin = self.extension

        extout = extout + '_' + mmf
        if not hotpix_mask is None:
            extout = '_hp_' + mmf

        if not self.imcomb:
            extin_noexist, extout_noexist = self.check_existence(extin, extout)
            for i, fitsid in enumerate(tqdm.tqdm(extin_noexist)):
                filen = self.anadir/extin_noexist[i]
                hdu = pyf.open(filen)[0]
                im = hdu.data
                header = hdu.header
                rawspec, pixcoord, _, _, _, _ = flatten(
                    im, trace_legendre, y0, xmin, xmax, coeff, inst=self.inst, width=width)
                rsd = multiorder_to_rsd(rawspec, pixcoord)
                if not hotpix_mask is None:
                    save_path = self.anadir/('hotpix_%s_%s.fits'%(self.band,self.trace.mmf))
                    rsd = apply_hotpixel_mask(hotpix_mask, rsd, y0, xmin, xmax, coeff, save_path=save_path)
                hdux = pyf.PrimaryHDU(rsd, header)
                hdulist = pyf.HDUList([hdux])
                hdulist.writeto(extout_noexist[i], overwrite=True)
        else:
            save_path = self.anadir/('%s_%s_%s.fits'%(self.streamid,self.band,mmf))
            if not os.path.exists(save_path):
                median_image=self.immedian()
                if not hotpix_mask is None:
                    median_image=median_image*(~hotpix_mask)
                hdu = pyf.open(self.path()[0])[0]
                header = hdu.header
                rawspec, pixcoord, _, _, _, _ = flatten(
                    median_image, trace_legendre, y0, xmin, xmax, coeff, inst=self.inst, width=width)
                rsd = multiorder_to_rsd(rawspec, pixcoord)
                hdux = pyf.PrimaryHDU(rsd, header)
                hdulist = pyf.HDUList([hdux])
                hdulist.writeto(save_path, overwrite=True)

        if self.info:
            print('flatten (+ hotpix mask): output extension=', extout)

        self.fitsdir = self.anadir
        self.extension = extout
        os.chdir(currentdir)

    def flatten_check(self, trace_path=None, extin=None, wavcal_path=None, hotpix_mask=None):
        """
        Args:
           trace_path: trace file to be used in flatten
           extin: input extension
            wavcal_path: path to master file of the wavelength calibration
           hotpix_mask: hotpix mask
        """
        from pyird.image.oned_extract import flatten
        from pyird.image.trace_function import trace_legendre
        from pyird.io.iraf_trace import read_trace_file
        from pyird.spec.rsdmat import multiorder_to_rsd
        from pyird.image.hotpix import apply_hotpixel_mask

        def im_to_rsd(im, hotpix_mask=None,wavcal_path=None):
            rawspec, pixcoord, rotim, tl, iys_plot, iye_plot = flatten(
                im, trace_legendre, y0, xmin, xmax, coeff, self.inst)
            mask_shape = (2048,2048)
            mask = np.zeros(mask_shape)
            for i in range(len(y0)):
                mask[i,xmin[i]:xmax[i]+1] = tl[i]
            rsd = multiorder_to_rsd(rawspec, pixcoord)
            if not hotpix_mask is None:
                rsd = apply_hotpixel_mask(hotpix_mask, rsd, y0, xmin, xmax, coeff)
            if not wavcal_path is None:
                hdu = pyf.open(wavcal_path)
                wav = hdu[0].data
                hdu.close()
            else:
                wav = []
            return rsd,wav,mask,pixcoord,rotim,iys_plot,iye_plot

        currentdir = os.getcwd()
        os.chdir(str(self.anadir))

        global y0, xmin, xmax, coeff
        if trace_path is None:
            y0, xmin, xmax, coeff = self.trace.y0, self.trace.xmin, self.trace.xmax, self.trace.coeff
            mmf = self.trace.mmf
        else:
            y0, interp_function, xmin, xmax, coeff = read_trace_file(trace_path)

        if extin is None:
            extin = self.extension

        if not self.imcomb:
            print(self.extpath(extin, string=False, check=False))
            #extin_noexist, extout_noexist = [str(self.extpath(extin, string=False, check=False)[0].name)], [str(self.extpath(extout, string=False, check=False)[0].name)]#self.check_existence(extin, extout)
            for i, fitsid in enumerate(self.fitsid):
                filen = self.anadir/self.extpath(extin)[i].name #extin_noexist[i]
                print(filen)
                hdu = pyf.open(filen)
                im = hdu[0].data
                #header = hdu.header
                hdu.close()
                rsd,wav,mask,pixcoord,rotim,iys_plot,iye_plot = im_to_rsd(im,hotpix_mask=hotpix_mask,wavcal_path=wavcal_path)
        else:
            median_image=self.immedian()
            imall = []
            rsdall = []
            for i, fitsid in enumerate(self.fitsid):
                filen = self.anadir/self.extpath(extin)[i].name
                hdu = pyf.open(filen)
                im = hdu[0].data
                imall.append(im)
                hdu.close()
            imall = np.array(imall)
            median_image = np.nanmedian(imall, axis=0)
            rsd,wav,mask,pixcoord,rotim,iys_plot,iye_plot = im_to_rsd(median_image,hotpix_mask=hotpix_mask,wavcal_path=wavcal_path)

        self.fitsdir = self.anadir
        os.chdir(currentdir)
        return rsd,wav,mask,pixcoord,rotim,iys_plot,iye_plot#,xmin,xmax

    def apnormalize(self,rsd=None,hotpix_mask=None,ignore_orders=None):
        """normalize 2D apertures by 1D functions

        Returns:
            dictionary of pandas DataFrames of extracted and normalized spectra in each pixel
        """
        from pyird.image.oned_extract import flatten
        from pyird.image.trace_function import trace_legendre
        from pyird.spec.continuum import continuum_rsd
        from pyird.image.hotpix import apply_hotpixel_mask

        if rsd is None:
            flatfile = self.anadir/('%s_%s_%s.fits'%(self.streamid,self.band,self.trace.mmf))
            hdu=pyf.open(flatfile)[0]
            rsd=hdu.data
            header = hdu.header

        if not hotpix_mask is None:
            save_path = self.anadir/('hotpix_%s_%s.fits'%(self.band,self.trace.mmf))
            rsd = apply_hotpixel_mask(hotpix_mask, rsd, self.trace.y0, self.trace.xmin, self.trace.xmax, self.trace.coeff, save_path=save_path)

        df_continuum = continuum_rsd(rsd,ignore_orders=ignore_orders)

        flat_median = self.immedian('_cp')
        if not hotpix_mask is None:
            flat_median=flat_median*(~hotpix_mask)
        df_onepix = flatten(flat_median, trace_legendre, self.trace.y0, self.trace.xmin, self.trace.xmax, self.trace.coeff, inst=self.inst,onepix=True)#,width=[3,5])
        apertures = [int(i.split('ec')[-1]) for i in df_onepix.keys()]
        df_flatn = {}
        for i in apertures:
            df_flatn_tmp = df_onepix['ec%d'%(i)]/(df_continuum/len(apertures))
            df_flatn['ec%d'%(i)] = df_flatn_tmp
        return df_flatn

    def apext_flatfield(self, df_flatn, extout='_fln', extin=None, hotpix_mask=None, width=None):
        """aperture extraction and flat fielding (c.f., hdsis_ecf.cl)

        Args:
            df_flatn:
            extout: extension of output files
            extin: extension of input files
            hotpix_mask: hotpix masked spectrum ('extout' to be automatically '_flnhp')
            width: list of aperture widths ([width_start,width_end])
        """
        from pyird.image.hotpix import apply_hotpixel_mask
        def sum_weighted_apertures(im,df_flatn):
            from pyird.image.oned_extract import flatten
            from pyird.image.trace_function import trace_legendre
            df_onepix = flatten(im, trace_legendre, self.trace.y0, self.trace.xmin, self.trace.xmax, self.trace.coeff, inst=self.inst,onepix=True, width=self.trace.width)
            apertures = [int(i.split('ec')[-1]) for i in df_onepix.keys()]
            df_ecf = {}
            for i in apertures:
                #flatn_mean = np.nanmean(df_flatn['ec%d'%(i)].loc[2:len(df_flatn['ec%d'%(i)])-1].values.astype(float)) #cf) hdsis_ecf
                flatn_mean = np.nanmedian(df_flatn['ec%d'%(i)].loc[2:len(df_flatn['ec%d'%(i)])-1].values.astype(float))
                print("pixel = %d, Mean = %.5f"%(i,flatn_mean))
                df_ecf['ec%d'%(i)] = (df_onepix['ec%d'%(i)]/df_flatn['ec%d'%(i)])*flatn_mean #cf) hdsis_ecf
            for i,key in enumerate(df_ecf.keys()):
                if i==0:
                    df_sum_wap = df_ecf[key]
                else:
                    df_sum_wap += df_ecf[key]
            return df_sum_wap

        def rsd_order_medfilt(rsd,kernel_size=9):
            from scipy.signal import medfilt
            rsd_filtered = []
            for i in range(len(rsd[0])):
                filtered_data = medfilt(rsd[:,i],kernel_size=kernel_size)
                rsd_filtered.append(filtered_data)
            rsd_filtered = np.array(rsd_filtered).T
            return rsd_filtered

        if extin is None:
            extin = self.extension

        if hotpix_mask is None:
            extout = extout + '_' + self.trace.mmf
        else:
            extout = extout + 'hp_' + self.trace.mmf
        if not self.imcomb:
            extin_noexist, extout_noexist = self.check_existence(extin, extout)
            for i, fitsid in enumerate(tqdm.tqdm(extin_noexist)):
                filen = self.anadir/extin_noexist[i]
                hdu = pyf.open(filen)[0]
                im = hdu.data
                header = hdu.header
                df_sum_wap = sum_weighted_apertures(im,df_flatn)
                rsd = df_sum_wap.values.T.astype(float)
                if not hotpix_mask is None:
                    save_path = self.anadir/('hotpix_%s_%s.fits'%(self.band,self.trace.mmf))
                    rsd = apply_hotpixel_mask(hotpix_mask, rsd, self.trace.y0, self.trace.xmin, self.trace.xmax, self.trace.coeff, save_path=save_path)
                hdux = pyf.PrimaryHDU(rsd, header)
                hdulist = pyf.HDUList([hdux])
                hdulist.writeto(self.anadir/extout_noexist[i], overwrite=True)
        else:
            save_path = self.anadir/('n%s_%s_%s.fits'%(self.streamid,self.band,self.trace.mmf))
            if not os.path.exists(save_path):
                median_image=self.immedian('_cp')
                if not hotpix_mask is None:
                    median_image=median_image*(~hotpix_mask)
                df_sum_wap = sum_weighted_apertures(median_image,df_flatn)
                rsd = df_sum_wap.values.T.astype(float)
                rsd = rsd_order_medfilt(rsd)
                hdux = pyf.PrimaryHDU(data=rsd,header=None)
                hdulist = pyf.HDUList([hdux])
                hdulist.writeto(save_path, overwrite=True)

    def immedian(self, extension=None):
        """take image median.

        Return:
           median image
        """
        e = self.extension
        if not extension is None:
            self.extension = extension
        print('median combine: ',self.extension)
        imall = []
        for path in tqdm.tqdm(self.path(string=True, check=True)):
            imall.append(pyf.open(path)[0].data)
        imall = np.array(imall)
        # np.nansum(corrected_im_all,axis=0)#
        median_image = np.nanmedian(imall, axis=0)
        self.extension = e
        return median_image

    def calibrate_wavelength(self, trace_file_path=None, maxiter=30, stdlim=0.001, npix=2048, width=None):
        """wavelength calibration usgin Th-Ar.

        Args:
           trace_file_path: path to the trace file
           maxiter: maximum number of iterations
           stdlim: When the std of fitting residuals reaches this value, the iteration is terminated.
           npix: number of pixels
           width: list of aperture widths ([width_start,width_end])

        """
        from pyird.spec.wavcal import wavcal_thar
        from pyird.io.iraf_trace import read_trace_file
        from pyird.image.oned_extract import flatten
        from pyird.image.trace_function import trace_legendre
        from pyird.spec.rsdmat import multiorder_to_rsd

        def make_weight(): #REVIEW: there may be other appropreate weights
            path = (pkg_resources.resource_filename('pyird', 'data/IP_fwhms_h.dat'))
            fwhms = np.loadtxt(path)
            calc_wtmp = lambda fwhm: 1/(fwhm)
            w = []
            for order in range(len(fwhms)):
                w_ord = []
                for part in range(len(fwhms[0])):
                    if part != len(fwhms[0])-1:
                        w_tmp = [calc_wtmp(fwhms[order][part])]*108
                    else:
                        w_tmp = [calc_wtmp(fwhms[order][part])]*104
                    w_ord.extend(w_tmp)
                w.append(w_ord)
            w = np.array(w).T
            return w

        currentdir = os.getcwd()
        os.chdir(str(self.anadir))

        median_image = self.immedian()

        if trace_file_path is None:
            y0, xmin, xmax, coeff = self.trace.y0, self.trace.xmin, self.trace.xmax, self.trace.coeff
            mmf = self.trace.mmf
            width = self.trace.width
        else:
            y0, interp_function, xmin, xmax, coeff = read_trace_file(trace_file_path)

        master_path = 'thar_%s_%s.fits'%(self.band,mmf)
        if not os.path.exists(master_path):
            filen = self.path()[0] #header of the first file
            hdu = pyf.open(filen)[0]
            im = hdu.data
            header = hdu.header
            nord = len(y0)
            rawspec, pixcoord, _, _, _, _ = flatten(
                median_image, trace_legendre, y0, xmin, xmax, coeff, inst=self.inst, width=width)
            rsd = multiorder_to_rsd(rawspec, pixcoord)
            print(np.shape(rsd.T)[0])
            ## set weights
            if self.band=='h':
                #w = make_weight()
                w = np.ones(rsd.shape) #XXX!!
            else: #TODO: for y band
                w = np.ones(rsd.shape)
            wavsol, data = wavcal_thar(rsd.T, w, maxiter=maxiter, stdlim=stdlim)
            #np.save('thar_%s_%s_final.npy'%(self.band,mmf),data)
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
        inst = self.inst

        flatmedian=self.immedian()
        if nap==42 or nap==21:
            flatmedian = flatmedian[::-1,::-1]
        y0, xmin, xmax, coeff = aptrace(flatmedian,cutrow,nap)

        return TraceAperture(trace_legendre, y0, xmin, xmax, coeff, inst)

    def dispcor(self, extin='_fl', prefix='w', master_path=None, blaze=True):
        """dispersion correct and resample spectra

        Args:
            extin: extension of input files
            prefix: prefix for output files
            master: master file for the wavelength calibrated ThAr file
            master_path: path to the directory containing the master ThAr file

        """
        from pyird.plot.showspec import show_wavcal_spectrum
        def mkwspec(spec_m2,reference,save_path):
            wspec = pd.DataFrame([],columns=['wav','order','flux'])
            for i in range(len(reference[0])):
                wav = reference[:,i]
                order = np.ones(len(wav))
                order[:] = i+1
                data_order = [wav,order,spec_m2[:,i]]
                df_order = pd.DataFrame(data_order,index=['wav','order','flux']).T
                wspec = pd.concat([wspec,df_order])
            wspec = wspec.fillna(0)
            wspec.to_csv(save_path,**self.tocsvargs)
            return wspec

        if master_path==None:
            master_path = self.anadir.joinpath('..','thar').resolve()/('thar_%s_%s.fits'%(self.band,self.trace.mmf))
        elif not str(master_path).endswith('.fits'):
            master_path = master_path/('thar_%s_%s.fits'%(self.band,self.trace.mmf))

        extin = extin + '_' + self.trace.mmf

        if not self.imcomb:
            inputs = self.extpath(extin,string=False, check=True)
            for j,input in enumerate(inputs):
                hdu = pyf.open(input)[0]
                spec_m12 = hdu.data
                spec_m2 = spec_m12#[:,::2] # choose mmf2 (star fiber)

                hdu = pyf.open(master_path)[0]
                reference = hdu.data

                id = self.fitsid[j]
                save_path = self.anadir/('%s%d_%s.dat'%(prefix,id,self.trace.mmf)) ##e.g. w41511_m2.dat
                wspec = mkwspec(spec_m2,reference,save_path)
                if self.info:
                    print('dispcor: output spectrum= %s%d_%s.dat'%(prefix,id,self.trace.mmf))
                #plot
                show_wavcal_spectrum(wspec,alpha=0.5)
        else:
            def save_wspec(hdu,save_path):
                spec_m2 = hdu.data
                hdu_ref = pyf.open(master_path)[0]
                reference = hdu_ref.data
                wspec = mkwspec(spec_m2,reference,save_path)
                if self.info:
                    outpath = str(save_path).split('/')[-1]
                    print('dispcor: output spectrum= ', outpath)
                #plot
                show_wavcal_spectrum(wspec,alpha=0.5)

            if blaze:
                hdu = pyf.open(self.anadir/('n%s_%s_%s.fits'%(self.streamid,self.band,self.trace.mmf)))[0]
                save_path = self.anadir/('%sblaze_%s_%s.dat'%(prefix,self.band,self.trace.mmf))
                save_wspec(hdu,save_path)
            else:## need for fringe removal?
                hdu = pyf.open(self.anadir/('%s_%s_%s.fits'%(self.streamid,self.band,self.trace.mmf)))[0]
                save_path = self.anadir/('%s%s_%s_%s.dat'%(prefix,self.streamid,self.band,self.trace.mmf))
                save_wspec(hdu,save_path)


    def normalize1D(self,flatid='blaze',master_path=None,skipLFC=False):
        """combine orders and normalize spectrum

        Args:
            flatid: streamid for flat data
            master_path: path to the directory containing the calibrated flat file

        """
        from pyird.spec.continuum import comb_norm
        from pyird.plot.showspec import show_wavcal_spectrum

        if master_path==None:
            flatfile = self.anadir.joinpath('..','flat').resolve()/('w%s_%s_%s.dat'%(flatid,self.band,self.trace.mmf))
        elif not str(master_path).endswith('.dat'):
            flatfile = master_path/('w%s_%s_%s.dat'%(flatid,self.band,self.trace.mmf))

        if self.imcomb:
            fits_range = self.fitsid[:1]
        else:
            fits_range = self.fitsid

        for i,id in enumerate(fits_range):
            if self.imcomb:
                wfile = self.anadir/('w%s_%s_%s.dat'%(self.streamid,self.band,self.trace.mmf))
                nwsave_path = self.anadir/('nw%s_%s_%s.dat'%(self.streamid,self.band,self.trace.mmf))
                ncwsave_path = self.anadir/('ncw%s_%s_%s.dat'%(self.streamid,self.band,self.trace.mmf))
                if not skipLFC:
                    LFC_path = self.anadir/('w%s_y_m1.dat'%(self.streamid,))
                else:
                    LFC_path = None
            else:
                wfile = self.anadir/('w%d_%s.dat'%(id,self.trace.mmf))
                nwsave_path = self.anadir/('nw%d_%s.dat'%(id,self.trace.mmf))
                ncwsave_path = self.anadir/('ncw%d_%s.dat'%(id,self.trace.mmf))
                if not skipLFC:
                    if self.band == 'h':
                        LFC_path = self.anadir/('w%d_%s.dat'%(id-1,'m1'))
                    elif self.band == 'y':
                        LFC_path = self.anadir/('w%d_%s.dat'%(id,'m1'))
                else:
                    LFC_path = None
            df_continuum, df_interp = comb_norm(wfile,flatfile,LFC_path,blaze=(flatid=='blaze'))
            df_continuum_save = df_continuum[['wav','order','nflux','sn_ratio','uncertainty']]
            df_interp_save = df_interp[['wav','nflux','sn_ratio','uncertainty']]
            df_continuum_save.to_csv(nwsave_path,**self.tocsvargs)
            df_interp_save.to_csv(ncwsave_path,**self.tocsvargs)
            if self.info and ~self.imcomb:
                print('normalize1D: output normalized 1D spectrum= nw%d_%s.dat'%(id,self.trace.mmf))
            #plot
            show_wavcal_spectrum(df_continuum_save,alpha=0.5)
            show_wavcal_spectrum(df_interp_save,alpha=0.5)

class Stream1D(DatSet):
    def __init__(self, streamid, rawdir, anadir, fitsid=None, prefix='', extension='', inst='IRD'):
        """initialization
        Args:
           streamid: ID for stream
           rawdir: directory where the raw data are
           anadir: directory in which the processed file will put
           fitsid: fitsid

        """
        super(Stream1D, self).__init__(rawdir, prefix=prefix, extension=extension)
        self.streamid = streamid
        self.rawdir = rawdir
        self.anadir = anadir
        self.unlock = False
        self.info = False
        self.inst = inst
        self.readargs = {'header':None,'delim_whitespace':True,'names':['wav','order','flux','sn_ratio','uncertainty']}
        self.tocsvargs = {'header':False,'index':False,'sep':' '}
        if fitsid is not None:
            print('fitsid:', fitsid)
            self.fitsid = fitsid
            if (fitsid[0]%2==0):
                self.band = 'y'
            else:
                self.band = 'h'
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
        self.band = 'h'

    def fitsid_decrement(self):
        """Decrease fits id +1."""
        for i in range(0, len(self.fitsid)):
            self.fitsid[i] = self.fitsid[i]-1
        self.rawpath = self.path(string=False, check=True)

    ############################################################################################
    def extpath(self, extension, string=False, check=True):
        """path to file with extension.

        Args:
           extension: extension

        Returns:
           path array of fits files w/ extension
        """
        f = self.dir
        e = self.extension
        self.dir = self.anadir
        self.extension = extension
        path_ = self.path(string, check)
        self.dir = f
        self.extension = e
        return path_

    def specmedian(self,method='mean'):
        """take spectrum mean or median.

        Return:
           dataframe including median spectrum
        """
        specall,suffixes = [],[]
        for i,path in enumerate(tqdm.tqdm(self.path(string=True, check=True))):
            data = pd.read_csv(path,**self.readargs)
            suffix = 'flux_%d'%(i)
            suffixes.append(suffix)
            data = data.rename(columns={'flux':'flux_%d'%(i),'sn_ratio':'sn_ratio_%d'%(i),'uncertainty':'uncertainty_%d'%(i)})
            if i==0:
                df_merge = data
            else:
                df_merge = pd.merge(df_merge,data,on=['wav','order'])
        if method=='median':
            df_merge['flux_median'] = df_merge[suffixes].median(axis=1)
        elif method=='mean':
            df_merge['flux_median'] = df_merge[suffixes].mean(axis=1)
            df_merge['flux_err'] = np.sqrt(np.nansum(df_merge.filter(like='uncertainty') **2,axis=1))/len(df_merge.filter(like='uncertainty').columns)
        return df_merge

    def remove_fringe(self):
        """removing periodic noise (fringe) in REACH spectrum
        """
        from pyird.utils.remove_fringe import remove_fringe_order
        flatpath = self.flatpath(string=True,check=True)
        df_flat = pd.read_csv(flatpath,**self.readargs)

        df_targetmed = self.specmedian()

        orders = df_targetmed['order'].unique()
        rmfringe_all = []
        err_all = []
        for order in tqdm.tqdm(orders):
            flux_rmfringe, flux_err_rmfringe = remove_fringe_order(df_flat,df_targetmed,order,mask=True)
            rmfringe_all.extend(flux_rmfringe.values)
            err_all.extend(flux_err_rmfringe.values)
        df_targetmed['flux_rmfringe'] = rmfringe_all
        df_targetmed['flux_err_rmfringe'] = err_all
        save_path = self.anadir/('rfnw%s_%s_%s%s.dat'%(self.streamid,self.date,self.band,self.extension))
        df_targetmed_save = df_targetmed[['wav','order','flux_rmfringe','flux_err_rmfringe']]
        df_targetmed_save.to_csv(save_path,**self.tocsvargs)
        if self.info:
            print('removing fringe: output = rfnw%s_%s_%s%s.dat'%(self.streamid,self.date,self.band,self.extension))
