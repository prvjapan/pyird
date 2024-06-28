"""File stream for IRD analysis."""

from pyird.utils.fitsset import FitsSet
from pyird.utils.datset import DatSet
from pyird.image.trace_function import trace_legendre
from pyird.image.aptrace import aptrace
from pyird.utils.aperture import TraceAperture
import astropy.io.fits as pyf
import numpy as np
import pandas as pd
import tqdm
import os

import pkg_resources

__all__ = ["Stream1D", "Stream2D"]


class Stream2D(FitsSet):
    def __init__(
        self,
        streamid,
        rawdir,
        anadir,
        fitsid=None,
        rawtag="IRDA000",
        inst="IRD",
        rotate=False,
        inverse=False,
        detector_artifact=False,
    ):
        """initialization
        Args:
            streamid: ID for stream
            rawdir: directory where the raw data are
            anadir: directory in which the processed file will put
            fitsid: fitsid list, such as [10301,10303]
            rawtag:
            rotate (boolen): If True, the image is rotated in 90 deg (for old detector). See #80 in GitHub
            inverse (boolen): If True, the image is inversed along y axis. See #80 in GitHub
            detector_artifact (boolen): If True, fill the gaps seen in the old detector. See #80 in GitHub
        """
        super(Stream2D, self).__init__(rawtag, rawdir, extension="")
        self.streamid = streamid
        self.rawdir = rawdir
        self.anadir = anadir
        self.unlock = False
        self.info = False
        self.imcomb = False
        self.inst = inst
        self.rotate = rotate
        self.inverse = inverse
        self.detector_artifact = detector_artifact

        self.tocsvargs = {"header": False, "index": False, "sep": " "}
        if fitsid is not None:
            print("fitsid:", fitsid)
            self.fitsid = fitsid
            if (rawtag == "IRDA000" and fitsid[0] % 2 == 0) or (rawtag == "IRDBD000"):
                self.band = "y"
            else:
                self.band = "h"
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
            self.fitsid[i] = self.fitsid[i] + 1
        self.rawpath = self.path(string=False, check=True)
        self.band = "h"

    def fitsid_decrement(self):
        """Decrease fits id +1."""
        for i in range(0, len(self.fitsid)):
            self.fitsid[i] = self.fitsid[i] - 1
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

        if extension == "":
            print("extclean cannot clean w/o extension fits.")
            return

        for i in range(len(self.path())):
            if self.path(check=False)[i].exists():
                os.remove(self.extpath(extension, check=False, string=True)[i])
                print("rm old " + self.extpath(extension, check=False, string=True)[i])

    def remove_bias(self, rot=None):
        self.print_if_info_is_true("remove_bias: files="+str(self.rawpath))
        if rot == "r":
            print("180 degree rotation applied.")
        self.fitsdir = self.anadir
        self.extension = "_rb"

    def clean_pattern(
        self, trace_mask=None, hotpix_mask=None, extout="_cp", extin=None
    ):
        """

        Args:
            trace_mask: trace mask (from iraf, use image.mask.trace_from_iraf_trace_file)
            hotpix_mask: hot pixel mask
            extout: output extension
            extin: input extension

        """
        from pyird.image.pattern_model import median_XY_profile

        currentdir = os.getcwd()
        os.chdir(str(self.anadir))
        self.print_if_info_is_true("clean_pattern: output extension="+extout)

        if trace_mask is None:
            trace_mask = self.trace.mask()

        if extin is None:
            extin = self.extension

        extin_noexist, extout_noexist = self.check_existence(extin, extout)
        for i, fitsid in enumerate(tqdm.tqdm(extin_noexist)):
            filen = self.rawdir / extin_noexist[i]
            im, header = self.load_image(filen)
            calim = np.copy(im)  # image for calibration
            if self.band == "y":
                trace_mask = trace_mask.reshape(1, int(trace_mask.size))
                trace_mask = trace_mask[0][::-1].reshape(
                    im.shape[0], im.shape[1]
                )  # rotate mask matrix 180 degree
            calim[trace_mask] = np.nan
            if hotpix_mask is not None:
                calim[hotpix_mask] = np.nan

            model_im = median_XY_profile(calim)
            corrected_im = im - model_im

            self.write_image(extout_noexist[i], header, corrected_im)

        self.fitsdir = self.anadir
        self.extension = extout
        os.chdir(currentdir)

    def print_if_info_is_true(self, msg):
        if self.info:
            print(msg)

    def load_image(self, filen):
        """
        loads a fits image, if self.rotate=True, the image is rotated in 90 deg, self.inverse=True, the image is inversed along the 0-th-axis

        Args:
            filen (str): filename

        Returns:
            image, header
        """

        img, header = self.load_fits_data_header(filen)
        if self.rotate:
            img = np.rot90(img)
            self.print_if_info_is_true("rotates the image in 90 deg when loading")
            
        if self.inverse:
            img = img[::-1, :]
            self.print_if_info_is_true("inverse the image along the 0th-axis when loading")
            
        return img, header

    def load_rsd(self, filen):
        """
        loads Raw Spectral Detector matrix (RSD matrix)

        Args:
            filen (str): filename

        Returns:
            rsd, header
        """
        return self.load_fits_data_header(filen)

    def load_fits_data_header(self, filen):
        """
        loads data and header from fits

        Args:
            filen (str): filename

        Returns:
            data, header
        """
        hdu = pyf.open(filen)[0]
        data = hdu.data
        header = hdu.header
        hdu._close()
        return data, header

    def write_image(self, extout_noexist, header, img):
        """
        write an image to fits format, self.inverse=True, the image is inversed along the 0-th-axis, if self.rotate=True, the image is rotated in - 90 deg,

        Args:
            extout_noexist (_type_): _description_
            header (_type_): _description_
            img (_type_): _description_
        """
        if self.inverse:
            img = img[::-1, :]
        if self.rotate:
            img = np.rot90(img, k=-1)
        return self.write_fits_data_header(extout_noexist, header, img)

    def write_rsd(self, extout_noexist, header, rsd):
        """
        write Raw Spectral Detector matrix (RSD matrix) to fits format

        Args:
            extout_noexist (_type_): _description_
            header (_type_): _description_
            rsd (_type_): _description_
        """
        self.write_fits_data_header(extout_noexist, header, rsd)

    def write_fits_data_header(self, extout_noexist, header, rsd):
        hdu = pyf.PrimaryHDU(rsd, header)
        hdulist = pyf.HDUList([hdu])
        hdulist.writeto(self.anadir / extout_noexist, overwrite=True)
        hdu._close()

    def extract_trace_info(self,trace_path=None):
        """
        extract parameters for trace

        Args:
            trace_path: path to trace file

        Returns:
            data and parameters for aperture trace
        """
        if trace_path is None:
            y0, xmin, xmax, coeff = (
                self.trace.y0,
                self.trace.xmin,
                self.trace.xmax,
                self.trace.coeff,
            )
        else:
            y0, interp_function, xmin, xmax, coeff = read_trace_file(trace_path)
        mmf = self.trace.mmf
        width = self.trace.width
        return y0, xmin, xmax, coeff, mmf, width

    def flatten(
        self, trace_path=None, extout="_fl", extin=None, hotpix_mask=None, width=None, check=False, master_path=None
    ):
        """
        Args:
            trace_path: trace file to be used in flatten
            extout: output extension
            extin: input extension
            hotpix_mask: hotpix masked spectrum ('extout' to be automatically '_hp')
            width: list of aperture widths ([width_start,width_end])
            check: if True, return the extracted spectrum
            master_path: if this set the path to wavelength reference file with check=True, return wavelength allocated spectrum

        """
        from pyird.image.oned_extract import flatten
        from pyird.image.trace_function import trace_legendre
        from pyird.io.iraf_trace import read_trace_file
        from pyird.spec.rsdmat import multiorder_to_rsd
        from pyird.image.hotpix import apply_hotpixel_mask

        currentdir = os.getcwd()
        os.chdir(str(self.anadir))

        y0, xmin, xmax, coeff, mmf, width = self.extract_trace_info(trace_path)

        if extin is None:
            extin = self.extension

        extout = extout + "_" + mmf
        if not hotpix_mask is None:
            extout = "_hp_" + mmf

        if self.imcomb:
            out_path = self.anadir / (
                "%s_%s_%s.fits" % (self.streamid, self.band, mmf)
            )   
            extout_noexist = [out_path] * (not os.path.exists(out_path))
        else:
            extin_noexist, extout_noexist = self.check_existence(extin, extout)

        for i in tqdm.tqdm(range(len(extout_noexist))):
            if self.imcomb:
                im = self.immedian()  # median combine
                filen = self.path()[0]  # header of the first file
                _, header = self.load_image(filen)
            else:
                filen = self.anadir / extin_noexist[i]
                im, header = self.load_image(filen)
            rawspec, pixcoord, rotim, tl, iys_plot, iye_plot = flatten(
                im,
                trace_legendre,
                y0,
                xmin,
                xmax,
                coeff,
                inst=self.inst,
                width=width,
            )
            rsd = multiorder_to_rsd(rawspec, pixcoord)
            if not hotpix_mask is None:
                save_path = self.anadir / (
                    "hotpix_%s_%s.fits" % (self.band, self.trace.mmf)
                )
                rsd = apply_hotpixel_mask(
                    hotpix_mask, rsd, y0, xmin, xmax, coeff, save_path=save_path
                )
            if not check:
                self.write_rsd(extout_noexist[i], header, rsd)
            else:
                mask_shape = (2048, 2048)
                trace_mask = np.zeros(mask_shape)
                for i in range(len(y0)):
                    trace_mask[i, xmin[i] : xmax[i] + 1] = tl[i]
                if not master_path is None:
                    wav, _ = self.load_image(master_path)
                else:
                    wav = []
                return rsd,wav,trace_mask,pixcoord,rotim,iys_plot,iye_plot
        
        self.print_if_info_is_true("flatten (+ hotpix mask): output extension="+extout)
        self.fitsdir = self.anadir
        self.extension = extout
        os.chdir(currentdir)

    def apnormalize(self, rsd=None, hotpix_mask=None, ignore_orders=None):
        """normalize 2D apertures by 1D functions

        Returns:
            dictionary of pandas DataFrames of extracted and normalized spectra in each pixel
        """
        from pyird.image.oned_extract import flatten
        from pyird.image.trace_function import trace_legendre
        from pyird.spec.continuum import continuum_rsd
        from pyird.image.hotpix import apply_hotpixel_mask

        if rsd is None:
            flatfile = self.anadir / (
                "%s_%s_%s.fits" % (self.streamid, self.band, self.trace.mmf)
            )
            rsd, header = self.load_rsd(flatfile)

        if not hotpix_mask is None:
            save_path = self.anadir / (
                "hotpix_%s_%s.fits" % (self.band, self.trace.mmf)
            )
            rsd = apply_hotpixel_mask(
                hotpix_mask,
                rsd,
                self.trace.y0,
                self.trace.xmin,
                self.trace.xmax,
                self.trace.coeff,
                save_path=save_path,
            )

        df_continuum = continuum_rsd(rsd, ignore_orders=ignore_orders)

        flat_median = self.immedian("_cp")
        if not hotpix_mask is None:
            flat_median[hotpix_mask] = np.nan
        df_onepix = flatten(
            flat_median,
            trace_legendre,
            self.trace.y0,
            self.trace.xmin,
            self.trace.xmax,
            self.trace.coeff,
            inst=self.inst,
            onepix=True,
        )  # ,width=[3,5])
        apertures = [int(i.split("ec")[-1]) for i in df_onepix.keys()]
        df_flatn = {}
        for i in apertures:
            df_flatn_tmp = df_onepix["ec%d" % (i)] / (df_continuum / len(apertures))
            df_flatn["ec%d" % (i)] = df_flatn_tmp
        return df_flatn

    def apext_flatfield(
        self, df_flatn, extout="_fln", extin=None, hotpix_mask=None, width=None
    ):
        """aperture extraction and flat fielding (c.f., hdsis_ecf.cl)

        Args:
            df_flatn:
            extout: extension of output files
            extin: extension of input files
            hotpix_mask: hotpix masked spectrum ('extout' to be automatically '_flnhp')
            width: list of aperture widths ([width_start,width_end])
        """
        from pyird.image.hotpix import apply_hotpixel_mask
        from pyird.spec.rsdmat import rsd_order_medfilt
        from pyird.image.oned_extract import sum_weighted_apertures

        y0, xmin, xmax, coeff, mmf, width = self.extract_trace_info()

        if extin is None:
            extin = self.extension

        if hotpix_mask is None:
            extout = extout + "_" + self.trace.mmf
        else:
            extout = extout + "hp_" + self.trace.mmf

        if self.imcomb:
            out_path = self.anadir / (
                "n%s_%s_%s.fits" % (self.streamid, self.band, self.trace.mmf)
            )
            extout_noexist = [out_path] * (not os.path.exists(out_path))
        else:
            extin_noexist, extout_noexist = self.check_existence(extin, extout)
        
        for i in tqdm.tqdm(range(len(extout_noexist))):
            if self.imcomb:
                im = self.immedian("_cp")
                header = None
            else:
                filen = self.anadir / extin_noexist[i]
                im, header = self.load_image(filen)
            df_sum_wap = sum_weighted_apertures(im, df_flatn, y0, xmin, xmax, coeff, width, self.inst)
            rsd = df_sum_wap.values.T.astype(float)
            if not hotpix_mask is None:
                save_path = self.anadir / (
                    "hotpix_%s_%s.fits" % (self.band, self.trace.mmf)
                )
                rsd = apply_hotpixel_mask(
                    hotpix_mask,
                    rsd,
                    self.trace.y0,
                    self.trace.xmin,
                    self.trace.xmax,
                    self.trace.coeff,
                    save_path=save_path,
                )
            if self.imcomb:
                rsd = rsd_order_medfilt(rsd)
            self.write_rsd(extout_noexist[i], header, rsd)

    def immedian(self, extension=None):
        """take image median.

        Return:
           median image
        """
        e = self.extension
        if not extension is None:
            self.extension = extension
        print("median combine: ", self.extension)
        imall = []
        for path in tqdm.tqdm(self.path(string=True, check=True)):
            img, header = self.load_image(path)
            imall.append(img)
        imall = np.array(imall)
        # np.nansum(corrected_im_all,axis=0)#
        median_image = np.nanmedian(imall, axis=0)
        self.extension = e
        return median_image

    def calibrate_wavelength(
        self,
        trace_path=None,
        maxiter=30,
        std_threshold=0.001,
        npix=2048,
        width=None,
        force_rotate=False,
    ):
        """wavelength calibration usgin Th-Ar.

        Args:
            trace_path: path to the trace file
            maxiter: maximum number of iterations
            std_threshold: When the std of fitting residuals reaches this value, the iteration is terminated.
            npix: number of pixels
            width: list of aperture widths ([width_start,width_end])
            force_rotate: forces rotating the detector, when the number of the order is not standard values (i.e. 21 or 51)


        """
        from pyird.spec.wavcal import wavcal_thar
        from pyird.io.iraf_trace import read_trace_file
        from pyird.image.oned_extract import flatten
        from pyird.image.trace_function import trace_legendre
        from pyird.spec.rsdmat import multiorder_to_rsd

        currentdir = os.getcwd()
        os.chdir(str(self.anadir))

        median_image = self.immedian()

        y0, xmin, xmax, coeff, mmf, width = self.extract_trace_info(trace_path)

        master_path = "thar_%s_%s.fits" % (self.band, mmf)
        if not os.path.exists(master_path):
            filen = self.path()[0]  # header of the first file
            _, header = self.load_image(filen)
            nord = len(y0)
            rawspec, pixcoord, _, _, _, _ = flatten(
                median_image,
                trace_legendre,
                y0,
                xmin,
                xmax,
                coeff,
                inst=self.inst,
                width=width,
                force_rotate=force_rotate,
            )
            rsd = multiorder_to_rsd(rawspec, pixcoord)
            w = np.ones(rsd.shape) ##weights for identifying thar
            wavsol, data = wavcal_thar(
                rsd.T,
                w,
                maxiter=maxiter,
                std_threshold=std_threshold,
            )
            # np.save('thar_%s_%s_final.npy'%(self.band,mmf),data)
            wavsol_2d = wavsol.reshape((npix, nord))

            hdux = pyf.PrimaryHDU(wavsol_2d, header)
            hdulist = pyf.HDUList([hdux])
            hdulist.writeto(master_path, overwrite=True)

        os.chdir(currentdir)

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
                self.print_if_info_is_true("Ignore "+str(ext[i].name)+" -> "+str(extfi.name))
                skip = skip + 1

        if skip > 1:
            print("Skipped " + str(skip) + " files because they already exists.")

        return ext_noexist, extf_noexist

    def aptrace(self, cutrow=1000, nap=42):
        """extract aperture of trace from a median image of current fitsset

        Args:
            cutrow: cutting criterion
            nap: number of apertures, nap = 42 ## 42 for H band, 102 for YJ band

        Returns:
            TraceAperture instance

        """
        inst = self.inst

        flatmedian = self.immedian()
        if nap == 42 or nap == 21:
            flatmedian = flatmedian[::-1, ::-1]

        if self.detector_artifact:
            for i in range(0, 16):
                flatmedian[63 + i * 128 : 63 + i * 128 + 1, :] = flatmedian[
                    63 + i * 128 - 1 : 63 + i * 128, :
                ]
                flatmedian[63 + i * 128 + 1 : 63 + i * 128 + 2, :] = flatmedian[
                    63 + i * 128 + 2 : 63 + i * 128 + 3, :
                ]

        y0, xmin, xmax, coeff = aptrace(flatmedian, cutrow, nap)

        return TraceAperture(trace_legendre, y0, xmin, xmax, coeff, inst)

    def write_df_spec_wav(self, spec_m2, reference, save_path):
        """
        write spectrum and wavelengths to pandas.DataFrame format

        Args:
            spec_m2:    spectrum not yet allocated wavelength
            reference:  wavelength reference
            save_path:  path to .csv file 
        """
        wspec = pd.DataFrame([], columns=["wav", "order", "flux"])
        for i in range(len(reference[0])):
            wav = reference[:, i]
            order = np.ones(len(wav))
            order[:] = i + 1
            data_order = [wav, order, spec_m2[:, i]]
            df_order = pd.DataFrame(data_order, index=["wav", "order", "flux"]).T
            wspec = pd.concat([wspec, df_order])
        wspec = wspec.fillna(0)
        wspec.to_csv(save_path, **self.tocsvargs)
        return wspec

    def dispcor(self, extin="_fl", prefix="w", master_path=None, blaze=True):
        """dispersion correct and resample spectra

        Args:
            extin: extension of input files
            prefix: prefix for output files
            master: master file for the wavelength calibrated ThAr file
            master_path: path to the directory containing the master ThAr file

        """
        from pyird.plot.showspec import show_wavcal_spectrum
        from pyird.utils.datset import DatSet

        if master_path == None:
            master_path = self.anadir.joinpath("..", "thar").resolve() / (
                "thar_%s_%s.fits" % (self.band, self.trace.mmf)
            )
        elif not str(master_path).endswith(".fits"):
            master_path = master_path / (
                "thar_%s_%s.fits" % (self.band, self.trace.mmf)
            )

        extin = extin + "_" + self.trace.mmf

        hdu = pyf.open(master_path)[0]
        reference = hdu.data

        datset = DatSet(dir=self.anadir)
        datset.prefix = prefix
        if self.imcomb:
            datset.extension = "_%s_%s" % (self.band, self.trace.mmf)
            if blaze:
                inputs = [self.anadir / ("n%s_%s_%s.fits" % (self.streamid, self.band, self.trace.mmf))]
                datset.fitsid = ["blaze"]
            else:   ## need for fringe removal?
                inputs = [self.anadir / ("%s_%s_%s.fits" % (self.streamid, self.band, self.trace.mmf))]
                datset.fitsid = self.streamid
        else:
            inputs = self.extpath(extin, string=False, check=True)
            datset.extension = "_%s" % (self.trace.mmf)
            datset.fitsid = self.fitsid
        save_path = datset.path(string=True,check=False)

        for i, input in enumerate(inputs):
            hdu = pyf.open(input)[0]
            spec_m2 = hdu.data

            wspec = self.write_df_spec_wav(spec_m2, reference, save_path[i])

            outpath = str(save_path[i]).split("/")[-1]
            self.print_if_info_is_true(
                "dispcor: output spectrum= %s"
                % (outpath)
            )
            # plot
            show_wavcal_spectrum(wspec, title="Extracted spectrum: %s"%(outpath),alpha=0.5)

    def normalize1D(self, flatid="blaze", master_path=None, skipLFC=False):
        """combine orders and normalize spectrum

        Args:
            flatid: streamid for flat data
            master_path: path to the directory containing the calibrated flat file

        """
        from pyird.spec.continuum import comb_norm
        from pyird.plot.showspec import show_wavcal_spectrum
        from pyird.utils.datset import DatSet

        if master_path == None:
            flatfile = self.anadir.joinpath("..", "flat").resolve() / (
                "w%s_%s_%s.dat" % (flatid, self.band, self.trace.mmf)
            )
        elif not str(master_path).endswith(".dat"):
            flatfile = master_path / (
                "w%s_%s_%s.dat" % (flatid, self.band, self.trace.mmf)
            )

        datset = DatSet(dir=self.anadir)
        if self.imcomb:
            datset.fitsid = [self.streamid]
            datset.extension = "_%s_%s" % (self.band, self.trace.mmf)
            if not skipLFC:
                LFC_path = [self.anadir / ("w%s_y_m1.dat" % (self.streamid,))]
        else:
            datset.fitsid = self.fitsid
            datset.extension = "_%s" % (self.trace.mmf)
            if not skipLFC:
                LFC_path = []
                for id in self.fitsid:
                    if self.band == "h":
                        LFC_path.append(self.anadir / ("w%d_%s.dat" % (id - 1, "m1")))
                    elif self.band == "y":
                        LFC_path.append(self.anadir / ("w%d_%s.dat" % (id, "m1")))
        if skipLFC:
            LFC_path = [None]
        datset.prefix = "w"
        wfile = datset.path(string=True,check=False)
        datset.prefix = "nw"
        nwsave_path = datset.path(string=True,check=False)
        datset.prefix = "ncw"
        ncwsave_path = datset.path(string=True,check=False)

        for i in range(len(wfile)):
            df_continuum, df_interp = comb_norm(
                wfile[i], flatfile, LFC_path[i], blaze=(flatid == "blaze")
            )
            df_continuum_save = df_continuum[
                ["wav", "order", "nflux", "sn_ratio", "uncertainty"]
            ]
            df_interp_save = df_interp[["wav", "nflux", "sn_ratio", "uncertainty"]]
            df_continuum_save.to_csv(nwsave_path[i], **self.tocsvargs)
            df_interp_save.to_csv(ncwsave_path[i], **self.tocsvargs)
            if self.info and ~self.imcomb:
                print(
                    "normalize1D: output normalized 1D spectrum= nw%d_%s.dat"
                    % (self.fitsid[i], self.trace.mmf)
                )
            # plot
            show_wavcal_spectrum(df_continuum_save,title="Normalized spectrum: nw%d_%s.dat"
                                % (self.fitsid[i], self.trace.mmf),alpha=0.5)
            show_wavcal_spectrum(df_interp_save, title="Normalized & Order combined spectrum: ncw%d_%s.dat"
                                % (self.fitsid[i], self.trace.mmf),alpha=0.5)


class Stream1D(DatSet):
    def __init__(
        self, streamid, rawdir, anadir, fitsid=None, prefix="", extension="", inst="IRD"
    ):
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
        self.readargs = {
            "header": None,
            "delim_whitespace": True,
            "names": ["wav", "order", "flux", "sn_ratio", "uncertainty"],
        }
        self.tocsvargs = {"header": False, "index": False, "sep": " "}
        if fitsid is not None:
            print("fitsid:", fitsid)
            self.fitsid = fitsid
            if fitsid[0] % 2 == 0:
                self.band = "y"
            else:
                self.band = "h"
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
            self.fitsid[i] = self.fitsid[i] + 1
        self.rawpath = self.path(string=False, check=True)
        self.band = "h"

    def fitsid_decrement(self):
        """Decrease fits id +1."""
        for i in range(0, len(self.fitsid)):
            self.fitsid[i] = self.fitsid[i] - 1
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

    def specmedian(self, method="mean"):
        """take spectrum mean or median.

        Return:
           dataframe including median spectrum
        """
        specall, suffixes = [], []
        for i, path in enumerate(tqdm.tqdm(self.path(string=True, check=True))):
            data = pd.read_csv(path, **self.readargs)
            suffix = "flux_%d" % (i)
            suffixes.append(suffix)
            data = data.rename(
                columns={
                    "flux": "flux_%d" % (i),
                    "sn_ratio": "sn_ratio_%d" % (i),
                    "uncertainty": "uncertainty_%d" % (i),
                }
            )
            if i == 0:
                df_merge = data
            else:
                df_merge = pd.merge(df_merge, data, on=["wav", "order"])
        if method == "median":
            df_merge["flux_median"] = df_merge[suffixes].median(axis=1)
        elif method == "mean":
            df_merge["flux_median"] = df_merge[suffixes].mean(axis=1)
            df_merge["flux_err"] = np.sqrt(
                np.nansum(df_merge.filter(like="uncertainty") ** 2, axis=1)
            ) / len(df_merge.filter(like="uncertainty").columns)
        return df_merge

    def remove_fringe(self):
        """removing periodic noise (fringe) in REACH spectrum"""
        from pyird.utils.remove_fringe import remove_fringe_order

        flatpath = self.flatpath(string=True, check=True)
        df_flat = pd.read_csv(flatpath, **self.readargs)

        df_targetmed = self.specmedian()

        orders = df_targetmed["order"].unique()
        rmfringe_all = []
        err_all = []
        for order in tqdm.tqdm(orders):
            flux_rmfringe, flux_err_rmfringe = remove_fringe_order(
                df_flat, df_targetmed, order, mask=True
            )
            rmfringe_all.extend(flux_rmfringe.values)
            err_all.extend(flux_err_rmfringe.values)
        df_targetmed["flux_rmfringe"] = rmfringe_all
        df_targetmed["flux_err_rmfringe"] = err_all
        save_path = self.anadir / (
            "rfnw%s_%s_%s%s.dat" % (self.streamid, self.date, self.band, self.extension)
        )
        df_targetmed_save = df_targetmed[
            ["wav", "order", "flux_rmfringe", "flux_err_rmfringe"]
        ]
        df_targetmed_save.to_csv(save_path, **self.tocsvargs)
        if self.info:
            print(
                "removing fringe: output = rfnw%s_%s_%s%s.dat"
                % (self.streamid, self.date, self.band, self.extension)
            )
