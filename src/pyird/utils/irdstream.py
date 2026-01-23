"""File stream for IRD analysis."""

from pyird.utils.fitsset import FitsSet
from pyird.utils.datset import DatSet
from pyird.image.trace_function import trace_legendre
from pyird.image.oned_extract import flatten
from pyird.spec.rsdmat import multiorder_to_rsd
from pyird.image.hotpix import apply_hotpixel_mask
from pyird.utils.class_handler import update_attr
from pyird.utils.fitsutils import load_fits_data_header, write_fits_data_header
from pyird.plot.showspec import show_wavcal_spectrum
import numpy as np
import pandas as pd
import tqdm
import os
import warnings

from pyird.utils.decorators import deprecate_kwargs, rename_kwargs
import functools
rename_12_20 = functools.partial(deprecate_kwargs, since="v1.2", remove_in="v2.0")

__all__ = ["StreamCommon", "Stream1D", "Stream2D"]


class StreamCommon:
    """Common Stream class for 1D and 2D.

    Attributes:
        streamid: ID for stream
        rawdir: directory where the raw data are
        anadir: directory in which the processed file will put
        fitsid: fitsid list, such as [10301,10303]
        rawtag: prefix of file name, such as "IRDA000", "IRDAD000", or "IRDBD000"
        rotate (boolen): If True, the image is rotated in 90 deg (for old detector). See #80 in GitHub
        inverse (boolen): If True, the image is inversed along y axis. See #80 in GitHub
        detector_artifact (boolen): If True, fill the gaps seen in the old detector. See #80 in GitHub
        band: band of the data, "y" or "h"
        imcomb: If True, the median combine is applied to the data
        tocvsargs: arguments for to_csv
        info: If True, print the information
        inst: instrument name, "IRD" or "REACH"
        fitsid_already_increment: If True, the fitsid is already incremented

    """

    def __init__(self, streamid, rawdir, anadir, inst):
        self.streamid = streamid
        self.rawdir = rawdir
        self.anadir = anadir
        self.unlock = False
        self.info = True
        self.inst = inst
        self.fitsid_already_increment = False

    @property
    def fitsid(self):
        return self._fitsid

    @fitsid.setter
    def fitsid(self, fitsid):
        self._fitsid = fitsid
        self.rawpath = self.path(string=False, check=True)

    def fitsid_increment(self):
        """Increase fits id +1. Only one time increment is allowed."""
        if self.fitsid_already_increment:
            raise ValueError("fitsid is already incremented.")
        else:
            for i in range(0, len(self.fitsid)):
                self.fitsid[i] = self.fitsid[i] + 1
            self.rawpath = self.path(string=False, check=True)
            self.band = "h"
            self.fitsid_already_increment = True
            print("fitsid is incremented.")

    def fitsid_decrement(self):
        """Decrease fits id -1."""
        if self.fitsid_already_increment:
            for i in range(0, len(self.fitsid)):
                self.fitsid[i] = self.fitsid[i] - 1
            self.rawpath = self.path(string=False, check=True)
            self.fitsid_already_increment = False
            print("fitsid is decremented.")
        else:
            raise ValueError("fitsid is not incremented yet.")

    def extpath(self, extension, string=False, check=True):
        """path to file with extension.

        Args:
            extension: extension

        Returns:
            path array of fits files w/ extension
        """
        if hasattr(self, "fitsdir"):
            f = self.fitsdir
            self.fitsdir = self.anadir

        elif hasattr(self, "dir"):
            f = self.dir
            self.dir = self.anadir

        e = self.extension
        self.extension = extension
        path_ = self.path(string, check)
        self.dir = f
        self.extension = e
        return path_

    def extclean(self, extension):
        """Clean i.e. remove fits files if exists.

        Args:
            extension: extension of which files to be removed
        """
        if extension == "":
            print("extclean cannot clean w/o extension fits.")
            return

        for i in range(len(self.path())):
            if self.path(check=False)[i].exists():
                os.remove(self.extpath(extension, check=False, string=True)[i])
                print("rmoved old " + self.extpath(extension, check=False, string=True)[i])

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
                self.print_if_info_is_true(
                    "Ignore " + str(ext[i].name) + " -> " + str(extfi.name)
                )
                skip = skip + 1

        if skip > 1:
            print("Skipped " + str(skip) + " files because they already exists.")

        return ext_noexist, extf_noexist

    def remove_bias(self, rot=None):
        """set extension when loading an image whose bias is removed by IRD_bias_sub.py"""
        self.print_if_info_is_true("Performing remove_bias: files=" + str(self.rawpath))
        if rot == "r":
            print("180 degree rotation applied.")
        self.fitsdir = self.anadir
        self.extension = "_rb"

    def print_if_info_is_true(self, msg):
        if self.info:
            print(msg)


class Stream2D(FitsSet, StreamCommon):
    """Class for processing 2D spectral images and reducing them to 1D spectra.
    """
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
        band=None,
        not_ignore_warning = True    
    ):
        """initialization
        Args:
            streamid: ID for stream
            rawdir: directory where the raw data are
            anadir: directory in which the processed file will put
            fitsid: fitsid list, such as [10301,10303]
            rawtag: prefix of file name, such as "IRDA000", "IRDAD000", or "IRDBD000"
            rotate (boolen): If True, the image is rotated in 90 deg (for old detector). See #80 in GitHub
            inverse (boolen): If True, the image is inversed along y axis. See #80 in GitHub
            detector_artifact (boolen): If True, fill the gaps seen in the old detector. See #80 in GitHub
            band: band of the data, "y" or "h", requires fitsid
            not_ignore_warning: If True and files for fitsid do not exist, raise a ValueError (default: True)
        """
        FitsSet.__init__(self, rawtag, rawdir, extension="")
        StreamCommon.__init__(self, streamid, rawdir, anadir, inst)
        self.not_ignore_warning = not_ignore_warning

        self.imcomb = False
        self.rotate = rotate
        self.inverse = inverse
        self.detector_artifact = detector_artifact

        self.tocsvargs = {"header": False, "index": False, "sep": " "}

        self.fitsid = fitsid
        self.band = None

        if inst in ["IRD", "REACH"]:
            self.init_band_IRD(rawtag, band)
        elif inst == "IRCS":
            # Not the band name — "h" is assumed as the detector configuration for IRCS data in PyIRD.
            # This naming convention may change in the future.
            self.band = "h"
        print(f"Processing {self.band} band")

        if self.fitsid is None:
            print("No fitsid yet.")
        else:
            print("Processing fitsid:", self.fitsid)

    def init_band_IRD(self, rawtag, band):
        """initialize band for IRD/REACH data
        
        Args:
            rawtag: prefix of file name, such as "IRDA000", "IRDAD000", or "IRDBD000
            band: band of the data, "y" or "h", requires fitsid
        """
        if rawtag == "IRDBD000":
            if band is not None and band != "y":
                raise ValueError("band must be 'y' for IRDBD000")
            self.band = "y"

        elif rawtag == "IRDAD000":
            if band is not None and band != "h":
                raise ValueError("band must be 'y' for IRDAD000")
            self.band = "h"

        elif rawtag == "IRDA000":
            if self.fitsid is None:
                if band is not None:
                    raise ValueError("band option requires fitsid")
            else:
                first = int(self.fitsid[0]) # fitsid of the first file
                if band is None:
                    self.band = "y" if (first % 2 == 0) else "h"
                else:
                    if band not in ("y", "h"):
                        raise AttributeError("band must be 'y' or 'h'")
                    self.band = band
                    expected_parity = 0 if self.band == "y" else 1
                    if (first % 2) != expected_parity:
                        if self.band == "y" and (first % 2 == 1):
                            raise ValueError(
                                "For rawtag='IRDA000', band='y' corresponds to fitsid with EVEN numbers.\n"
                                f"You set band={self.band} and fitsid starting with {first} (ODD number).\n"
                                "- To analyze h band: set band='h'.\n"
                                "- To analyze y band: specify fitsid ±1 accordingly.\n"
                                "- Does your filename start with IRDAD or IRDBD? Then set rawtag='IRDAD' or 'IRDBD'."
                            )
                        if self.band == "h" and (first % 2 == 0):
                            self.fitsid_increment()

    def detector_handling(self, img, mode="load"):
        """apply rotation and/or inverse to image

        Args:
            img: image data
            mode: load or write

        Returns:
            rotated and/or inversed image
        """
        if mode not in ("load", "write"):
            raise AttributeError("mode must be 'load' or 'write'")

        if mode == "load":
            if self.rotate:
                img = np.rot90(img, k=1)
                self.print_if_info_is_true("rotates the image in 90 deg when loading")
            if self.inverse:
                img = img[::-1, :]
                self.print_if_info_is_true(
                    "inverse the image along the 0th-axis when loading"
                )
        elif mode == "write":
            if self.inverse:
                img = img[::-1, :]
                self.print_if_info_is_true(
                    "inverse the image along the 0th-axis when writing"
                )
            if self.rotate:
                img = np.rot90(img, k=-1)
                self.print_if_info_is_true("rotates the image in 90 deg when writing")

        return img

    def load_fitsset(self):
        """Load fitsset and make imcube.

        Returns:
            imcube
        """
        imcube = []
        for data in tqdm.tqdm(self.rawpath):
            im, _ = load_fits_data_header(str(data))
            im = self.detector_handling(im, mode="load")
            imcube.append(im)
        return np.array(imcube)

    def immedian(self, extension=None):
        """take image median.

        Return:
           median image
        """
        e = self.extension
        if not extension is None:
            self.extension = extension
            self.print_if_info_is_true(
                f"Performing median combine for images with extension '{self.extension}'."
                )
        else:
            self.print_if_info_is_true("Performing median combine for raw images.")

        imall = []
        for path in tqdm.tqdm(self.path(string=True, check=True)):
            img, _ = load_fits_data_header(path)
            img = self.detector_handling(img, mode="load")
            imall.append(img)
        imall = np.array(imall)
        # np.nansum(corrected_im_all,axis=0)#
        median_image = np.nanmedian(imall, axis=0)
        self.extension = e
        return median_image

    def extract_trace_info(self, trace_path=None):
        """extract parameters for trace

        Args:
            trace_path: path to trace file

        Returns:
            data and parameters for aperture trace
        """
        from pyird.io.iraf_trace import read_trace_file

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

    ############################################################################################

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

        self.print_if_info_is_true("[STEP] Performing clean_pattern: output extension=" + extout)

        if trace_mask is None:
            trace_mask = self.trace.mask()

        if extin is None:
            extin = self.extension

        extin_noexist, extout_noexist = self.check_existence(extin, extout)
        for i, fitsid in enumerate(tqdm.tqdm(extin_noexist)):
            filen = self.rawdir / extin_noexist[i]
            im, header = load_fits_data_header(filen)
            im = self.detector_handling(im, mode="load")
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

            corrected_im = self.detector_handling(corrected_im, mode="write")
            write_fits_data_header(
                self.anadir / extout_noexist[i], header, corrected_im
            )

        self.fitsdir = self.anadir
        self.extension = extout

    def flatten(
        self,
        trace_path=None,
        extout="_fl",
        extin=None,
        hotpix_mask=None,
        width=None,
        check=False,
        master_path=None,
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

        self.print_if_info_is_true("[STEP] Performing flatten...")

        y0, xmin, xmax, coeff, mmf, width = self.extract_trace_info(trace_path)

        if extin is None:
            extin = self.extension

        extout = extout + "_" + mmf
        if not hotpix_mask is None:
            extout = "_hp_" + mmf

        if self.imcomb:
            out_path = self.anadir / ("%s_%s_%s.fits" % (self.streamid, self.band, mmf))
            extout_noexist = [out_path] * (not os.path.exists(out_path))
        else:
            extin_noexist, extout_noexist = self.check_existence(extin, extout)
        if check:
            extin_noexist = self.extpath(extin, string=False, check=False)
            extout_noexist = self.extpath(extout, string=False, check=False)

        for i, extout_i in tqdm.tqdm(enumerate(extout_noexist)):
            if self.imcomb:
                im = self.immedian()  # median combine
                filen = self.path()[0]  # header of the first file
                _, header = load_fits_data_header(filen)
            else:
                filen = self.anadir / extin_noexist[i]
                im, header = load_fits_data_header(filen)
                im = self.detector_handling(im, mode="load")
            rawspec, pixcoord, rotim, tl, iys_plot, iye_plot = flatten(
                im, trace_legendre, y0, xmin, xmax, coeff, inst=self.inst, width=width
            )
            rsd = multiorder_to_rsd(rawspec, pixcoord)
            if not hotpix_mask is None:
                save_path = self.anadir / ("hotpix_%s_%s.fits" % (self.band, mmf))
                rsd = apply_hotpixel_mask(
                    hotpix_mask, rsd, y0, xmin, xmax, coeff, save_path=save_path
                )
            if not check:
                write_fits_data_header(self.anadir / extout_noexist[i], header, rsd)
            else:
                mask_shape = (2048, 2048)
                trace_mask = np.zeros(mask_shape)
                for i in range(len(y0)):
                    trace_mask[i, xmin[i] : xmax[i] + 1] = tl[i]
                if not master_path is None:
                    wav, _ = load_fits_data_header(master_path)
                else:
                    wav = []
                return rsd, wav, trace_mask, pixcoord, rotim, iys_plot, iye_plot

            self.print_if_info_is_true("Created " + str(extout_i))
        self.fitsdir = self.anadir
        self.extension = extout

    def apnormalize(
        self, rsd=None, hotpix_mask=None, ignore_orders=None, **kwargs_continuum
    ):
        """normalize 2D apertures by 1D functions

        Returns:
            dictionary of pandas DataFrames of extracted and normalized spectra in each pixel
        """
        from pyird.spec.continuum import ContinuumFit

        self.print_if_info_is_true("[STEP] Performing apnormalize...")

        y0, xmin, xmax, coeff, mmf, width = self.extract_trace_info()

        if rsd is None:
            flatfile = self.anadir / ("%s_%s_%s.fits" % (self.streamid, self.band, mmf))
            rsd, _ = load_fits_data_header(flatfile)

        if not hotpix_mask is None:
            save_path = self.anadir / ("hotpix_%s_%s.fits" % (self.band, mmf))
            rsd = apply_hotpixel_mask(
                hotpix_mask, rsd, y0, xmin, xmax, coeff, save_path=save_path
            )

        continuum_fit = ContinuumFit()
        update_attr(continuum_fit, **kwargs_continuum)
        df_continuum = continuum_fit.continuum_rsd(rsd, ignore_orders=ignore_orders)

        flat_median = self.immedian("_cp")
        if not hotpix_mask is None:
            flat_median[hotpix_mask] = np.nan
        df_onepix = flatten(
            flat_median,
            trace_legendre,
            y0,
            xmin,
            xmax,
            coeff,
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
        from pyird.spec.rsdmat import rsd_order_medfilt
        from pyird.image.oned_extract import sum_weighted_apertures

        self.print_if_info_is_true("[STEP] Performing apext_flatfield...")

        y0, xmin, xmax, coeff, mmf, width = self.extract_trace_info()

        if extin is None:
            extin = self.extension

        if hotpix_mask is None:
            extout = extout + "_" + mmf
        else:
            extout = extout + "hp_" + mmf

        if self.imcomb:
            out_path = self.anadir / (
                "n%s_%s_%s.fits" % (self.streamid, self.band, mmf)
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
                im, header = load_fits_data_header(filen)
                im = self.detector_handling(im, mode="load")
            df_sum_wap = sum_weighted_apertures(
                im, df_flatn, y0, xmin, xmax, coeff, width, self.inst
            )
            rsd = df_sum_wap.values.T.astype(float)
            if not hotpix_mask is None:
                save_path = self.anadir / ("hotpix_%s_%s.fits" % (self.band, mmf))
                rsd = apply_hotpixel_mask(
                    hotpix_mask, rsd, y0, xmin, xmax, coeff, save_path=save_path
                )
            if self.imcomb:
                rsd = rsd_order_medfilt(rsd)
            # save as fits file without applying rotation/inverse
            write_fits_data_header(self.anadir / extout_noexist[i], header, rsd)

    def calibrate_wavelength(
        self,
        trace_path=None,
        channelfile_path=None,
        ign_ord=[],
        maxiter=30,
        std_threshold=0.001,
        npix=2048,
        width=None,
        force_rotate=False,
    ):
        """wavelength calibration usgin Th-Ar.

        Args:
            trace_path: path to the trace file
            channelfile_path: path to the channel file
            ign_ord: orders to be ignored
            maxiter: maximum number of iterations
            std_threshold: When the std of fitting residuals reaches this value, the iteration is terminated.
            npix: number of pixels
            width: list of aperture widths ([width_start,width_end])
            force_rotate: forces rotating the detector, when the number of the order is not standard values (i.e. 21 or 51)


        """
        from pyird.spec.wavcal import wavcal_thar

        self.print_if_info_is_true("[STEP] Performing calibrate_wavelength...")

        median_image = self.immedian()

        y0, xmin, xmax, coeff, mmf, width = self.extract_trace_info(trace_path)

        master_path = self.anadir / ("thar_%s_%s.fits" % (self.band, mmf))
        if not os.path.exists(master_path):
            filen = self.path()[0]  # header of the first file
            _, header = load_fits_data_header(filen)
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
            w = np.ones(rsd.shape)  ##weights for identifying thar
            wavsol, data = wavcal_thar(
                rsd.T,
                w,
                maxiter=maxiter,
                std_threshold=std_threshold,
                channelfile_path=channelfile_path,
                ign_ord=ign_ord,
            )
            # np.save('thar_%s_%s_final.npy'%(self.band,mmf),data)
            wavsol_2d = wavsol.reshape((npix, nord))
            write_fits_data_header(master_path, header, wavsol_2d)
            self.print_if_info_is_true(
                "Created a new file of the ThAr spectrum with a wavelength solution: " + str(master_path)
                )

    @rename_12_20({"cutrow":"search_start_row", "nap":"num_aperture"})
    def aptrace(self, search_start_row=1000, num_aperture=42, ign_ord=[]):
        """extract aperture of trace from a median image of current fitsset

        Args:
            search_start_row: starting row number to search apertures
            num_aperture: number of apertures, num_aperture = 42 ## 42 for H band, 102 for YJ band

        Returns:
            TraceAperture instance

        """
        from pyird.image.aptrace import aptrace
        from pyird.utils.aperture import TraceAperture

        self.print_if_info_is_true("[STEP] Performing aptrace...")

        inst = self.inst
        band = self.band

        flatmedian = self.immedian()
        if band=='h':
            flatmedian = flatmedian[::-1, ::-1]

        if self.detector_artifact:
            for i in range(0, 16):
                flatmedian[63 + i * 128 : 63 + i * 128 + 1, :] = flatmedian[
                    63 + i * 128 - 1 : 63 + i * 128, :
                ]
                flatmedian[63 + i * 128 + 1 : 63 + i * 128 + 2, :] = flatmedian[
                    63 + i * 128 + 2 : 63 + i * 128 + 3, :
                ]

        y0, xmin, xmax, coeff = aptrace(flatmedian, search_start_row, num_aperture, ign_ord)

        return TraceAperture(trace_legendre, y0, xmin, xmax, coeff, inst)

    def dispcor(self, extin="_fl", prefix="w", master_path=None, blaze=True):
        """dispersion correct and resample spectra

        Args:
            extin: extension of input files
            prefix: prefix for output files
            master_path: path to the directory containing the master ThAr file
            blaze: set True for blaze function

        """
        from pyird.utils.datset import DatSet

        self.print_if_info_is_true("[STEP] Performing dispcor...")

        if master_path == None:
            master_path = self.anadir.joinpath("..", "thar").resolve() / (
                "thar_%s_%s.fits" % (self.band, self.trace.mmf)
            )
        elif not str(master_path).endswith(".fits"):
            master_path = master_path / (
                "thar_%s_%s.fits" % (self.band, self.trace.mmf)
            )
        self.print_if_info_is_true("Allocate wavelengths based on the ThAr file: " + str(master_path))

        extin = extin + "_" + self.trace.mmf

        reference, _ = load_fits_data_header(master_path)

        datset = DatSet(dir=self.anadir)
        datset.prefix = prefix
        if self.imcomb:
            datset.extension = "_%s_%s" % (self.band, self.trace.mmf)
            if blaze:
                inputs = [
                    self.anadir
                    / ("n%s_%s_%s.fits" % (self.streamid, self.band, self.trace.mmf))
                ]
                datset.fitsid = ["blaze"]
            else:  ## need for fringe removal?
                inputs = [
                    self.anadir
                    / ("%s_%s_%s.fits" % (self.streamid, self.band, self.trace.mmf))
                ]
                datset.fitsid = [self.streamid]
        else:
            inputs = self.extpath(extin, string=False, check=True)
            datset.extension = "_%s" % (self.trace.mmf)
            datset.fitsid = self.fitsid
        save_path = datset.path(string=True, check=False)

        for i, input in enumerate(inputs):
            spec_m2, _ = load_fits_data_header(input)

            wspec = self.write_df_spec_wav(spec_m2, reference, save_path[i])

            outpath = str(save_path[i]).split("/")[-1]
            print("Created 1D spectrum: %s" % (outpath))
            # plot
            show_wavcal_spectrum(
                wspec, title="Extracted spectrum: %s" % (outpath), alpha=0.5
            )

    def normalize1D(
        self,
        flatid="blaze",
        master_path=None,
        readout_noise_mode="default",
        skipLFC=None,
        **kwargs_normalize,
    ):
        """combine orders and normalize spectrum

        Args:
            flatid: streamid for flat data
            master_path: path to the directory containing the calibrated flat file
            readout_noise_mode: 'real' or 'default'. 'real' calculates the readout noise based on the observed LFC spectrum.

        Note:
            If readout_noise_mode='real', mmf1 data for Y/J band should be reduced at first.

        """
        from pyird.spec.normalize import SpectrumNormalizer

        self.print_if_info_is_true("[STEP] Performing normalize1D...")

        if skipLFC is not None:
            warn_msg = "The 'skipLFC' parameter is deprecated and will be removed in the next release. Please use 'readout_noise_mode' instead."
            warnings.warn(warn_msg, FutureWarning)
            if skipLFC:
                readout_noise_mode = "default"
            else:
                readout_noise_mode = "real"

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
            if readout_noise_mode == "real":
                LFC_path = [self.anadir / ("w%s_y_m1.dat" % (self.streamid,))]
        else:
            datset.fitsid = self.fitsid
            datset.extension = "_%s" % (self.trace.mmf)
            if readout_noise_mode == "real":
                LFC_path = []
                for id in self.fitsid:
                    if self.band == "h":
                        LFC_path.append(self.anadir / ("w%d_%s.dat" % (id - 1, "m1")))
                    elif self.band == "y":
                        LFC_path.append(self.anadir / ("w%d_%s.dat" % (id, "m1")))
        if readout_noise_mode == "default":
            LFC_path = [None] * len(self.fitsid)
        if readout_noise_mode not in ["default", "real"]:
            warnings.warn(
                f"'{readout_noise_mode}' is not a recognized mode. Using 'default' instead.",
                UserWarning,
            )
            # Set a default behavior or handle the unexpected input accordingly
            readout_noise_mode = "default"
        datset.prefix = "w"
        wfile = datset.path(string=True, check=False)
        datset.prefix = "nw"
        nwsave_path = datset.path(string=True, check=False)
        datset.prefix = "ncw"
        ncwsave_path = datset.path(string=True, check=False)

        for i in range(len(wfile)):
            spectrum_normalizer = SpectrumNormalizer(combfile=LFC_path[i])
            update_attr(spectrum_normalizer, **kwargs_normalize)
            df_continuum, df_interp = spectrum_normalizer.combine_normalize(
                wfile[i], flatfile, blaze=(flatid == "blaze")
            )
            df_continuum_save = df_continuum[
                ["wav", "order", "nflux", "sn_ratio", "uncertainty"]
            ]
            df_interp_save = df_interp[["wav", "nflux", "sn_ratio", "uncertainty"]]
            df_continuum_save.to_csv(nwsave_path[i], **self.tocsvargs)
            df_interp_save.to_csv(ncwsave_path[i], **self.tocsvargs)
            if not self.imcomb:
                print(
                    "Created normalized 1D spectrum: nw%d_%s.dat"
                    % (self.fitsid[i], self.trace.mmf)
                )
                print(
                    "Created normalized & order-combined 1D spectrum: ncw%d_%s.dat"
                    % (self.fitsid[i], self.trace.mmf)
                )
            # plot
            show_wavcal_spectrum(
                df_continuum_save,
                title="Normalized spectrum: nw%d_%s.dat"
                % (self.fitsid[i], self.trace.mmf),
                alpha=0.5,
            )
            show_wavcal_spectrum(
                df_interp_save,
                title="Normalized & Order combined spectrum: ncw%d_%s.dat"
                % (self.fitsid[i], self.trace.mmf),
                alpha=0.5,
            )


class Stream1D(DatSet):
    """Class for post-processing 1D spectra.
    """
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
        DatSet.__init__(self, rawdir, prefix=prefix, extension=extension)
        StreamCommon.__init__(self, streamid, rawdir, anadir, inst)

        self.readargs = {
            "header": None,
            "sep": "\s+",
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

    ############################################################################################

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
