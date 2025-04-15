from pyird.spec.continuum import ContinuumFit
from pyird.spec.uncertainty import FluxUncertainty   
from scipy import optimize
import numpy as np
import pandas as pd
from astropy.stats import sigma_clip

__all__ = ["SpectrumNormalizer"]

class SpectrumNormalizer(ContinuumFit, FluxUncertainty):
    """Spectrum Normalization Class"""

    def __init__(self, 
                combfile = None,
                interp_nonzero_ind = [2, 9], 
                nsigma_sigmaclip = [1, 2], 
                maxiter_sigmaclip = 3
                ):
        """initialize SpectrumNormalizer

        Args:
            interp_nonzero_ind: pixels to cut the zero flux region at both ends of the order
            nsigma_sigmaclip: sigma clipping threshold: tuple (low, high) or list
            maxiter_sigmaclip: the maximum number of iteration for sigma clipping

        """
        ContinuumFit.__init__(self)
        FluxUncertainty.__init__(self,combfile=combfile)

        self.interp_nonzero_ind = interp_nonzero_ind 
        self.nsigma_sigmaclip = nsigma_sigmaclip 
        self.maxiter_sigmaclip = maxiter_sigmaclip 
    
    def determine_scale_continuum(self, wdata, flat, standard_order):
        """determine the scaling factor of the blaze function

        Args:
            wdata: the wavelength calibrated target spectrum
            flat: the wavelength calibrated FLAT
            standard_order: Once an order number is set, the blaze functions are standardized based on that order

        Returns:
            scaling factor of the blaze function

        """
        sigma_lower, sigma_upper = self.nsigma_sigmaclip
        if standard_order is not None:
            wav_tmp, flux_tmp, flatflux_tmp, useind, continuum = self.continuum_oneord(wdata, flat, standard_order)
            args = (flux_tmp[useind], continuum, sigma_lower, sigma_upper, self.maxiter_sigmaclip)
            scale = optimize.brent(self.objective_function, args=args)
        else:
            scale = 1
        return scale

    def objective_function(self, scale, flux, continuum, sigma_lower, sigma_upper, maxiter_sigmaclip):
        """objective function to optimize the scaling factor of the blaze function

        Args:
            scale: scaling factor of continuum
            flux: flux
            continuum: continuum
            sigma_lower: sigma clipping threshold (lower)
            sigma_upper: sigma clipping threshold (upper)
            maxiter_sigmaclip: the maximum number of iterations for sigma clipping

        """
        res = flux - continuum * scale
        res_clipped = sigma_clip(res, sigma_lower=sigma_lower, sigma_upper=sigma_upper, maxiters=maxiter_sigmaclip)
        abs_sum = np.abs(np.sum(res_clipped))
        return abs_sum

    def combine_normalize(self, wfile, flatfile, blaze=True):
        """read .dat file and make 1D normalized spectrum

        Args:
            wfile: path to the wavelength calibrated 1D target spectrum
            flatfile: path to the wavelength calibrated 1D FLAT or blaze file
            blaze: True when self.apext_flatfield() is used. if False, blaze function will be created from 1D FLAT

        Returns:
            pandas.DataFrame of 1D normalized spectrum
        """

        read_args = {'header': None, 'sep': r'\s+', 'names':['wav','order','flux']}
        wdata = pd.read_csv(wfile, **read_args)
        flat = pd.read_csv(flatfile, **read_args)

        if blaze:
            df_continuum = self.blaze_to_df(wdata,flat)
        else:
            df_continuum = self.make_blaze(wdata,flat)

        df_interp = self.normalize(df_continuum)
        return df_continuum, df_interp

    def blaze_to_df(self, wdata, blaze):
        """divide target flux by the blaze function

        Args:
            wdata: the wavelength calibrated 1D target spectrum
            blaze: blaze function

        Returns:
            normalized spectrum of each order
        """
        blaze = blaze.assign(blaze_continuum=0)
        blaze = blaze.rename(columns={'flux':'continuum'})
        df_continuum = pd.merge(wdata,blaze,on=['wav','order'])
        df_continuum['nflux'] = df_continuum['flux']/df_continuum['continuum']

        df_continuum = self.calc_uncertainty(df_continuum)
        
        return df_continuum

    def make_blaze(self, wdata, flat, standard_order=None):
        """extract blaze function for target based on FLAT

        Args:
            wdata: the wavelength calibrated target spectrum
            flat: the wavelength calibrated FLAT
            standard_order: Once an order number is set, the blaze functions are standardized based on that order

        Return:
            the blaze function created by scaling FLAT by a constant
        """
        scale = self.determine_scale_continuum(wdata, flat, standard_order)

        orders = wdata['order'].unique()
        df_continuum = pd.DataFrame([])
        for order in orders:
            wav_tmp, flux_tmp, flatflux_tmp, useind, continuum = self.continuum_oneord(wdata, flat, order)
            order_tmp = np.ones(len(continuum)) * order
            data = np.array([wav_tmp[useind], order_tmp, flux_tmp[useind], flatflux_tmp[useind], continuum*scale]).T
            df_tmp = pd.DataFrame(data, columns=['wav','order','flux','flat','continuum'])
            df_continuum = pd.concat([df_continuum, df_tmp],ignore_index=True)

        df_continuum['nflux'] = df_continuum['flux']/df_continuum['continuum']

        df_continuum = self.calc_uncertainty(df_continuum)

        return df_continuum

    def normalize(self, df_continuum):
        """normalize flux after combining all orders

        Args:
            df_continuum: pandas.DataFrame that contain the blaze function

        Returns:
            pandas.DataFrame of 1D normalized spectrum
        """

        orders = df_continuum['order'].unique()

        df_interp = pd.DataFrame([],columns=['wav','flux','continuum','sn_ratio','tmp_uncertainty'])
        for order in orders:
            df_former, df_latter = self.define_former_and_latter(df_continuum, order, max(orders))

            if order == max(orders):
                df_interp = self.concat_nonoverlap_region_last_order(df_former, df_interp)
                continue
            else:
                df_interp = self.concat_nonoverlap_region(df_former, df_latter, df_interp, order, min(orders))

            df_interp = self.concat_overlap_region(df_former, df_latter, df_interp)

        df_interp = self.divide_by_continuum(df_interp)
        return df_interp

    def define_former_and_latter(self, df_continuum, order, max_order):
        """define data of the former order (order N) and the latter order (order N+1)

        Args:
            df_continuum: pandas.DataFrame that contain the blaze function
            order: order number (N)
            max_order: the maximum number of the orders

        Returns:
            pandas.DataFrames of the former order and the latter order

        """
        df_former = df_continuum[df_continuum['order'] == order]
        df_latter = df_continuum[df_continuum['order'] == (int(order) + 1)]

        df_former = self.trim_nonzero_flux(df_former)
        if order != max_order:
            df_latter = self.trim_nonzero_flux(df_latter)
            return df_former, df_latter
        else:
            return df_former, None

    def trim_nonzero_flux(self, df):
        """cut pixels in the zero flux region at both ends of the order

        Args:
            df: pandas.DataFrame that contains 'flux' column

        Returns:
            pandas.DataFrame which is cut by +/- interp_nonzero_ind pixels 

        """
        interp_nonzero_beginning, interp_nonzero_ending = self.interp_nonzero_ind
        ind_nonzero_flux = np.where(df.flux != 0)[0]
        return df[min(ind_nonzero_flux) + interp_nonzero_beginning : max(ind_nonzero_flux) - interp_nonzero_ending]

    def concat_nonoverlap_region_last_order(self, df_former, df_interp):
        """concatenate the data of the nonoverlap region of df_former and df_interp for the last order

        Args:
            df_former: pandas.DataFrame of the former order
            df_interp: pandas.DataFrame of the interpolated data

        Returns:
            concateneted DataFrame

        """
        add_ind = (df_former['wav'] > df_interp['wav'].values[-1]) # non-overlapping region
        df_interp = pd.concat([df_interp, df_former[add_ind][df_interp.columns]])
        return df_interp

    def concat_nonoverlap_region(self, df_former, df_latter, df_interp, order, min_order):
        """concatenate the data of the nonoverlap region of df_former and df_latter

        Args:
            df_former: pandas.DataFrame of the former order
            df_latter: pandas.DataFrame of the latter order
            df_interp: pandas.DataFrame of the interpolated data
            order: number of the using order
            min_order: number of the initial order

        Returns:
            concateneted DataFrame

        """
        wav_former = df_former['wav'].values
        wav_latter = df_latter['wav'].values
        wav_tail_tip = wav_former[-1]
        wav_head_tip = wav_latter[0]

        ## non-overlapping region
        if order==min_order:
            add_ind = wav_former < wav_head_tip
        else:
            add_ind = (wav_former > df_interp['wav'].values[-1]) & (wav_former < wav_head_tip)
        
        df_interp = pd.concat([df_interp, df_former[add_ind][df_interp.columns]])
        return df_interp

    def concat_overlap_region(self, df_former, df_latter, df_interp):
        """concatenate the data of the orverlapping region of df_former and df_latter 

        Args:
            df_former: pandas.DataFrame of the former order
            df_latter: pandas.DataFrame of the latter order
            df_interp: pandas.DataFrame of the interpolated data

        Returns:
            concatenated DataFrame
            
        """
        wav_former = df_former['wav'].values
        wav_latter = df_latter['wav'].values
        wav_tail_tip = wav_former[-1]
        wav_head_tip = wav_latter[0]

        ## overlapping region
        ## head and tail
        df_head = df_latter[wav_latter <= wav_tail_tip]
        wav_head_tip_in_former = wav_former[wav_former < wav_head_tip][-1]
        df_tail = df_former[wav_former >= wav_head_tip_in_former]

        ## interp and add
        wav_head = df_head['wav'].values
        wav_tail = df_tail['wav'].values
        flux_tail_interp = np.interp(wav_head, wav_tail, df_tail['flux'].values)
        continuum_tail_interp = np.interp(wav_head, wav_tail, df_tail['continuum'].values)
        #
        sn_ratio_interp, tmp_uncertainty = self.calc_uncertainty_overlap_region(df_head, df_tail)
        #

        flux_sum = df_head['flux'].values + flux_tail_interp
        continuum_sum = df_head['continuum'].values + continuum_tail_interp

        uncertainty_ratio_square_sum = np.sqrt(df_head['tmp_uncertainty'].values**2 + tmp_uncertainty**2)
        sn_ratio_sum = flux_sum / uncertainty_ratio_square_sum

        data = np.array([df_head['wav'].values, flux_sum, continuum_sum, sn_ratio_sum, uncertainty_ratio_square_sum]).T
        df_tmp = pd.DataFrame(data, columns=df_interp.columns)
        df_interp = pd.concat([df_interp, df_tmp])
        return df_interp

    def divide_by_continuum(self, df_interp):
        """normalizing by deviding by the continuum

        Args:
            df_interp: pandas.DataFrame that contains flux, continuum, and tmp_uncertainty

        Returns:
            add the normalized flux and the corresponding uncertainty to the DataFrame.
        """
        zeroind = df_interp['continuum']==0
        df_interp = df_interp.assign(nflux=np.nan,uncertainty=np.nan)
        df_interp.loc[~zeroind,'nflux'] = df_interp['flux'][~zeroind]/df_interp['continuum'][~zeroind]
        df_interp.loc[~zeroind,'uncertainty'] = df_interp['tmp_uncertainty'][~zeroind]/df_interp['continuum'][~zeroind]
        return df_interp
    
if __name__ == '__main__':
    wfile = '/Users/yuikasagi/IRD/PhDwork/pyird/data/20210317/target/w41511_m2.dat'
    flatfile = '/Users/yuikasagi/IRD/PhDwork/pyird/data/20210317/flat/wflat_m2.dat'
    df_interp = SpectrumNormalizer.combine_normalize(wfile, flatfile)
    print(df_interp)
